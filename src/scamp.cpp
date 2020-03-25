//This file is part of scamp, selective clustering annotated using modes of projections.
      
//scamp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//scamp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with scamp.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <Rcpp.h>
#include "scmp.h"
#include <unordered_map>

#ifndef _VERBOSE_NOISY_SCAMP_COMPILATION_
#define _VERBOSE_NOISY_SCAMP_COMPILATION_ false
#endif

#ifndef _NOISY2_DEBUG_USER_ANNOTATIONS_
#define _NOISY2_DEBUG_USER_ANNOTATIONS_ false
#endif

std::vector<std::vector<int>> scamp(const std::vector<std::vector<double>>& dataValsIn,
				    double dipT,
				    int clusterLB,
				    bool repeatsAllowed,
				    int maxSearchDepth,
				    long maxClusterNum,
				    bool debugScamp,
				    std::vector<std::string> clusterAnnotations,
				    long maxNumberOfGates,
				    bool randomSearch,
				    bool randomResidualSearch,
				    const std::vector<double>& finalAnnotationQs,
				    unsigned long numThreadsToUse,
				    bool useRestrictedValue,
				    const std::vector<std::vector<int>>& restrictedVals,
				    bool useForestValues,
				    Rcpp::List forestValues,
				    const double& maxAllowedTime,
				    double gaussianScaleParam,
				    unsigned long long& randomSeed,
				    double depthScoreThreshold,
				    unsigned long subSampleThreshold,
				    unsigned long subSampleSize,
				    unsigned long subSampleIterations)
{
  auto colNum = dataValsIn.size();

  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Adding noise to data matrix." << std::endl;
  }
  //To begin, add shape-preserving Uniform and Gaussian noise to input data.
  std::vector<std::vector<double>> dataVals = addNoiseToDataMatrix(dataValsIn,numThreadsToUse,gaussianScaleParam,randomSeed);

  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Task complete. Noise added to data matrix." << std::endl;
  }

  //Next, we search for candidate clusters.
  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Starting search for candidate clusters." << std::endl;
  }

  searchResults candClusters = findCandidateClusters(dataVals,
						     restrictedVals,
						     dipT,
						     clusterLB,
						     repeatsAllowed,
						     maxSearchDepth,
						     maxClusterNum,
						     maxNumberOfGates,
						     useRestrictedValue,
						     randomSearch,
						     maxAllowedTime,
						     debugScamp,
						     false,
						     numThreadsToUse,
						     randomSeed,
						     subSampleThreshold,
						     subSampleSize,
						     subSampleIterations);


  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Task complete. Have completed search for candidate clusters." << std::endl;
  }

  if (candClusters.abortIteration) {
    //if the search for candidate clusters has run too long, signal
    //this to cppNoisyScamp with a special vector.
    //NOTE: a successful annotatedClustering has minimum value 0,
    //so cppNoisy scamp can check for this condition.
    std::cout << "Initial search abort." << std::endl;
    std::vector<int> failVec = {-1,-1};
    std::vector<std::vector<int>> failContainer(1,failVec);
    return failContainer;
  }

  //record gate locations found by the candidate cluster search
  std::vector<gateInfo> gatePlacements = candClusters.gateLocations;
  
  //Now we have found the requisitne number of candidate clusters.
  //An initial clustering is determined by selecting among them to parititon the data set.
  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Starting initial clustering selection." << std::endl;
  }
  std::vector<int> selectedClustering = selectCandidateClusters((candClusters.candidates),dataVals,debugScamp);
  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Task complete. An initial clustering of the data matrix has been determined." << std::endl;
  }

  //once an initial clustering has been selected, examine those observations still unclassified.
  //these observations have cluster number 0.
  //we will test these residual observations for unimodality (according to the dip test).
  //if they appear unimodal, they will be lumped together. Otherwise, the procedure will be repeated
  //on the subset of unclassified observations. This will continue until a final lumping is possible.

  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Starting to cluster residual observations." << std::endl;
  }

  int cMaxClusNum = *std::max_element(selectedClustering.begin(),selectedClustering.end());
  auto numObs = selectedClustering.size();
  long long zeroCount = 0;
  std::vector<bool> zeroActive(numObs,false);
  std::vector<long> zeroIndex;
  
  for (auto i = 0; i != numObs; ++i) {
    if (selectedClustering[i] == 0) {
      zeroCount += 1;
      zeroActive[i] = true;
      zeroIndex.push_back(i);
    }
  }

  //if the zeroCount contains more observations than an admissible number of cluster
  //we examine the observations with a class label of 0 for multimodality.
  bool stillViable;
  std::vector<std::vector<double>> dataValSubset, sortedSubset;
  std::vector<std::vector<int>> restrictedValSubset;
  searchResults residClusters;
  std::vector<int> residClustering;
  int tmpClusterNum;
  bool residualBelowMaxTime = true;
  auto lastSize = (dataVals[0]).size();
  auto currentSize = (dataVals[0]).size();
  while ((zeroCount > clusterLB) && residualBelowMaxTime){
    //begin by collecting the unclassified subset of observations.
    //by construction, unclassified equivalent to cluster label of 0.
    dataValSubset.clear();
    dataValSubset.resize(colNum);
    restrictedValSubset.clear();
    restrictedValSubset.resize(colNum);
    
    for (auto i = 0; i != numObs; ++i) {
      if (zeroActive[i]) {
	for (auto j = 0; j != colNum; ++j) {
	  (dataValSubset[j]).push_back((dataVals[j])[i]);
	  if (useRestrictedValue) {
	    (restrictedValSubset[j]).push_back((restrictedVals[j])[i]);
	  }
	}
      }
    }
    //create a sorted version of the subset.
    sortedSubset = dataValSubset;
    for (auto it = sortedSubset.begin(); it != sortedSubset.end(); ++it){
      std::sort((*it).begin(),(*it).end());
    }

    //check for pathological selection path
    //abort and try again if detected
    currentSize = (sortedSubset[0]).size();
    /*
    if ((lastSize - currentSize) < (2*clusterLB)) {
      if (debugScamp) {
	std::cout << "Pathological selection: step-size. Abort clustering and re-noise." << std::endl;
      }
      std::vector<int> failVec = {-1,-1};
      std::vector<std::vector<int>> failContainer(1,failVec);
      return failContainer;
    }
    */
    //update for next iteration
    lastSize=currentSize;
    
    if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
      std::cout << "Residual population consists of " << currentSize << " observations." << std::endl;
    }

    //using the sorted subset, check if the residual cluster exhibits multimodality.
    stillViable = false;
    for (auto i = 0; i != colNum; ++i) {
      if (singleDip(sortedSubset[i]) < dipT) {
	stillViable = true;
	break;
      }
    }
    if (stillViable) {
      if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
	std::cout << "Begin residual search for candidate clusters." << std::endl;
      }
      //residual multimodality is found so we repeat the procedure on the observations in the residual cluster.
      residClusters = findCandidateClusters(dataValSubset,restrictedValSubset,dipT,clusterLB,
					    repeatsAllowed,maxSearchDepth,maxClusterNum,maxNumberOfGates,
					    useRestrictedValue,randomSearch,maxAllowedTime,debugScamp,false,
					    numThreadsToUse,randomSeed,
					    subSampleThreshold,subSampleSize,subSampleIterations);
      if (residClusters.abortIteration) {
	std::cout << "Residual search abort." << std::endl;
	std::vector<int> failVec = {-1,-1};
	std::vector<std::vector<int>> failContainer(1,failVec);
	return failContainer;
      }
      
      if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
	std::cout << "Search complete." << std::endl;
	std::cout << "Begin residual selection of candidate clusters." << std::endl;
      }

      //select among the residual candidates
      residClustering = selectCandidateClusters((residClusters.candidates),dataVals,debugScamp);
      
      if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
	std::cout << "Selection complete." << std::endl;
	std::cout << "Update main clustering." << std::endl;
      }

      //use the residual clustering to update the initial cluster assignment.
      for (auto k = 0; k != residClustering.size(); ++k) {
	tmpClusterNum = residClustering[k];
	if (tmpClusterNum > 0) {
	  selectedClustering[zeroIndex[k]] = tmpClusterNum + cMaxClusNum;
	  zeroActive[zeroIndex[k]] = false;
	}
      }
      
      //finally update cluster num and collect new zero indices/zero count.
      cMaxClusNum = *std::max_element(selectedClustering.begin(),selectedClustering.end());
      zeroIndex.clear();
      zeroCount = 0;
      for (auto k = 0; k != numObs; ++k) {
	if (selectedClustering[k] == 0) {
	  zeroCount += 1;
	  zeroIndex.push_back(k);
	}
      }
      if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
	std::cout << "Update complete." << std::endl;
      }
    }
    else {
      //multimodality is not detected. terminate procedure.
      zeroCount = 0;
    }
  }

  
  //at this point, the data has been clustered. we now proceed to annotate numerical clusters with
  //labels that attempt to describe their key features. Since each candidate cluster obeys the definition
  //of an alpha-m cluster at this point, we will label clusters by comparing their user-selected quantiles to the median
  //gate location by channel in the dataset. 
  
  // note: the annotationGates store *ONLY* the phenotypically relevant gate data.
  // the column min & max are **NOT** stored, contrary to the return value of tsGates and tsModeEstimate

  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Starting to annotate clusters." << std::endl;
  }

  std::vector<double> tmpBnd;
  std::vector<std::vector<double>> annotationGates;
  annotationGates.resize(colNum);
  std::vector<bool> annCols;
  if ((useForestValues) || (randomSearch)) {
    if (useForestValues) {
      for (auto i = 0; i != colNum; ++i){
	annCols.push_back(true);
      }
    }
    else {
      annCols = determineAnnotationColumns(dataVals,colNum,dipT,restrictedVals,useRestrictedValue);
    }
    for (auto i = 0; i != colNum; ++i){
      //if we are using forestValues, use them instead of sample-specific boundaries.
      if (useForestValues) {
	tmpBnd = Rcpp::as<std::vector<double>>(forestValues[i]);
	for (auto j = 0; j != tmpBnd.size(); ++j) {
	  (annotationGates[i]).push_back(tmpBnd[j]);
	}
      }
      else if (annCols[i] == true) {
	tmpBnd = determineAnnotationBoundaries(gatePlacements,i,debugScamp);
	if (tmpBnd.size()) {
	  for (auto j = 0; j != tmpBnd.size(); ++j) {
	    (annotationGates[i]).push_back(tmpBnd[j]);
	  }
	}
      }
    }
  }
  
  else {
    //the candidate cluster search was exhaustive.
    columnSummary defaultSummary, currentSummary;
    std::vector<columnSummary> dataSummary(colNum,defaultSummary);
    int rootsInForest = 1; 
    
    //parse the gate placements for depth calculations.
    gateInfo currentGateInfo;
    int gateColumnNum;
    int currentGateDepth;
    unsigned long currentNumGates;
    std::vector<double> gateLocs;
    for (auto gatePlaceNum = 0; gatePlaceNum != gatePlacements.size(); ++gatePlaceNum) {
      currentGateInfo = gatePlacements[gatePlaceNum];
      gateColumnNum = currentGateInfo.colNumber;
      gateLocs = currentGateInfo.gates;
      currentNumGates = ((currentGateInfo.numGates)-2);
      currentGateDepth = currentGateInfo.gateDepth;
      if (currentNumGates > maxNumberOfGates) {
	continue; //user supplies the upper bound to ignore.
      }
      else {
	dataSummary[gateColumnNum].depthMap[currentNumGates].push_back(currentGateDepth);
	for (int gateLocNum = 1; gateLocNum != (gateLocs.size()-1); ++gateLocNum) 
	  ((dataSummary[gateColumnNum].gateMaps[currentNumGates])[gateLocNum]).push_back(gateLocs[gateLocNum]);
      }
    }

    if (debugScamp) {
      //check the parse.
      std::unordered_map<int, std::vector<double>> tmpMap;
      for (auto i = 0; i != dataSummary.size(); ++i) {
	currentSummary = dataSummary[i];
	std::cout << "Column " << (i+1) << " summary." << std::endl;
	for (auto p : currentSummary.gateMaps) {
	  //std::cout << p.first  << "-gates " << " contribution " << ((p.second)[0]).size() << std::endl;
	  std::cout << p.first << ":-gates located at ";
	  tmpMap = p.second;
	  for (auto v : tmpMap){
	    for (auto q : v.second){
	      std::cout << q << " ";
	    }
	    std::cout << std::endl;
	  }
	}
	std::cout << std::endl;
	std::cout << "Gate depths." << std::endl;
	for (auto p : currentSummary.depthMap) {
	  std::cout << p.first << "-gates at depth: ";;
	  for (auto v : p.second) {
	    std::cout << v << " ";
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;
      }
    }

    //add exhaustive labeling
    exhaustiveAnnotations exAnnVs= determineAnnotationsEx(dataVals,colNum,dipT,
							  restrictedVals,useRestrictedValue,
							  dataSummary,rootsInForest,
							  depthScoreThreshold);

    
    //copy out cut points for annotation
    annCols = exAnnVs.columnIsAnnotated;
    for (auto gn = 0; gn != exAnnVs.cutPointLocations.size(); ++gn) {
      if ((exAnnVs.cutPointLocations)[gn].size()) {
	if (debugScamp) 
	  std::cout << "Cut points for column " << (gn + 1) << " follow." << std::endl;
	for (auto jn = 0; jn != (exAnnVs.cutPointLocations)[gn].size(); ++jn) {
	  (annotationGates[gn]).push_back(((exAnnVs.cutPointLocations)[gn])[jn]);
	  if (debugScamp) 
	    std::cout << (((exAnnVs.cutPointLocations)[gn])[jn]) << " ";
	}
	if (debugScamp)
	  std::cout << std::endl; 
      }
    }
  }

  
  
  //finally, label the selectedClustering with annotations.
  std::vector<std::vector<int>> annotatedClustering;
  if (numThreadsToUse > 1) {
    //try passing the raw data after clustering, so that the labels match the observed values
    annotatedClustering = parallelAnnotateCluster(dataValsIn, //had passed dataVals 
						  selectedClustering,
						  annCols,
						  annotationGates,
						  clusterAnnotations,
						  finalAnnotationQs,
						  numThreadsToUse);
  }
  else {
    annotatedClustering = annotateCluster(dataValsIn, //had passed dataVals
					  selectedClustering,
					  annCols,
					  annotationGates,
					  clusterAnnotations,
					  finalAnnotationQs);
  }
  /*
  std::vector<std::vector<int>> annotatedClustering;
  annotatedClustering = annotateCluster(dataVals,
					selectedClustering,
					annCols,
					annotationGates,
					clusterAnnotations,
					finalAnnotationQs);
  */
  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Task complete. The clusters are labeled. Returning clustering." << std::endl;
  }
  return annotatedClustering;
}
