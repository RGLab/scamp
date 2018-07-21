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
#include <sys/resource.h>
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
				    bool useFixedAnnotationBoundaries,
				    Rcpp::List fixedAnnBdrys,
				    const double& maxAllowedTime,
				    double gaussianScaleParam,
				    unsigned long long& randomSeed)
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

  searchResults candClusters = findCandidateClusters(dataVals,restrictedVals,dipT,clusterLB,repeatsAllowed,
						     maxSearchDepth,maxClusterNum,maxNumberOfGates,
						     useRestrictedValue,randomSearch,maxAllowedTime,debugScamp,false,
						     numThreadsToUse,randomSeed);

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
    if ((lastSize - currentSize) < (2*clusterLB)) {
      if (debugScamp) {
	std::cout << "Pathological selection: step-size. Abort clustering and re-noise." << std::endl;
      }
      std::vector<int> failVec = {-1,-1};
      std::vector<std::vector<int>> failContainer(1,failVec);
      return failContainer;
    }
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
					    numThreadsToUse,randomSeed);

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

  
  //the data has been clustered. we now proceed to annotate numerical clusters with
  //labels that attempt to describe their key features. Since each candidate cluster obeys the definition
  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Starting to annotate clusters." << std::endl;
  }

  std::vector<double> tmpBnd;
  std::vector<std::vector<double>> annotationGates;
  annotationGates.resize(colNum);
  std::vector<bool> annCols;

  if (useFixedAnnotationBoundaries) {
    for (auto i = 0; i != colNum; ++i){
      annCols.push_back(true);
    }
  }
  else {
    annCols = determineAnnotationColumns(dataVals,colNum,dipT,restrictedVals,useRestrictedValue);
  }
  for (auto i = 0; i != colNum; ++i){
    //if we are using fixedAnnBdrys, use them instead of sample-specific boundaries.
    if (useFixedAnnotationBoundaries) {
      tmpBnd = Rcpp::as<std::vector<double>>(fixedAnnBdrys[i]);
      if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
	std::cout << "Number of boundary values for column " << (i+1) << ": " << tmpBnd.size() << std::endl;
      }
      for (auto j = 0; j != tmpBnd.size(); ++j) {
	if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
	  std::cout << "Boundary value: " << tmpBnd[j] << std::endl;
	}
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

  

  
  
  //finally, label the selectedClustering with annotations.
  std::vector<std::vector<int>> annotatedClustering;
  if (numThreadsToUse > 1) {
    annotatedClustering = parallelAnnotateCluster(dataVals,
						  selectedClustering,
						  annCols,
						  annotationGates,
						  clusterAnnotations,
						  finalAnnotationQs,
						  numThreadsToUse);
  }
  else {
    annotatedClustering = annotateCluster(dataVals,
					  selectedClustering,
					  annCols,
					  annotationGates,
					  clusterAnnotations,
					  finalAnnotationQs);
  }
  if (_VERBOSE_NOISY_SCAMP_COMPILATION_ || debugScamp) {
    std::cout << "Task complete. The clusters are labeled. Returning clustering." << std::endl;
  }
  return annotatedClustering;
}
