//This file is part of scamp, selective clustering using modes of projects.
      
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

std::vector<std::vector<int>> annotateCluster(const std::vector<std::vector<double>>& dataValsIn,
					 const std::vector<int>& clusteringIn,
					 std::vector<bool>& annotatingCols,
					 std::vector<std::vector<double>>& annotationBoundaries,
					 std::vector<std::string>& annotationsIn,
					 const std::vector<double>& annComparisonQuantiles){
  //first determine the unique cluster labels in the input clustering
  std::vector<int> uniqClusterLabels = clusteringIn;
  std::sort(uniqClusterLabels.begin(),uniqClusterLabels.end());
  auto uit = std::unique(uniqClusterLabels.begin(),uniqClusterLabels.end());
  uniqClusterLabels.resize(std::distance(uniqClusterLabels.begin(),uit));
  auto numberOfColumns = dataValsIn.size();
  auto numberOfObservations = clusteringIn.size();

  std::vector<double> quantileProbs = annComparisonQuantiles;
     
  std::vector<std::vector<double>> dvSubset;  
  int currentClusterLabel;
  
  std::vector<bool> inCluster(numberOfObservations,false);
  std::vector<unsigned long> clusterLookup;
  std::vector<double> qnt = {0.0,0.0,0.0};
  std::vector<double> clusterLowerQs, clusterMedians, clusterUpperQs;
  double cqLower,cq50,cqUpper,cDataValue,cAnnotationValue;
  std::vector<int> newLabel;
  std::vector<double> colAnnBdry;
  int colAnnScore = 0;
  int numAnnPts;
  std::vector<int> labelTemplate((2*numberOfColumns),0);
  std::vector<std::vector<int>> annotatedClusters(numberOfObservations,labelTemplate);

  for (auto clusterNum = 0; clusterNum != uniqClusterLabels.size(); ++clusterNum) {
    //inCluster flags if an observation of dataValsIn belongs to the cluster. reset at each loop iteration.
    std::fill(inCluster.begin(),inCluster.end(),false); 
    currentClusterLabel = uniqClusterLabels[clusterNum];
    
    //pick out observations in the current cluster.
    //clusteringIn.size() == dataValsIn[0].size().
    //that is, the clustering equals the number of rows in the initial data set.
    for (auto currentObs = 0; currentObs != clusteringIn.size(); ++currentObs) { 
      if (clusteringIn[currentObs] == currentClusterLabel)  {
	inCluster[currentObs] = true;
      }
    }
    
    //using the subset flags, get the correspond data from dataValsIn
    dvSubset.clear();
    dvSubset.resize(numberOfColumns);
    clusterLookup.clear();
    for (auto currentObsIndex = 0; currentObsIndex != numberOfObservations; ++currentObsIndex)  {
      if (inCluster[currentObsIndex]) {
	clusterLookup.push_back(currentObsIndex);
	for (auto currentColumn = 0; currentColumn != numberOfColumns; ++currentColumn) {
	  (dvSubset[currentColumn]).push_back((dataValsIn[currentColumn])[currentObsIndex]);
	}
      }
    }
    //compute cluster quantiles for all columns in the dataSet..
    clusterLowerQs.clear();
    clusterMedians.clear();
    clusterUpperQs.clear();
    
    for (auto currentColumn = 0; currentColumn != numberOfColumns; ++currentColumn) {
      qnt = rQuantile(dvSubset[currentColumn],quantileProbs);
      clusterLowerQs.push_back(qnt[0]);
      clusterMedians.push_back(qnt[1]);
      clusterUpperQs.push_back(qnt[2]);
    }
    
    for (auto clusterObsIndex = 0; clusterObsIndex != clusterLookup.size(); ++clusterObsIndex) {
      newLabel.clear();
      //determine new cluster label based off median of cluster in annotating columns.
      for (auto currentColNum = 0; currentColNum != numberOfColumns; ++currentColNum) {
	if (annotatingCols[currentColNum]) {
	  cqLower = clusterLowerQs[currentColNum];
	  cq50 = clusterMedians[currentColNum];
	  cqUpper = clusterUpperQs[currentColNum];
	  
	  colAnnBdry = annotationBoundaries[currentColNum];
	  numAnnPts = colAnnBdry.size();
	  colAnnScore = 0;
	  if (numAnnPts == 1) {
	    //if there is only one annotation boundary, 
	    //the decision to split the cluster follows by compare the lower and upper cluster quantiles to
	    //the the columns annotation boundary. 
	    cDataValue = (dataValsIn[currentColNum])[(clusterLookup[clusterObsIndex])];
	    cAnnotationValue = colAnnBdry[0]; //isCustomLabel == true implies this is the unique annotation value.
	    if (cqLower >= cAnnotationValue) {
	      //if the lower quantile exceeds the single annotation boundary, all obs in cluster == high
	      colAnnScore = 1;
	    }
	    else if (cqUpper <= cAnnotationValue) {
	      //if the upper quantile falls below the single annotation boundary, all obs in cluster == low
	      colAnnScore = 0;
	    }
	    else if (cDataValue >= cAnnotationValue) {
	      //if both previous conditions fail, split the cluster.
	      colAnnScore = 1;
	    }
	    else
	      colAnnScore = 0;
	  }
	  else {
	    //with multiple annotation boundaries, simply use the median.
	    for (auto j = 0; j != numAnnPts; ++j) {
	      if (cq50 >= colAnnBdry[j]){
		colAnnScore += 1;
	      }
	    }
	  }
	  //derive annotations for output.
	  newLabel.push_back(colAnnScore);
	  newLabel.push_back(numAnnPts);
	}
	else {
	  newLabel.push_back(0);
	  newLabel.push_back(0);
	}
      }
      //update annotations with the new label
      annotatedClusters[(clusterLookup[clusterObsIndex])] = newLabel;
    }
  }
  return annotatedClusters;
}
