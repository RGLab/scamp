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
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <Rcpp.h>
#include "scmp.h"

bool anyIntersections(std::vector<long>& v1, std::vector<long>& v2) {
  bool anyIntersection = false;
  //v1, v2 are vectors of indices: natural numbers.
  //by construction, both are in sort order.
  //so, we can quickly determine if intersections are possible by checking the boundary points.
  long v1Min = v1[0];
  long v1Max = v1[(v1.size()-1)];
  long v2Min = v2[0];
  long v2Max = v2[(v2.size()-1)];

  if ((v1Min > v2Max) || (v1Max < v2Min))
    return anyIntersection;

  //some overlap occurs, so now must determine if any elements
  long value1, value2;

  for (auto i = 0; i != v1.size(); ++i) {
    value1 = v1[i];
    for (auto j = 0; j != v2.size(); ++j) {
      value2 = v2[j];
      if (value2 == value1) {
	anyIntersection = true;
	break;
      }
      //if value2 exceeds value1, so will all indices above it
      if (value2 > value1) {
	break;
      }
    }
    if (anyIntersection) {
      break;
    }
  }
  return anyIntersection;
}

double sumLmomentVector(const std::vector<double>& Lmomvec) {
  double totVal = 0.0;
  for (auto value : Lmomvec)
      totVal += value;
  return totVal; 
}


std::vector<double> scoreCandidateClusters(std::vector<std::vector<long>>& candClusters,
					   std::vector<std::vector<double>>& dataMat)
{
  int _LmomentNumber = 10; //fix the number of L-momnts
  int _LmomentTrimming = 2; //fix the symmetric trimming value
  auto numCandidates = candClusters.size();
  auto colNum = dataMat.size();
  std::vector<double> combinedScore(numCandidates,0.0);
  std::vector<double> sizeScore(numCandidates,0.0);
  double tmpSizeScore;
  std::vector<double> dipScore(numCandidates,0.0);
  std::vector<double> tmpDip(colNum,0.0);
  std::vector<double> lmomLinfScore2(numCandidates,0.0);
  std::vector<double> lmomLinfScore3(numCandidates,0.0);
  std::vector<double> lmomLinfScore4(numCandidates,0.0);

  std::vector<double> lmomTotScore2(numCandidates,0.0);
  std::vector<double> lmomTotScore3(numCandidates,0.0);
  std::vector<double> lmomTotScore4(numCandidates,0.0);
  
  std::vector<double> tmpLmom2(colNum,0.0);
  std::vector<double> tmpLmom3(colNum,0.0);
  std::vector<double> tmpLmom4(colNum,0.0);
  std::vector<double> tmpLmomVec;

  std::vector<long> currentCluster;
  std::vector<double> dataV;
  long maxSize = 0;
  long cSize = 0;
  //first, determine the maximum cluster size among the candidates
  for (auto i = 0; i != numCandidates; ++i) {
    cSize = (candClusters[i]).size();
    if (cSize > maxSize)
      maxSize = cSize;
  }

  //next iterate across the candiate clusters to score according to size, dip test p-value, and sum of standardized l-moments.
  for (auto i = 0; i != numCandidates; ++i) {
    currentCluster.clear();
    currentCluster = candClusters[i];
    //the size score of given candidate is its magnitude relative to the maximum cluster.
    //we assume larger clusters are more desirable than smaller clusters,
    //with the exception of the largest clusters: we assume the largest are found early in the search process.
    //to bias against them, we down-weight their score to 0.5
    tmpSizeScore = (currentCluster.size())/static_cast<double>(maxSize);
    if (tmpSizeScore < 0.95)
      sizeScore[i] = tmpSizeScore;
    else
      sizeScore[i] = 0.5; 

    for (auto j =0; j != colNum; ++j) {
      dataV.clear();
      for (auto k = 0; k != currentCluster.size(); ++k) {
	dataV.push_back(((dataMat[j])[(currentCluster[k])]));
      }
      //both singleDip and absSampleLmoments assume dataV is sorted.
      std::sort(dataV.begin(),dataV.end());
      tmpDip[j] = singleDip(dataV);
      tmpLmomVec = absSampleLmoments(dataV,_LmomentNumber,_LmomentTrimming);
      tmpLmom2[j] = tmpLmomVec[1];
      tmpLmom3[j] = tmpLmomVec[2];
      tmpLmom4[j] = tmpLmomVec[3];
    }
    //the column with minimum dip statistic is the dip score.
    //higher dip score == better shape.
    dipScore[i] = *std::min_element(tmpDip.begin(),tmpDip.end());
    //compute the column with the maximum L-moment values.
    lmomLinfScore2[i] = *std::max_element(tmpLmom2.begin(),tmpLmom2.end());
    lmomLinfScore3[i] = *std::max_element(tmpLmom3.begin(),tmpLmom3.end());
    lmomLinfScore4[i] = *std::max_element(tmpLmom4.begin(),tmpLmom4.end());

    //compute the sum of L-moment values across all columns
    lmomTotScore2[i] = sumLmomentVector(tmpLmom2);
    lmomTotScore3[i] = sumLmomentVector(tmpLmom3);
    lmomTotScore4[i] = sumLmomentVector(tmpLmom4);
  }

  //normalize the lmoment score to lie between 0 and 1, then flip so close to 1 is desirable
  double maxLmomLinfScore2 = *std::max_element(lmomLinfScore2.begin(),lmomLinfScore2.end()); 
  double maxLmomLinfScore3 = *std::max_element(lmomLinfScore3.begin(),lmomLinfScore3.end()); 
  double maxLmomLinfScore4 = *std::max_element(lmomLinfScore4.begin(),lmomLinfScore4.end()); 

  double maxLmomTotScore2 = *std::max_element(lmomTotScore2.begin(),lmomTotScore2.end()); 
  double maxLmomTotScore3 = *std::max_element(lmomTotScore3.begin(),lmomTotScore3.end()); 
  double maxLmomTotScore4 = *std::max_element(lmomTotScore4.begin(),lmomTotScore4.end()); 

  for (auto v = 0; v != lmomLinfScore2.size(); ++v) {
    lmomLinfScore2[v] = (1.0 - (lmomLinfScore2[v]/maxLmomLinfScore2)); // worst-case measure of L-moment "variance" standardized into (0,1) and inverted -- 1 most desirable, 0 least
    lmomLinfScore3[v] = (1.0 - (lmomLinfScore3[v]/maxLmomLinfScore3)); // worst-case measure of L-moment "skewness" standardized into (0,1) and inverted -- 1 most desirable, 0 least
    lmomLinfScore4[v] = (1.0 - (lmomLinfScore4[v]/maxLmomLinfScore4)); // worst-case measure of L-moment "kurtosis" standardized into (0,1) and inverted -- 1 most desirable, 0 least

    lmomTotScore2[v] = (1.0 - (lmomTotScore2[v]/maxLmomTotScore2)); // overall measure of L-moment "variance" for a candidate cluster
    lmomTotScore3[v] = (1.0 - (lmomTotScore3[v]/maxLmomTotScore3)); // overall measure of L-moment "skewness" for a candidate cluster
    lmomTotScore4[v] = (1.0 - (lmomTotScore4[v]/maxLmomTotScore4)); // overall measure of L-moment "kurtosis" for a candidate cluster
  }
  //finally compute a combined score for each cluster
  for (auto i = 0; i != numCandidates; ++i) {
    combinedScore[i] = ((dipScore[i] +
			 ((lmomLinfScore2[i] + lmomTotScore2[i])/(2.0)) +
			 ((lmomLinfScore3[i] + lmomTotScore3[i])/(2.0)) +
			 ((lmomLinfScore4[i] + lmomTotScore4[i])/(2.0)))); // * sizeScore[i]);

  }
  return combinedScore;
}

std::vector<int> selectCandidateClusters(std::vector<std::vector<long>>& candClusters,
					 std::vector<std::vector<double>>& dataMat,
					 const bool& debugSelection)
{
  //note that the number of candidate clusters in candClusters does not necessarily equal the length of a vector in dataMat.
  //furthermore, the convention of dataMat is each  vector<double> has the same size.
  if (debugSelection) {
    std::cout << "Starting to score candidate clusters." << std::endl;
  }
  std::vector<double> clusterScores = scoreCandidateClusters(candClusters,dataMat);
  if (debugSelection) {
    std::cout << "Candidate clusters have been scored." << std::endl;
  }

  auto candidateNumber = candClusters.size();
  std::vector<int> clustering(dataMat[0].size(),0);
  std::vector<bool> clusterActive(candidateNumber,true);
  double maxActiveScore, curScore;
  auto maxActiveIndex = candidateNumber;
  std::vector<long> selectedCluster;
  std::vector<long> compCluster;
  int currentClusterNum = 1;
  bool stillClustering = true;
  auto numActive = candidateNumber;


  while (stillClustering) {
    if (debugSelection) {
      std::cout << "Selection iteration. numActive: " << numActive << std::endl;
    }

    maxActiveScore = (-1.0); //start negative to select among candidate clusters with score of zero.
    maxActiveIndex = 0;
    //iterate over candidate clusters
    for (auto i = 0; i != candidateNumber; ++i) {
      //subset to clusters that can still be selected (the active clusters)
      if ((clusterActive[i]) == true) { //note these are boolean values, which we emphasize with comparison to true.
	//if a candidate cluster can still be selected, check if it has the best score
	curScore = clusterScores[i];
	if (curScore > maxActiveScore) {
	  maxActiveScore = curScore;
	  maxActiveIndex = i;
	}
      }
    }

    if (debugSelection) {
      std::cout << "maxActiveScore: " << maxActiveScore << std::endl;
      std::cout << "maxActiveIndex: " << maxActiveIndex << std::endl;
      std::cout << "currentClusterNum: " << currentClusterNum << std::endl;
    }

    
    selectedCluster = candClusters[maxActiveIndex]; //we select the cluster associate with the max available score
    clusterActive[maxActiveIndex] = false; //can not select this cluster again. 

    //iterate over the indices; update the clustering with the current number.
    for (auto &ind : selectedCluster) {
      clustering[ind] = currentClusterNum;
    }
    //increment the cluster assignment
    currentClusterNum += 1;
    //iterate over the clusters again
    for (auto i = 0; i != candidateNumber; ++i) {
      //again subset to active clusters
      if ((clusterActive[i]) == true) {
	compCluster = candClusters[i];
	//if any indices match between the selected cluster and an active cluster
	//remove the active cluster from contention. This avoids multiple clustering of a single observation.
	if (anyIntersections(selectedCluster,compCluster)) {
	  clusterActive[i] = false;
	}
      }
    }

    //finally count the number of active clusters that remain
    numActive = 0;
    for (bool bv : clusterActive)
      if (bv == true)
	numActive += 1;
    //if no more clusters are active, we are done.
    if (numActive == 0) {
      stillClustering = false;      
    }
  }
  return clustering;
}
