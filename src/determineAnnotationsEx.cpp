      
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
#include <math.h> 
#include <algorithm>
#include <numeric>
#include <Rcpp.h>
#include <unordered_map>
#include "scmp.h"

exhaustiveAnnotations determineAnnotationsEx(const std::vector<std::vector<double>>& dataValsIn,
					     unsigned long colNum,
					     double dipThreshold,
					     const std::vector<std::vector<int>>&resValsIn,
					     const bool useRestrictedVals,
					     const std::vector<columnSummary>& dataSummary,
					     int& rootsInForest,
					     double depthScoreThreshold) {
  std::vector<std::vector<double>> sortedDataVals;
  std::vector<int> restrictObs;
  std::vector<double> valsForP;
  if (useRestrictedVals) {
    // if values are restricted (and so have a score of 1 or 2), omit them from the sorted data.
    // assumption: phenotypical boundaries are determined only by antimodes of **non-zero** cytof observations.
    for (auto i = 0; i != resValsIn.size(); ++i) {
      valsForP.clear();
      restrictObs = resValsIn[i];
      for (auto j = 0; j != restrictObs.size(); ++j) {
	if (restrictObs[j]==0) {
	  valsForP.push_back(((dataValsIn[i])[j]));
	}
      }
      sortedDataVals.push_back(valsForP);
    }
  }
  else {
    //if no values are restricted, use all available data.
    sortedDataVals = dataValsIn;
  }

  //copy & sort the raw data
  for (auto it = sortedDataVals.begin(); it != sortedDataVals.end(); ++it) {
    std::sort((*it).begin(),(*it).end());
  }
  std::vector<bool> labelColumn(colNum,false);
  std::vector<std::vector<double>> cutPointsOut;
  cutPointsOut.resize(colNum);

  std::vector<double> curSortData, depthNum, depthDenom;
  columnSummary curSum;
  
  //first, count how many roots there are in the forest.
  for (int i = 0; i != colNum; ++i) {
    curSortData.clear();
    curSortData = sortedDataVals[i];
    if (singleDip(curSortData) < dipThreshold) 
      ++rootsInForest;
  }
  
  double tmpDub = std::ceil((rootsInForest/2.0));
  rootsInForest = static_cast<int>(tmpDub);
  std::vector<int> currentDepths;
  int cDepth = 0, bestCutByDepthScore = 0;
  double overallDepthScore = 0, bestDepthScore = 0, currentDepthScore = 0;
  double offsetDepth, depthPenalty, cutEstimate, minDepth, dwSum;
  std::unordered_map<int, std::vector<double>> selCuts;
  std::vector<int> selDepths;
  std::vector<double> obsCutVec, cutEstVec, cutsOut, depthWeight;
  for (int numCurCol = 0; numCurCol != colNum; ++numCurCol) {
    overallDepthScore = 0;
    bestDepthScore = 0;
    curSum = dataSummary[numCurCol];
    for (auto p : curSum.depthMap) {
      //depthMap key is the number of gates, value is a vector of int depths.
      depthNum.clear();
      depthDenom.clear();
      currentDepths = p.second;
      for (auto depthIndex = 0; depthIndex != currentDepths.size(); ++depthIndex) {
	depthNum.push_back(1.0);
	cDepth = currentDepths[depthIndex];
	offsetDepth = (static_cast<double>(cDepth)-1.0);
	depthPenalty = std::pow(2.0,offsetDepth);
	if (cDepth != 1) 
	  depthDenom.push_back((depthPenalty*rootsInForest));
	else 
	  depthDenom.push_back(depthPenalty);
      }
      currentDepthScore = 0;
      for (auto dRatioIndex = 0; dRatioIndex != depthDenom.size(); ++dRatioIndex) 
	currentDepthScore += (depthNum[dRatioIndex]/depthDenom[dRatioIndex]);
      overallDepthScore += currentDepthScore;
      if (currentDepthScore > bestDepthScore) {
	bestDepthScore = currentDepthScore;
	bestCutByDepthScore = p.first;
      }
    }
    if (overallDepthScore > depthScoreThreshold) {
      cutEstVec.clear();
      labelColumn[numCurCol]=true;
      selCuts = curSum.gateMaps[bestCutByDepthScore];
      selDepths = curSum.depthMap[bestCutByDepthScore];
      minDepth = *std::min_element(selDepths.begin(),selDepths.end());
      depthWeight.clear();
      for (auto selectedIndex = 0; selectedIndex != selDepths.size(); ++selectedIndex) {
	cDepth = selDepths[selectedIndex];
	offsetDepth = (static_cast<double>(cDepth)-minDepth);
	depthPenalty = (1.0/std::pow(2.0,offsetDepth));
	depthWeight.push_back(depthPenalty);
      }
      dwSum = 0;
      for (auto v : depthWeight)
	dwSum += v;
      for (auto p : selCuts) {
	obsCutVec = p.second;
	cutEstimate = 0;
	for (auto j = 0; j != obsCutVec.size(); ++j) 
	  cutEstimate += (obsCutVec[j] * ((depthWeight[j])/dwSum));
	cutEstVec.push_back(cutEstimate);
      }
      //sort since using unordered_map -- cut-points in increasing order.
      std::sort(cutEstVec.begin(),cutEstVec.end());
      for (auto j = 0; j != cutEstVec.size(); ++j) { 
	(cutPointsOut[numCurCol]).push_back(cutEstVec[j]);
      }
    }
  }
  exhaustiveAnnotations retVal = {labelColumn, cutPointsOut};
  return retVal;
}
