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
#include <map>
#include "scmp.h"
#include <fstream>
#include <cstdio> //for std::rename
#include <thread>

#ifndef _NOISY_DEBUG_USER_ANNOTATIONS_
#define _NOISY_DEBUG_USER_ANNOTATIONS_ false
#endif

int localGCD(int n1, int n2) {
  //Euclid...
  int temp;
  while (n2 != 0) {
    temp = n1;
    n1 = n2;
    n2 = temp % n2;
  }
  return n1;
}

std::vector<int> topTwoLabel(const std::vector<int>& vecOne,
			     const std::vector<int>& vecTwo){
  //given two label vectors of ints, vecOne and vecTwo,
  //combine them into a single labeling vector.
  std::vector<int> resVec;
  int tempInt1,tempInt2,tempInt3,tempInt4, tempNum, tempDenom, tempGCD;
  for (auto n = 0; n != vecOne.size(); n+=2) {
    tempInt1 = vecOne[n];
    tempInt2 = vecOne[(n+1)];
    tempInt3 = vecTwo[n];
    tempInt4 = vecTwo[(n+1)];
    if (tempInt2 == 0) {
      //plural vote did not annotate: do no add.
      resVec.push_back(0);
      resVec.push_back(0);
    }
    else if (tempInt4 == 0) {
      //secondary vote did not annotate: defer to plurality.
      resVec.push_back(tempInt1);
      resVec.push_back(tempInt2);
    }
    else {
      //disagreement between primary and secondary vote.
      //combine the labels.
      tempNum = ((tempInt1*tempInt4)+(tempInt2*tempInt3));
      tempDenom = (2 * tempInt2 * tempInt4);
      tempGCD = localGCD(tempNum,tempDenom);
      tempNum /= tempGCD;
      tempDenom /= tempGCD;
      resVec.push_back(tempNum);
      resVec.push_back(tempDenom);
    }
  }
  return resVec;
}

std::vector<int> topThreeLabel(const std::vector<int>& vecOne,
				  const std::vector<int>& vecTwo,
				  const std::vector<int>& vecThree){
  //given three label vectors of ints, vecOne, vecTwo, and vecThree,
  //combine them into a single labeling vector.
  std::vector<int> resVec;
  int tempInt1,tempInt2,tempInt3,tempInt4,tempInt5,tempInt6, tempNum, tempDenom, tempGCD;
  for (auto n = 0; n != vecOne.size(); n+=2) {
    tempInt1 = vecOne[n];
    tempInt2 = vecOne[(n+1)];
    tempInt3 = vecTwo[n];
    tempInt4 = vecTwo[(n+1)];
    tempInt5 = vecThree[n];
    tempInt6 = vecThree[(n+1)];
    if ((tempInt2 == 0) && (tempInt4 == 0) && (tempInt6 == 0)) {
      resVec.push_back(0);
      resVec.push_back(0);
    }
    else if ((tempInt2 == 0) && (tempInt4 == 0)){
      //third-place only label -- do not add.
      resVec.push_back(0);
      resVec.push_back(0);
    }
    else if ((tempInt2 == 0) && (tempInt6 == 0)){
      //second-place only label -- do not add.
      resVec.push_back(0);
      resVec.push_back(0);
    }
    else if ((tempInt4 == 0) && (tempInt6 == 0)){
      //first-place only label -- defer to it.
      resVec.push_back(tempInt1);
      resVec.push_back(tempInt2);
    }
    else if (tempInt2 == 0) {
      //second- and third- place label -- combine.
      tempNum = ((tempInt4*tempInt5)+(tempInt3*tempInt6));
      tempDenom = (2 * tempInt4 * tempInt6);
      tempGCD = localGCD(tempNum,tempDenom);
      tempNum /= tempGCD;
      tempDenom /= tempGCD;
      resVec.push_back(tempNum);
      resVec.push_back(tempDenom);
    }
    else if (tempInt4 == 0) {
      //first- and third- place label -- combine.
      tempNum = ((tempInt1*tempInt6)+(tempInt2*tempInt5));
      tempDenom = (2 * tempInt2 * tempInt6);
      tempGCD = localGCD(tempNum,tempDenom);
      tempNum /= tempGCD;
      tempDenom /= tempGCD;
      resVec.push_back(tempNum);
      resVec.push_back(tempDenom);
    }
    else if (tempInt6 == 0) {
      //first- and second- place label -- combine.
      tempNum = ((tempInt1*tempInt4)+(tempInt2*tempInt3));
      tempDenom = (2 * tempInt2 * tempInt4);
      tempGCD = localGCD(tempNum,tempDenom);
      tempNum /= tempGCD;
      tempDenom /= tempGCD;
      resVec.push_back(tempNum);
      resVec.push_back(tempDenom);
    }
    else {
      //all votes label -- combine.
      tempNum = ((tempInt1*tempInt4*tempInt6)+(tempInt2*tempInt3*tempInt6)+(tempInt2*tempInt4*tempInt5));
      tempDenom = (3 * tempInt2 * tempInt4 * tempInt6);
      tempGCD = localGCD(tempNum,tempDenom);
      tempNum /= tempGCD;
      tempDenom /= tempGCD;
      resVec.push_back(tempNum);
      resVec.push_back(tempDenom);
    }
  }
  return resVec;
}

std::vector<std::string> assignLabelToScoreVector(const std::vector<std::vector<int>>& scoredClustering,
						  const std::vector<std::string>& clusterAnnotations,
						  const std::vector<std::vector<std::string>>& finalLabels){
  std::string annotationString = "";
  std::vector<int> tempIntVector;
  int tempColNum, setSize, setScore;
  auto numObservations = scoredClustering.size();
  std::vector<std::string> labelVector(numObservations,annotationString);
  for (auto i=0; i != numObservations; ++i) {
    annotationString = "";
    tempIntVector = scoredClustering[i];
    for (auto j=0; j != tempIntVector.size(); j+=2) {
      tempColNum = j/2;
      if ((tempIntVector[j]==0) && (tempIntVector[(j+1)]==0)) {
	continue;
      }
      else {
	setScore = tempIntVector[j];
	setSize = (tempIntVector[(j+1)]-1);
	if ((setSize <= 7) && (setScore <= (setSize+1))){
	  annotationString += clusterAnnotations[tempColNum] + "~" + ((finalLabels[setSize])[setScore]) + "~";
	}
	else {
	  annotationString += clusterAnnotations[tempColNum] + "~AnnotationError~"; 
	}
      }
    }
    if (annotationString != "")  {
      labelVector[i]=annotationString;
    }
    else {
      labelVector[i]="Uncertain";
    }
  }
  return labelVector;
}


std::vector<std::string> assignNumLabelToScoreVector(const std::vector<std::vector<int>>& scoredClustering,
						     const std::vector<std::string>& clusterAnnotations){
  std::string annotationString = "";
  std::vector<int> tempIntVector;
  int tempColNum, setSize, setScore;
  auto numObservations = scoredClustering.size();
  std::vector<std::string> labelVector(numObservations,annotationString);
  for (auto i=0; i != numObservations; ++i) {
    annotationString = "";
    tempIntVector = scoredClustering[i];
    for (auto j=0; j != tempIntVector.size(); j+=2) {
      tempColNum = j/2;
      if ((tempIntVector[j]==0) && (tempIntVector[(j+1)]==0)) {
	continue;
      }
      else {
	setScore = (tempIntVector[j]+1);
	setSize = (tempIntVector[(j+1)]+1);
	if (setScore <= setSize) {
	  annotationString += (clusterAnnotations[tempColNum] + "~" + std::to_string(setScore) + "~" + std::to_string(setSize) + "~");
	}
	else {
	  annotationString += clusterAnnotations[tempColNum] + "~AnnotationError~"; 
	}
      }
    }
    if (annotationString != "")  {
      labelVector[i]=annotationString;
    }
    else {
      labelVector[i]="Uncertain";
    }
  }
  return labelVector;
}




// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List cppNoisyScamp(Rcpp::NumericMatrix& rawDataMatrix,
			 double dipT,
			 int clusterLB,
			 bool repeatsAllowed,
			 int maxSearchDepth,
			 long maxClusterNum,
			 Rcpp::StringVector userAnnotations,
			 long maxNumberOfGates,
			 bool randomSearch,
			 bool randomResidualSearch,
			 const std::vector<double>& finalAnnotationQs,
			 int numThreadsRequested,
			 bool useRestrictedValue,
			 Rcpp::NumericMatrix& restrictedValueMatrix,
			 bool useForestValues,
			 Rcpp::List forestValues,
			 int numberOfScampIterations,
			 std::string outputDirectory,
			 bool verboseOutput,
			 double maxSearchTime,
			 bool debugScampRun,
			 double gaussianScale,
			 unsigned long long randomSeed,
			 double depthScoreThreshold)
{
  //Determine number of threads to use.
  unsigned long numThreadsToUse;
  unsigned long const maxNumThreadsPossible = std::thread::hardware_concurrency();
  if (numThreadsRequested < 1) {
    if (debugScampRun) {
      std::cout << "User requested max threads." << std::endl;
      std::cout << "Hardware supports maximum of " << maxNumThreadsPossible << " threads." << std::endl;
    }
    numThreadsToUse = maxNumThreadsPossible;
    if (debugScampRun) {
      std::cout << "Using " << numThreadsToUse << " threads." << std::endl;
    }
  }
  else {
    if (debugScampRun) {
      std::cout << "User requested " << numThreadsRequested << ". " << std::endl;
      std::cout << "Hardware supports maximum of " << maxNumThreadsPossible << " threads." << std::endl;
    }
    numThreadsToUse = numThreadsRequested;
    if (debugScampRun) {
      std::cout << "Using " << numThreadsToUse << " threads." << std::endl;
    }
  }
  
  //Copy R matrices into C++ data structures.
  std::vector<std::vector<double>> dataValsIn;
  std::vector<std::vector<int>> restrictedVals;
  std::vector<int> tmpFillVector(2,0);
    
  std::vector<double> tmpV;
  std::vector<int> tmpV2;
  auto lColNum = rawDataMatrix.ncol();
  for (auto i = 0; i != lColNum; ++i) {
    tmpV.clear();
    tmpV2.clear();
    Rcpp::NumericMatrix::Column tmp = rawDataMatrix(Rcpp::_,i);
    for (auto j : tmp)
      tmpV.push_back(j);
    dataValsIn.push_back(tmpV);
    if (useRestrictedValue) {
      Rcpp::NumericMatrix::Column tmp2 = restrictedValueMatrix(Rcpp::_,i);
      for (auto k : tmp2) {
	if (k == 1) {
	  tmpV2.push_back(1);
	}
	else if (k ==2) {
	  tmpV2.push_back(2);
	}
	else {
	  //default to unrestricted in the event a non-{0,1,2} int is in set.
	  //R wrapper is responsible for making sure this does not happen.
	  tmpV2.push_back(0);
	}
      }
      restrictedVals.push_back(tmpV2);
    }
    else {
      //pass along null restrictions to simplify later function calls.
      restrictedVals.push_back(tmpFillVector);
    }
  }

  //Next, copy user annotations into C++ vector
  std::vector<std::string> clusterAnnotations(userAnnotations.size());
  for (auto i = 0; i != userAnnotations.size(); ++i) {
    clusterAnnotations[i] = userAnnotations(i);
    //the clusterAnnotations are by default the column names of the R data matrix.
    //however, the user can provide custom labels. check those values here.
    if (_NOISY_DEBUG_USER_ANNOTATIONS_) {
      std::cout << "Annotation " << i << ": " << clusterAnnotations[i] << std::endl;
    }
  }

  //label containters
  std::vector<std::vector<std::string>> finalLabels(8);
  std::vector<std::string> twoLabels = {"Lowest","Highest"};
  std::vector<std::string> threeLabels = {"Lowest","Medium","Highest"};
  std::vector<std::string> fourLabels = {"Lowest","MediumLow","MediumHigh","Highest"};
  std::vector<std::string> fiveLabels = {"Lowest","MediumLow","Medium","MediumHigh","Highest"};
  std::vector<std::string> sixLabels = {"Lowest","Low","MediumLow","MediumHigh","High","Highest"};
  std::vector<std::string> sevenLabels = {"Lowest","Low","MediumLow","Medium","MediumHigh","High","Highest"};
  std::vector<std::string> eightLabels = {"Lowest","VeryLow","Low","MediumLow","MediumHigh","High","VeryHigh","Highest"};
  std::vector<std::string> nineLabels = {"Lowest","VeryLow","Low","MediumLow","Medium","MediumHigh","High","VeryHigh","Highest"};
  finalLabels[0]=twoLabels;
  finalLabels[1]=threeLabels;
  finalLabels[2]=fourLabels;
  finalLabels[3]=fiveLabels;
  finalLabels[4]=sixLabels;
  finalLabels[5]=sevenLabels;
  finalLabels[6]=eightLabels;
  finalLabels[7]=nineLabels;
  
  //containers for scamp
  int numIter = numberOfScampIterations;
  int elapsedIter;
  std::map<std::vector<int>, int> assignmentCount, tempMap;
  auto numberOfDataObs = dataValsIn[0].size();
  std::vector<std::map<std::vector<int>, int>> allScampResults(numberOfDataObs,assignmentCount);
  std::vector<std::vector<int>> scampResult, tmpStringVec;
  std::vector<int> uncertainVec((2*lColNum),0);
  std::vector<std::vector<int>> finalScampResult(numberOfDataObs,uncertainVec);
  std::vector<std::vector<int>> maxScampResult(numberOfDataObs,uncertainVec);
  std::vector<std::string> voteScampResult, mlabScampResult;
  std::vector<int> finalScampCounts(numberOfDataObs,0);
  std::string iterFile = (outputDirectory+"/scampResults/scampMetaData.txt");
  std::string mostRecent = (outputDirectory+"/scampResults/scampLatestClusteringVote.txt");
  std::string lastIter = (outputDirectory+"/scampResults/scampPreviousClusteringVote.txt");
  std::string currentClustering = (outputDirectory+"/scampResults/scampLatestClustering.txt");
  std::string previousClustering = (outputDirectory+"/scampResults/scampPreviousClustering.txt");
  std::string currentMax = (outputDirectory+"/scampResults/scampLatestMax.txt");
  std::string previousMax = (outputDirectory+"/scampResults/scampPreviousMax.txt");

  std::string iterString;
  std::ofstream ofile; 
  std::vector<std::pair<std::vector<int>, int>> votePairs;
  std::pair<std::vector<int>, int> firstPair, secondPair, thirdPair, tempPair;
  double firstRatio,secondRatio,thirdRatio;
  std::vector<int> tempIntVec;
  int returnCode;
  const char* cMR = mostRecent.c_str();
  const char* cLI = lastIter.c_str();
  const char* cCC = currentClustering.c_str();
  const char* cPC = previousClustering.c_str();
  const char* cCM = currentMax.c_str();
  const char* cPM = previousMax.c_str();
  bool firstWriteToFile = true;

  unsigned long long currentRandomSeed = randomSeed;
  
  //with the preparation complete, begin finding scamp clusterings.
  while (numIter) {
    scampResult = scamp(dataValsIn,dipT,clusterLB,repeatsAllowed,maxSearchDepth,
			maxClusterNum,debugScampRun,clusterAnnotations,
			maxNumberOfGates,randomSearch,randomResidualSearch,
			finalAnnotationQs,numThreadsToUse,
			useRestrictedValue,restrictedVals,useForestValues,
			forestValues,maxSearchTime,gaussianScale,currentRandomSeed,
			depthScoreThreshold);
    
    //check if scamp iteration was terminated early
    //if so, continue to next loop iteration
    tempIntVec = scampResult[0];
    if ((tempIntVec.size() == 2) && (tempIntVec[0] == -1) && (tempIntVec[1] == -1)) {
      if (debugScampRun) {
	std::cout << "Candidate cluster search failed. Abort iteration and re-noise." << std::endl;
      }
      continue;
    }

    //if scamp successful, store the iteration's clustering result in the map associated with each observation.
    for (auto i = 0; i != numberOfDataObs; ++i) {
      ++((allScampResults[i])[(scampResult[i])]); 
    }

    elapsedIter = (numberOfScampIterations-numIter+1);
    iterString = std::to_string(elapsedIter);

    //write number of scamp iterations to file.
    if (verboseOutput && ((elapsedIter % 10)==0)) { 
      ofile.open(iterFile);
      ofile << ("The clusterings are based off of " + iterString + " scamp iterations.") << std::endl;
      ofile.close();
    }

    //update clustering scores.
    for (auto i=0; i != numberOfDataObs; ++i) {
      tempMap = allScampResults[i];
      votePairs.clear();
      //copy an observation voting history in votePairs.
      for (auto itr = tempMap.begin(); itr != tempMap.end(); ++itr)
	votePairs.push_back(*itr);

      //sort vote pairs by decreasing ballot frequency.
      std::sort(votePairs.begin(), votePairs.end(),
		[=](std::pair<std::vector<int>, int>& a, std::pair<std::vector<int>, int>& b){return a.second > b.second;});

      if (verboseOutput) {
	//write the sorted votePairs to file
	if (i == 0) {
	  //if scamp has executed at least one time, store the previous result.
	  if (numIter != numberOfScampIterations) {
	    returnCode = std::rename(cMR, cLI); 
	    if (returnCode) { std::perror("Error renaming"); return 1;}
	  }
	
	  //if its the first votePair, clobber the last itertation.
	  ofile.open(mostRecent);
	  for (auto i=0; i != votePairs.size(); ++i) {
	    tempPair = votePairs[i];
	    tempIntVec = tempPair.first;
	    for (auto k=0; k != (tempIntVec.size()-1); ++k)
	      ofile << tempIntVec[k] << ":";
	    ofile << tempIntVec[(tempIntVec.size()-1)] << "," << tempPair.second << ((i == (votePairs.size()-1)) ? " " : ",");
	  }
	  ofile << std::endl;
	  ofile.close();
	}
	else {
	  //not first votePair, so append.
	  ofile.open(mostRecent, std::ofstream::app);
	  for (auto i=0; i != votePairs.size(); ++i) {
	    tempPair = votePairs[i];
	    tempIntVec = tempPair.first;
	    for (auto k=0; k != (tempIntVec.size()-1); ++k)
	      ofile << tempIntVec[k] << ":";
	    ofile << tempIntVec[(tempIntVec.size()-1)] << "," << tempPair.second << ((i == (votePairs.size()-1)) ? " " : ",");
	  }
	  ofile << std::endl;
	  ofile.close();
	}
      }
      
      //update clustering with maxvote/runoff vote.
      firstPair = votePairs[0];
      firstRatio = static_cast<double>(firstPair.second)/static_cast<double>(elapsedIter);
      maxScampResult[i] = firstPair.first;
      if (firstRatio >= 0.5) {
	finalScampResult[i] = firstPair.first;
      }
      else if (firstRatio >= 0.3) {
	secondPair = votePairs[1];
	secondRatio = static_cast<double>(secondPair.second)/static_cast<double>(elapsedIter);
	if (secondRatio >= 0.2) {
	  finalScampResult[i] = topTwoLabel(firstPair.first,secondPair.first);
	}
	else {
	  finalScampResult[i] = firstPair.first;
	}
      }
      else if (firstRatio >= 0.2) {
	secondPair = votePairs[1];
	secondRatio = static_cast<double>(secondPair.second)/static_cast<double>(elapsedIter);
	thirdPair = votePairs[2];
	thirdRatio = static_cast<double>(thirdPair.second)/static_cast<double>(elapsedIter);
	if (thirdRatio >= 0.15) {
	  finalScampResult[i] = topThreeLabel(firstPair.first,secondPair.first,thirdPair.first);
	}
	else if (secondRatio >= 0.175) {
	  finalScampResult[i] = topTwoLabel(firstPair.first,secondPair.first);
	}
	else {
	  finalScampResult[i] = firstPair.first;
	}
      }
      else {
	finalScampResult[i] = uncertainVec;
      }
    }

    if (useForestValues) {
      voteScampResult = assignNumLabelToScoreVector(finalScampResult,clusterAnnotations);
      mlabScampResult = assignNumLabelToScoreVector(maxScampResult,clusterAnnotations);
    }
    else {
      //transform the final scamp result into a vector of strings, and write to file.
      voteScampResult = assignLabelToScoreVector(finalScampResult,clusterAnnotations,finalLabels);
      mlabScampResult = assignLabelToScoreVector(maxScampResult,clusterAnnotations,finalLabels);
    }    
    //write the labeled scamp results to file
    if (verboseOutput && ((elapsedIter % 50)==0)) { 
      for (auto i=0; i != numberOfDataObs; ++i) {
	if (i == 0) {
	  //once again, store the previous result.
	  if (!firstWriteToFile) {
	    returnCode = std::rename(cCC, cPC); 
	    if (returnCode) { std::perror("Error renaming"); return 1;}
	    returnCode = std::rename(cCM, cPM); 
	    if (returnCode) { std::perror("Error renaming"); return 1;}
	  }
	  ofile.open(currentClustering);
	  ofile << voteScampResult[i] << std::endl;
	  ofile.close();
	  ofile.open(currentMax);
	  ofile << mlabScampResult[i] << std::endl;
	  ofile.close();
	}
	else {
	  //not first votePair, so append.
	  ofile.open(currentClustering, std::ofstream::app);
	  ofile << voteScampResult[i] << std::endl;
	  ofile.close();
	  ofile.open(currentMax, std::ofstream::app);
	  ofile << mlabScampResult[i] << std::endl;
	  ofile.close();
	}
      }
      if (firstWriteToFile) {
	firstWriteToFile = false;
      }
    }
    --numIter;
  }
  return Rcpp::List::create(Rcpp::Named("RunOffVote") = voteScampResult,
			    Rcpp::Named("MaxVote") = mlabScampResult);
}
