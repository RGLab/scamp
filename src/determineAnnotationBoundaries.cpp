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
#include <map>

std::vector<double> determineAnnotationBoundaries(const std::vector<gateInfo>& placements, unsigned long annCol, bool printDebugInfo){
  //first count the occurance of different gating numbers for the annotation column.
  std::map<unsigned long, unsigned long> gateCounter;
  gateInfo currentGI;
  for (auto i = 0; i != placements.size(); ++i) {
    currentGI = placements[i];
    if ((currentGI.colNumber)==annCol) {
      ++gateCounter[(currentGI.numGates)];
    }
  }
  //return a vacuous result.
  if (gateCounter.empty()) {
    std::cout << "No gates recorded for annotation column " << annCol << std::endl;
    std::cout << "Returning vacous annotation boundary of 0." << std::endl;
    std::vector<double> vacuousResult = {0.0};
    return vacuousResult;
  }
  
  //next, determine which number of gates occurs most frequently.
  unsigned long currentMax = 0, mostFreqGate = 0;
  for (auto p : gateCounter) {
    if (p.second > currentMax) {
      currentMax = p.second;
      mostFreqGate = p.first;
    }
  }
  //next, build up distribution of each gate type for most frequent gates.
  std::vector<double> currentGates;
  double gateVal;
  std::vector<std::vector<double>> selGates; 
  selGates.resize((mostFreqGate-2)); //gates record min and max too.
  for (auto i = 0; i != placements.size(); ++i) {
    currentGI = placements[i];
    if (((currentGI.colNumber)==annCol) && (currentGI.numGates == mostFreqGate)){
      currentGates = currentGI.gates;
      for (auto j = 1; j != (currentGates.size() - 1); ++j) {
	gateVal = currentGates[j];
	(selGates[(j-1)]).push_back(gateVal);
      }
    }
  }
  
  if (printDebugInfo) {
    std::cout << "annCol: " << annCol << std::endl;
    std::cout << "mostFreqGate: " <<  mostFreqGate << std::endl;
    std::cout << "numGates: " <<  (selGates[0]).size() << std::endl;
  }

  //finally, compute median of each type of selected gate, and return them as the annotation boundary.
  std::vector<double> gatesOut((selGates.size()),0);
  std::vector<double> probVal = {0.5};
  std::vector<double> qnt;
  for (auto k = 0; k != selGates.size(); ++k) {
    qnt = rQuantile(selGates[k],probVal);
    if (printDebugInfo) {
      std::cout << "qnt " << k << ": " << qnt[0] << std::endl;
    }
    gatesOut[k] = qnt[0];
  }
  return gatesOut;
}
