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


//This file is a simplified port of the function samlmu.s, found in the lmom library version 2.6.
//The author and maintainer of that package, J. R. M. Hosking <jrmhosking at gmail.com>, is not affiliated with the scamp project.
//C++ port written 8/7/2017 by Evan Greene <egreene@fredhutch.org>.
//Reduced functionality for use in SCMP scoring function. For trimming and other functionality, see the lmom package on CRAN.
//Following comment at the head of lmom.r file in "lmom" package
/***********************************************************************
 *                                                                     *
 *  R code for R package "lmom"                                        *
 *                                                                     *
 *  J. R. M. HOSKING <jrmhosking@gmail.com>                            *
 *                                                                     *
 *  Version 2.5    February 2015                                       *
 *                                                                     *
 ***********************************************************************/
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "scmp.h"

double _lMom_sumDV(const std::vector<double>& dVec) {
  double lmTmp = 0.0;
  for (auto myv : dVec)
    lmTmp += myv;
  return lmTmp;
}

std::vector<double> absSampleLmoments(std::vector<double> sortedDataVector, int numMoments, int trim)
{
  // this function assumes sortedDataVector is sorted data are sorted
  long numObs = sortedDataVector.size();
  
  int trimLow = trim, trimHigh = trim;

  //require numMoments+trimLow+trimHigh < n
  int maxmom = numMoments+trimLow+trimHigh;   // Number of untrimmed L-moments to compute
  //
  // Compute (untrimmed) sample L-moments.  The (j in 3:maxmom) loop computes
  // discrete Legendre polynomials recursively and uses them as weights for
  // the ordered observations.
  //
  std::vector<double> lMoments(maxmom,0.0);
  double tmp = _lMom_sumDV(sortedDataVector);
  
  lMoments[0] = tmp/(double(numObs));
  std::vector<double> symmetricMesh(numObs,0.0);
  double lb = -numObs-1.0;
  std::generate(symmetricMesh.begin(),symmetricMesh.end(), [&lb] { return lb += 2.0;});

  std::vector<double> legendrePoly = symmetricMesh;
  for (auto &v : legendrePoly) {
    v = v/(numObs-1);
  }

  std::vector<double> legendrePoly1(numObs,1.0);
  std::vector<double> legendrePoly2(numObs,0.0);

  std::vector<double> tmpDV(numObs,0.0);
  for (auto i = 0; i != tmpDV.size(); ++i) 
    tmpDV[i] = legendrePoly[i]*sortedDataVector[i];
  tmp = _lMom_sumDV(tmpDV);
  lMoments[1] = tmp/numObs;

  for (int j = 3; j != (maxmom+1); ++j) {
    legendrePoly2 = legendrePoly1;
    legendrePoly1 = legendrePoly;
    for (auto i = 0; i != legendrePoly.size(); ++i)  {
      legendrePoly[i] = (((2*j-3)*symmetricMesh[i]*legendrePoly1[i]-(j-2)*(numObs+j-2)*legendrePoly2[i])/((j-1.0)*(numObs-j+1.0)));
      tmpDV[i] = legendrePoly[i] * sortedDataVector[i];
    }
    tmp = _lMom_sumDV(tmpDV);
    lMoments[(j-1)] = tmp/(numObs);
  }

  //trim upper
  if (trim) {
    for (int j = 1; j != (trim+1); ++j) {
      --maxmom;
      for (int r = 1; r != (maxmom - j + 1); ++r) {
	lMoments[(r-1)] = ((lMoments[(r-1)]*(r+j) - lMoments[r]*(r+1))/static_cast<double>(((2*r)+j-1)));
      }
      lMoments.erase((lMoments.begin() + (maxmom - j)),lMoments.end());
    }
  }

  //trim lower. this implementation might introduce rounding errors in trimmed L-moments of order greater than 4,
  //and trim value >4. however, this function should only be called to procude second, third, and fourth order trimmed
  //L-moments. For more general functionality, use the R package lmom by Hoskings.
  if (trim) {
    for (int j = 1; j != (trim+1); ++j) {
      --maxmom;
      for (int r = 1; r != (lMoments.size() - 1); ++r) {
	lMoments[(r-1)] = ((lMoments[(r-1)]*(r+j+trim) + ((lMoments[r]*(r+1)*(r+trim))/static_cast<double>(r)))/static_cast<double>(((2*r)+j+trim-1)));
      }
      lMoments.erase((lMoments.begin() + (lMoments.size() - 1)),lMoments.end());
    }
  }

  //divide out second L-Moment;
  for (auto i = 0; i != lMoments.size(); ++i)
    if (i > 1)
      lMoments[i] = (lMoments[i]/lMoments[1]);


  //take absolue value -- scores will only be based on magnitude
  for (auto &lmv : lMoments)
    lmv = std::abs(lmv);

  return lMoments;
}



