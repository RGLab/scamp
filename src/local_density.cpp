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

/*
 *This file is a simple port of the fortran subroutine denlocal.f, taken from the ftnonpar library  by Laurie Davies and Arne Kovac.
 *The orginal work is found in version 0.1-88 of that package.
 *The maintainer of that library, Arne Kovac <A.Kovac at bristol.ac.uk>, is not affiliated with this package.
 *Minor changes have for use within scamp. 
 *--Evan Greene <egreene@fredhutch.org>, 9/19/2017
 */

#include <Rcpp.h>
#include <vector>
#include <iostream>
#include "scmp.h"


void local_density(const std::vector<double>& cf, std::vector<int>& kni, long N, long ICOMAX, long IRMAX) {
  long I = 1;
  long J = 1;
  long ICR = 0;
  long IC =0;
  bool localDensity = true;
  bool stillLooking = true;
  while (localDensity) {
    IC = 0;
    stillLooking = true;
    while (stillLooking) {
      if (I > N) {
	stillLooking = false;
      }
      else {
	if (cf[(I-1)] <= (double(J)/double(N))) {
	  if (ICR >= IRMAX) {
	    kni[(I-1)] = 1;
	    kni[(I-2)] = 1;
	  }
	  ICR = 0;
	  IC = IC + 1;
	  if (I < N) {
	    I = I+1;
	  }
	  else {
	    stillLooking = false;
	  }
	}
	else {
	  stillLooking = false;
	}
      }
    }
    if (IC == 0) {
      ICR = ICR + 1;
    }
    if (IC >= ICOMAX) {
      for (long m = (I-IC); m < I; m++) {
	kni[m]=1;
      }
    }
    if (I == N) {
      localDensity = false;
    }
    J = J+1;
  }
}



