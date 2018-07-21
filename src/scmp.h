#ifndef CPP_SCMP_H
#define CPP_SCMP_H

#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <Rcpp.h>
#include <map>
#include <unordered_map>

void cppdip(const double x[],
	    const int*,
	    double*,
	    std::vector<int>&,
	    int*,
	    int*,
	    int*,
	    int*,
	    int*,
	    const int*,
	    const int*);

double singleDip(const std::vector<double> &);
double doubleDip(const std::vector<double> &);
double tripleDip(const std::vector<double> &);


struct predictionLabelInfo {
  std::vector<double> labelGates;
  bool useCustomLabel;
};

struct stringInfo {
  std::vector<double> string;
  std::vector<int> knotsind;
  std::vector<double> knotst;
  std::vector<double> knotsy;
  int nknots;
  int nmax;
};


struct crParameters {
  bool clusterOverflow;
  int maxDepth;
  long maxClusters;
  double localDipThresh;
  int minClusterSize;
};


struct gateInfo {
  int colNumber;
  unsigned long numGates;
  std::vector<double> gates;
  int gateDepth;
};

struct searchResults {
  std::vector<std::vector<long>> candidates;
  std::vector<gateInfo> gateLocations;
  bool abortIteration;
};

void local_density(const std::vector<double>&, std::vector<int>&, long, long, long);

stringInfo tautString(const std::vector<double>&, 
		      const std::vector<double>&, 
		      const std::vector<double>&, 
		      const std::vector<double>&, 
		      double,
		      double,
		      long,
		      int); 

void easymax(std::vector<double>& , long ,long ,long ,long*,long*,double*);

void difficultmax(std::vector<double>&, 
		  long, 
		  long, 
		  long, 
		  long*, 
		  long*,
		  double*);

std::vector<double> kkuiper(std::vector<double>&, long, int);
std::vector<double> cppApprox(std::vector<double>&, std::vector<double>&, std::vector<double>&);

stringInfo cpPmden(const std::vector<double>&);
bool isLocalMinTS(double&,double&,double&);
bool isLocalMaxTS(double&,double&,double&);

std::vector<double> findKmedGates(const std::vector<double>&,
				  const std::vector<int>&,
				  const int);
std::vector<double> tsGates(const std::vector<double>&, int);
std::vector<double> tsModeEstimate(const std::vector<double>&);

std::vector<double> absSampleLmoments(std::vector<double>, int, int);

std::vector<int> selectCandidateClusters(std::vector<std::vector<long>>&,
					 std::vector<std::vector<double>>&,
					 const bool&);

std::vector<double> rQuantile(const std::vector<double>&, std::vector<double>);
double medianAbsoluteDeviation(const std::vector<double>&);

std::vector<bool> determineAnnotationColumns(const std::vector<std::vector<double>>&,
					     unsigned long,
					     double,
					     const std::vector<std::vector<int>>&,
					     const bool);

std::vector<double> determineAnnotationBoundaries(const std::vector<gateInfo>&, unsigned long, bool);

std::vector<std::vector<int>> annotateCluster(const std::vector<std::vector<double>>& ,
					      const std::vector<int>& ,
					      std::vector<bool>& ,
					      std::vector<std::vector<double>>&,
					      std::vector<std::string>&,
					      const std::vector<double>&);

std::vector<std::vector<int>> parallelAnnotateCluster(const std::vector<std::vector<double>>&,
						      const std::vector<int>&,
						      const std::vector<bool>&,
						      const std::vector<std::vector<double>>&,
						      const std::vector<std::string>&,
						      const std::vector<double>&,
						      const unsigned long);


std::vector<std::vector<double>> addNoiseToDataMatrix(const std::vector<std::vector<double>>&,
						      unsigned long,
						      double,
						      unsigned long long&);


Rcpp::List cppNoisyScamp(Rcpp::NumericMatrix&,
			 double,
			 int,
			 bool,
			 int,
			 long,
			 Rcpp::StringVector,
			 long,
			 bool,
			 bool,
			 const std::vector<double>&,
			 int,
			 bool,
			 Rcpp::NumericMatrix&,
			 bool,
			 Rcpp::List,
			 int,
			 std::string,
			 bool,
			 double,
			 bool,
			 double,
			 unsigned long long);
		

std::vector<std::vector<int>> scamp(const std::vector<std::vector<double>>&,
				    double,
				    int,
				    bool,
				    int,
				    long,
				    bool,
				    std::vector<std::string>,
				    long,
				    bool,
				    bool,
				    const std::vector<double>&,
				    unsigned long,
				    bool,
				    const std::vector<std::vector<int>>&,
				    bool,
				    Rcpp::List,
				    const double&,
				    double,
				    unsigned long long&);


//cleaning up the recursion; now stack implemented

searchResults findCandidateClusters(const std::vector<std::vector<double>>&,
				    const std::vector<std::vector<int>>&,
				    const double&,
				    const int&,
				    const bool&,
				    const int&,
				    const long&,
				    const long&,
				    const bool&,
				    const bool&,
				    const double&,
				    const bool&,
				    const bool&,
				    unsigned long,
				    unsigned long long&); 
 

searchResults candidateClusterSearch(const std::vector<std::vector<double>>&,
				     const std::vector<std::vector<int>>&,
				     double,
				     int,
				     bool,
				     int,
				     long,
				     long,
				     bool,
				     bool,
				     bool,
				     unsigned long long,
				     bool,
				     int);


#endif

