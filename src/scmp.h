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

struct columnSummary {
  //gateMaps key is number of cut points -- the cut-point-key.
  //value is dictionary with integer keys between 0 and ((cut-point-key)-1).
  //the vector<double> corresponding to values in nested map contain observed cut-point
  //location in sort order, with position indicating tuples of cut-points.
  //example of intended use: suppose the cut-point-key is 3. the nested map then contains
  //|  0: 0.25, 0.35 <-- corresponds to the smallest observed cut-point value.
  //|  1: 0.75, 0.68 <-- corresponds to the middle observed cut-points.
  //|  2: 3.99, 4.12 <-- corresponds to the largest cut-points.
  //this means that during the search for candidate clusters, this column was cut
  //into four sub-collections with boundaries (0.25,0.75,3.99) in one search tree
  //and (0.35,0.68,4.12) in another search tree.
  std::unordered_map<int, std::unordered_map<int, std::vector<double>>> gateMaps; 
  //depthMap key is ALSO number of cut points -- once again, the cut-point-key.
  //value is a vector of ints corresponding to the depth at which the cut-point-key was observed.
  //continuing with the preceding example, suppose the cut-point-key is 3. depthMaps then contains
  //|  3:  4, 7  <-- indicates the column was split at (0.25,0.75,3.99) at depth 4 of one search tree,
  //and split at values (0.35,0.68,4.12) at depth 7 in another search tree.
  //in particular, the depth position in depthMap describes the cut-point depth in gateMaps,
  //which is used to determine the columnScore -- the final parameter.
  std::unordered_map<int, std::vector<int>> depthMap;
  double columnScore;
};

struct exhaustiveAnnotations{
  std::vector<bool> columnIsAnnotated;
  std::vector<std::vector<double>> cutPointLocations;
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

exhaustiveAnnotations determineAnnotationsEx(const std::vector<std::vector<double>>&,
					     unsigned long,
					     double,
					     const std::vector<std::vector<int>>&,
					     const bool,
					     const std::vector<columnSummary>&,
					     int&,
					     double);


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
			 unsigned long long,
			 double,
			 unsigned long,
			 unsigned long,
			 unsigned long);


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
				    unsigned long long&,
				    double,
				    unsigned long,
				    unsigned long,
				    unsigned long);

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
				    unsigned long long&,
				    unsigned long,
				    unsigned long,
				    unsigned long);


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
				     int,
				     unsigned long,
				     unsigned long,
				     unsigned long);




#endif

