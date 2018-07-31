\name{cppNoisyScamp}
\alias{cppNoisyScamp}
\title{cppNoisyScamp}
\usage{
cppNoisyScamp(rawDataMatrix, dipT, clusterLB, repeatsAllowed,
maxSearchDepth, maxClusterNum, userAnnotations, maxNumberOfGates,
randomSearch, randomResidualSearch, finalAnnotationQs,
numThreadsRequested, useRestrictedValue, restrictedValueMatrix,
useFixedAnnotationBoundaries, fixedAnnBdrys, numberOfScampIterations,
outputDirectory, verboseOutput, maxSearchTime, debugScampRun,
gaussianScale, randomSeed)
}
\arguments{
  \item{rawDataMatrix}{Call from scamp function.}
  \item{dipT}{Call from scamp function.}
  \item{clusterLB}{Call from scamp function.}
  \item{repeatsAllowed}{Call from scamp function.}
  \item{maxSearchDepth}{Call from scamp function.}
  \item{maxClusterNum}{Call from scamp function.}
  \item{userAnnotations}{Call from scamp function.}
  \item{maxNumberOfGates}{Call from scamp function.}
  \item{randomSearch}{Call from scamp function.}
  \item{randomResidualSearch}{Call from scamp function.}
  \item{finalAnnotationQs}{Call from scamp function.}
  \item{numThreadsRequested}{Call from scamp function.}
  \item{useRestrictedValue}{Call from scamp function.}
  \item{restrictedValueMatrix}{Call from scamp function.}
  \item{useFixedAnnotationBoundaries}{Call from scamp function.}
  \item{fixedAnnBdrys}{Call from scamp function.}
  \item{numberOfScampIterations}{Call from scamp function.}
  \item{outputDirectory}{Call from scamp function.}
  \item{verboseOutput}{Call from scamp function.}
  \item{maxSearchTime}{Call from scamp function.}
  \item{debugScampRun}{Call from scamp function.}
  \item{gaussianScale}{Call from scamp function.}
  \item{randomSeed}{Call from scamp function.}
}
\value{
  This function is exposed by Rcpp. Called using the scamp function.
}
\description{
  This function is exposed by Rcpp. Called using the scamp function.
}
\examples{
clusterMatrix <- as.matrix(iris[,-5])
scampClustering <- scamp(dataSet=clusterMatrix,
                         numberIterations=20,
                         maximumClusterNum=20,
                         numberOfThreads=0) #0 automatically sets to
                                                         #hardware max.
table(scampClustering$RunOffVote)
table(scampClustering$MaxVote)
}