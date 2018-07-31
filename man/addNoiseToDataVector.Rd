\name{addNoiseToDataVector}
\alias{addNoiseToDataVector}
\title{addNoiseToDataVector}
\usage{
addNoiseToDataVector(dataVector, scaleGaussian, rSeed)
}
\arguments{
\item{dataVector}{A numeric R vector.}

\item{scaleGaussian}{Gaussian variance scaling parameter gamma.}

\item{rSeed}{Random seed.}
}
\value{
  Returns a numeric R vector with sturctured noise added.
}
\description{
  This function is used to add phase 1 structured noise to a data vector.
}
\examples{
inputDV <- rpois(1000,lambda=5)
noiseDV <- addNoiseToDataVector(inputDV,4,123)
par(mfrow=c(2,1))
plot(ecdf(inputDV),main="Input ECDF")
plot(ecdf(noiseDV),main="Noise ECDF")
}
