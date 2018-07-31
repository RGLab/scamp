\name{singleDip}
\alias{singleDip}
\title{singleDip}
\usage{
singleDip(x)
}
\arguments{
\item{x}{A sorted numeric R vector.}
}
\value{
  Returns the dip test p-value.
}
\description{
  Computes the dip test of Hartigan and Hartigan,
  returns the associated p-value.
}
\examples{
inputDV <- rnorm(300)
singleDip(sort(inputDV))

inputDV <- c(rnorm(150),rnorm(150,mean=5))
singleDip(sort(inputDV))

inputDV <- c(rnorm(100),rnorm(100,mean=5),rnorm(100,mean=10))
singleDip(sort(inputDV))
}
