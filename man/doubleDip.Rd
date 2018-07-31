\name{doubleDip}
\alias{doubleDip}
\title{doubleDip}
\usage{
doubleDip(x)
}
\arguments{
\item{x}{A sorted numeric R vector.}
}
\value{
  Returns the 2-dip test p-value.
}
\description{
  Computes the 2-dip test, returns the associated p-value.
}
\examples{
inputDV <- rnorm(300)
doubleDip(sort(inputDV))

inputDV <- c(rnorm(150),rnorm(150,mean=5))
doubleDip(sort(inputDV))

inputDV <- c(rnorm(100),rnorm(100,mean=5),rnorm(100,mean=10))
doubleDip(sort(inputDV))
}

