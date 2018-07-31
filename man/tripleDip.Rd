\name{tripleDip}
\alias{tripleDip}
\title{tripleDip}
\usage{
tripleDip(x)
}
\arguments{
\item{x}{A sorted numeric R vector.}
}
\value{
  Returns the 3-dip test p-value.
}
\description{
  Computes the 3-dip test, returns the associated p-value.
}
\examples{
inputDV <- rnorm(300)
tripleDip(sort(inputDV))

inputDV <- c(rnorm(150),rnorm(150,mean=5))
tripleDip(sort(inputDV))

inputDV <- c(rnorm(100),rnorm(100,mean=5),rnorm(100,mean=10))
tripleDip(sort(inputDV))
}

