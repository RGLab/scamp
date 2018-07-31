\name{tsGates}
\alias{tsGates}
\title{tsGates}
\usage{
tsGates(xVec,modePrior)
}
\arguments{
  \item{xVec}{A sorted numeric R vector.}

  \item{modePrior}{A non-negative integer indicating if the number of
    modal groups is known. Set to 0 if there is no prior estimate.}
}
\value{
  Returns vector of modal-group boundaries.
  First entry is always minimum value of xVec.
  Last entry is always maximum value of xVec.
  Intermediate values mark modal sub-collections.
}
\description{
  For xVec with more than 400 entries,
  uses the taut string to find modal sub-collections.
}
\examples{
inputDV <- c(rnorm(150),rnorm(150,mean=5))
tsGates(sort(inputDV),0)

inputDV <- c(rnorm(100),rnorm(100,mean=5),rnorm(100,mean=10))
tsGates(sort(inputDV),0)
}

