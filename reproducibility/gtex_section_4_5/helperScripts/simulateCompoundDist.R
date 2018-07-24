Posdef <- function (n, ev = runif(n, 0, 3))
{
    #function written by Ravi Varadhan, from r-help mailing list
    #Thu Feb 7, 20:02:30 CET 2008
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
}

genMvtSample <- function(numRows,sampleDim,tDOF=5){
    mu <- rep(2.5,sampleDim)
    sig <- Posdef(sampleDim)
    x <- rmvt(numRows,delta=mu,sigma=sig,df=tDOF)
    return(x)
}

absSqrt <- function(x) {return(sqrt(abs(x)))}
abs23 <- function(x) {return((abs(x))^(2/3))}
square <- function(x) {return((x^2))}
abs32 <- function(x) {return((abs(x))^(3/2))}
abs54 <- function(x) {return((abs(x))^(5/4))}
expb4 <- function(x) {return(exp((x/4)))}
logAbs <- function(x) {return(log(1+abs(x)))}

genSimulationSample <- function(sampleDim,sampleVec,
                                tRegime="none",sRegime="gaussianOnly",
                                fixedMeanMatrix=NA) {
    outData <- matrix(nrow=0,ncol=sampleDim)
    for (sampleNum in seq(length(sampleVec))) {
        currentMu <- as.numeric(fixedMeanMatrix[sampleNum,])
        print(currentMu)
        currentSigma<- Posdef(sampleDim)
        if (sRegime == "tPlusGauss") {
            if ((sampleNum %% 2) == 0) {
                currentSample <- mvrnorm(sampleVec[sampleNum],mu=currentMu,Sigma=currentSigma)
            }
            else {
                currentSample <- rmvt(sampleVec[sampleNum],delta=currentMu,sigma=currentSigma,df=5)
            }
        }
        else  {
            currentSample <- mvrnorm(sampleVec[sampleNum],mu=currentMu,Sigma=currentSigma)
        }
        if (tRegime == "four") {
            if ((sampleNum %% 4) == 0) {
                outData <- rbind(outData,absSqrt(currentSample))
            }
            else if ((sampleNum %% 4) == 2) {
                outData <- rbind(outData,exp(currentSample))
            }
            else if ((sampleNum %% 4) == 3) {
                outData <- rbind(outData,square(currentSample))
            }
            else {
                outData <- rbind(outData,currentSample)
            }
        }
        else if (tRegime == "none") {
            outData <- rbind(outData,currentSample)
        }
    }
    return(outData)
}
