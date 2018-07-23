suppressMessages(library(optparse))
suppressMessages(library(tmvtnorm))
suppressMessages(library(mvtnorm))
suppressMessages(library(MASS))
suppressMessages(library(cluster))
suppressMessages(library(mclust))
suppressMessages(library(mcclust))
suppressMessages(library(scamp))
suppressMessages(library(gmp))
suppressMessages(library(apcluster))
source("./helperScripts/clusterValidation.R")
source("./helperScripts/simulateCompoundDist.R")
option_list = list(make_option(c("-v", "--verbose"), action="store", default=FALSE,help="Should the program print extra stuff out? [default %default]"),
                   make_option(c("-a", "--regime_number"), action="store", default=NA, type="character", help="Subject to process"),
                   make_option(c("-b", "--amb_number"), action="store", default=NA, type="character", help="Subject to process"),
                   make_option(c("-c", "--dim_number"), action="store", default=NA, type="character", help="Subject to process"),
                   make_option(c("-d", "--trans_number"), action="store", default=NA, type="character", help="Subject to process"))

opt = parse_args(OptionParser(option_list=option_list))

#simulation settings
numIter <- 10
#simulation scenarios are controlled by the following three parameters.
#a fixed setting will produce half the jobs in a simulation scenario
#we describe the scenario as follows: {ssChoice, dsChoice, mvChoice}
#Scenario 1: {1,1,1} and {1,1,2}
#Scenario 2: {2,1,1} and {2,1,2}
#Secnario 3: {2,2,1} and {2,2,2}
ssChoice <- 1 #manually set to 1 or 2 -- sampleSize 3,000 == 1, sampleSize 30,000 ==2 
dsChoice <- 1 #manually set to 1 or 2 -- if ssChoice == 2, dsChoice == 1 -> numCluster = 30, dsChoice == 2 -> numCluster = 56
mvChoice <- 1 #manually set to 1 or 2 -- adjoin mean vector with {0,6} if 1, {0,3,6} if 2.

#remaining settings determined from command line
regimeChoice <- as.numeric(opt$regime_number)
transformationChoice <- as.numeric(opt$trans_number)
dimensionChoice <- as.numeric(opt$dim_number)
ambientChoice <- as.numeric(opt$amb_number)

#all simulation parameters
simDimension <- 20
sRegimes <- c("tPlusGauss","gaussianOnly")
transformationTypes <- c("none","four")
ambientDimensions <- c(0,20)
testDimensions <- c(0,20)
useBigSample <- c(FALSE,TRUE)
bigSampleDiverse <- c(FALSE,TRUE)
useTwoMV <- c(FALSE,TRUE)

#setup job
ambientDimension <- ambientDimensions[ambientChoice]
testDimension <- testDimensions[dimensionChoice]
transformationType <- transformationTypes[transformationChoice]
sRegime <- sRegimes[regimeChoice]
bigSample <- as.numeric(useBigSample[ssChoice])
hetBigSample <- as.numeric(bigSampleDiverse[dsChoice])
twoMV <- as.numeric(useTwoMV[mvChoice])

jobDesc <- paste0("aDim_",ambientDimension,
                  "_cDim_",testDimension,
                  "_tType_",transformationType,
                  "_sRegime_",sRegime)

#log job specs
print(paste0("Running job ",jobDesc))
remainingIter <- numIter

#set seed to reproduce
set.seed(91827)

#determine potential mean vectors
mv1Dimension <- 20
mvl <- list()
for (i in seq(mv1Dimension)) {
    mvl <- append(mvl,list(c(6,0))) 
}
allMeanVectors <- expand.grid(mvl)

mv2Dimension <- 10
mvl2 <- list()
for (i in seq(mv2Dimension)) {
    mvl2 <- append(mvl2,list(c(6,3,0)))
}
allMeanVectors2 <- expand.grid(mvl2)
allMeanVectors2 <- cbind(allMeanVectors2,allMeanVectors2)

#cluster specification
sampleSizes1  <- c(1000,500,250,250,250,250,125,125,125,125)
if (bigSample > 0) {
    if (hetBigSample) {
        sampleSizes1  <- c(10000,5000,2500,2500,1250,1250,1250,1250,500,500,500,500,
                           250,250,250,250,250,250,250,250,
                           100,100,100,100,100,100,100,100,100,100)
    }
    else {
        sampleSizes1  <- c(10000,5000,2500,2500,2500,2500,1250,1250,1250,1250)
    }
}
sampleSizes2  <- sampleSizes1
sampleSizes2 <- sampleSizes2[sample(seq(length(sampleSizes1)))]
getTV <- function(sampleSizes) {
    currentClusterNum <- 1
    truthVec <- c()
    for (i in seq(length(sampleSizes))) {
        currentClusterLength <- sampleSizes[i]
        truthVec <- append(truthVec,rep(currentClusterNum,currentClusterLength))
        currentClusterNum <- currentClusterNum + 1
    }
    return(truthVec)
}
tv1 <- as.character(getTV(sampleSizes1))
tv2 <- as.character(getTV(sampleSizes2))
truthVec <- as.numeric(as.factor(tv1))
if (testDimension > 0) {
    truthVec <- as.numeric(as.factor(paste0(tv1,tv2)))
}
print(table(truthVec))
while (remainingIter > 0) {
    print(paste0("Remaining iterations: ",remainingIter))
    mvChoices1 <- sample(seq(nrow(allMeanVectors)),length(sampleSizes1))
    meanVectors1 <- allMeanVectors[mvChoices1,]
    if (twoMV) {
        mvChoices2 <- sample(seq(nrow(allMeanVectors2)),length(sampleSizes1))
        meanVectors2 <- allMeanVectors2[mvChoices2,]
    }
    else {
        mvChoices2 <- sample(seq(nrow(allMeanVectors)),length(sampleSizes1))
        meanVectors2 <- allMeanVectors[mvChoices2,]
    }
    ssList <- list(sampleSizes1,sampleSizes2)
    mvList <- list(meanVectors1,meanVectors2)
    staticData <- genSimulationSample(sampleDim=simDimension,
                                      sampleVec=ssList[[1]],
                                      tRegime=transformationType,
                                      sRegime=sRegime,
                                      fixedMeanMatrix=mvList[[1]])
    if (testDimension > 0) {
        newData <- genSimulationSample(sampleDim=simDimension,
                                       sampleVec=ssList[[2]],
                                       tRegime=transformationType,
                                       sRegime=sRegime,
                                       fixedMeanMatrix=mvList[[2]])
        staticData <- cbind(staticData,newData)
    }
    if (ambientDimension > 0) {
        dp2 <- genMvtSample(nrow(staticData),ambientDimension,tDOF=5)
        staticData <- cbind(staticData,dp2)
    }
    print(paste0("Static data cols: ",ncol(staticData)))
    centeredData <- scale(staticData,center=TRUE,scale=FALSE)
    scaledData <- scale(staticData,center=TRUE,scale=TRUE)

    print("starting scamp")
    colnames(staticData) <- paste0("V",seq(ncol(staticData)))
    startTime <- proc.time()
    scampLabels <- scamp(dataSet=staticData,
                         numberIterations=1,
                         pValueThreshold=0.25,
                         clusterOutputString=paste0("./interScamp/",jobDesc),
                         randomCandidateSearch=TRUE,
                         randomResidualCandidateSearch=TRUE,
                         maximumClusterNum=(50*ncol(staticData)),
                         finalAnnQs = c(0.49, 0.5, 0.51),
                         getDebugInfo=FALSE,
                         numberOfThreads=16,
                         randomSeed=72313)

    scampTime <- proc.time()-startTime
    scampClustering <- scampLabels[[2]]
    resultsScamp <- myExternalValidation(truthVec,as.numeric(as.factor(scampClustering)),"SCAMP")
    
    resultsScamp <- rbind(resultsScamp,length(unique(scampClustering)))
    rownames(resultsScamp)[length(rownames(resultsScamp))] <- "Num Clusters"
    
    resultsScamp <- rbind(resultsScamp,as.numeric(scampTime[3]))
    rownames(resultsScamp)[length(rownames(resultsScamp))] <- "Elapsed Time"

    resultsScamp <- rbind(resultsScamp,length(which(table(scampClustering) > 24)))
    rownames(resultsScamp)[length(rownames(resultsScamp))] <- "Num Clusters > 25"
    print("ending scamp")

    print("starting scamp 20")
    colnames(staticData) <- paste0("V",seq(ncol(staticData)))
    startTime <- proc.time()
    scampLabels20 <- scamp(dataSet=staticData,
                         numberIterations=20,
                         pValueThreshold=0.25,
                         clusterOutputString=paste0("./interScamp/",jobDesc),
                         randomCandidateSearch=TRUE,
                         randomResidualCandidateSearch=TRUE,
                         maximumClusterNum=(50*ncol(staticData)),
                         finalAnnQs = c(0.49, 0.5, 0.51),
                         getDebugInfo=FALSE,
                         numberOfThreads=16,
                         randomSeed=23134)
    scampTime20 <- proc.time()-startTime
    scampClustering20 <- scampLabels20[[2]]
    
    resultsScamp20 <- myExternalValidation(truthVec,as.numeric(as.factor(scampClustering20)),"SCAMP_20")

    resultsScamp20 <- rbind(resultsScamp20,length(unique(scampClustering20)))
    rownames(resultsScamp20)[length(rownames(resultsScamp20))] <- "Num Clusters"

    resultsScamp20 <- rbind(resultsScamp20,as.numeric(scampTime20[3]))
    rownames(resultsScamp20)[length(rownames(resultsScamp20))] <- "Elapsed Time"

    resultsScamp20 <- rbind(resultsScamp20,length(which(table(scampClustering20) > 24)))
    rownames(resultsScamp20)[length(rownames(resultsScamp20))] <- "Num Clusters > 25"

    print("ending scamp 20")

    print("starting model based clustering")
    start.time <- proc.time()
    if (bigSample) {
        initVec <- sample(seq(nrow(staticData)),floor(0.1*nrow(staticData)))
        mres <- Mclust(staticData,
                       G=length(table(truthVec)),
                       initialization=list(subset=initVec))
    }
    else {
        mres <- Mclust(staticData,G=length(table(truthVec)))
    }
    mclustTime <- proc.time()-start.time
    mresLabels <- mres$classification
    resultsMclust <- myExternalValidation(truthVec,mresLabels,"Oracle_Mclust")
    resultsMclust <- rbind(resultsMclust,length(unique(mresLabels)))
    rownames(resultsMclust)[length(rownames(resultsMclust))] <- "Num Clusters"
    resultsMclust <- rbind(resultsMclust,as.numeric(mclustTime[3]))
    rownames(resultsMclust)[length(rownames(resultsMclust))] <- "Elapsed Time"
    resultsMclust <- rbind(resultsMclust,length(which(table(mresLabels) > 24)))
    rownames(resultsMclust)[length(rownames(resultsMclust))] <- "Num Clusters > 25"

    print("ending model based clustering")

    print("starting oracle kmeans")
    startTime <- proc.time()
    fit <- kmeans(scaledData,length(table(truthVec)),iter.max=10000,nstart=100)
    kmLabels <- fit$cluster
    kmeansTime <- proc.time()-startTime
    resultsKmeans <- myExternalValidation(truthVec,kmLabels,"Oracle_K_means")

    resultsKmeans <- rbind(resultsKmeans,length(table(truthVec)))
    rownames(resultsKmeans)[length(rownames(resultsKmeans))] <- "Num Clusters"
    
    resultsKmeans <- rbind(resultsKmeans,as.numeric(kmeansTime[3]))
    rownames(resultsKmeans)[length(rownames(resultsKmeans))] <- "Elapsed Time"

    resultsKmeans <- rbind(resultsKmeans,length(table(truthVec)))
    rownames(resultsKmeans)[length(rownames(resultsKmeans))] <- "Num Clusters > 25"
    
    print("ending oracle kmeans")

    if (bigSample == 0) {
        print("starting medoid clustering")
        startTime <- proc.time()
        orcmedFig <- pam(scaledData,k=length(table(truthVec)))
        kmedTime <- proc.time()-startTime
        orcmedLabels <- orcmedFig$clustering
        resultsKmed <- myExternalValidation(truthVec,orcmedLabels,"Oracle_K_Medoid")
        resultsKmed <- rbind(resultsKmed,length(table(truthVec)))
        rownames(resultsKmed)[length(rownames(resultsKmed))] <- "Num Clusters"
        resultsKmed <- rbind(resultsKmed,as.numeric(kmedTime[3]))
        rownames(resultsKmed)[length(rownames(resultsKmed))] <- "Elapsed Time"
        resultsKmed <- rbind(resultsKmed,length(table(truthVec)))
        rownames(resultsKmed)[length(rownames(resultsKmed))] <- "Num Clusters > 25"
        print("ending medoid clustering")
    }
    else {
        print("starting clara clustering")
        startTime <- proc.time()
        claraOracle <- clara(x=scaledData,
                             k=length(table(truthVec)),
                             samples=50,
                             sampsize=3000)
        claraTime <- proc.time()-startTime
        claraLabels <- claraOracle$clustering
        resultsClara <- myExternalValidation(truthVec,claraLabels,"Oracle_Clara")
        resultsClara <- rbind(resultsClara,length(table(truthVec)))
        rownames(resultsClara)[length(rownames(resultsClara))] <- "Num Clusters"
        resultsClara <- rbind(resultsClara,as.numeric(claraTime[3]))
        rownames(resultsClara)[length(rownames(resultsClara))] <- "Elapsed Time"
        resultsClara <- rbind(resultsClara,length(table(truthVec)))
        rownames(resultsClara)[length(rownames(resultsClara))] <- "Num Clusters > 25"
        print("ending clara clustering")
    }
    
    print("starting affinity propagation")
    startTime <- proc.time()
    if (bigSample) {
        print("ap is leveraged!")
        apResult <- apclusterL(negDistMat(r=2),
                               scaledData,
                               maxits=2000,
                               convits=2000,
                               frac=0.1,
                               sweeps=5,
                               q=0,
                               seed=76434)
    }
    else {
        print("full AP")
        apResult <- apclusterK(negDistMat(r=2),
                               scaledData,
                               bimaxit=100,
                               maxits=2000,
                               convits=2000,
                               K=length(table(truthVec)),
                               seed=541709)
    }
    apTime <- proc.time()-startTime
    apresLabels <- rep(0,nrow(staticData))
    apclusnum <- 1
    for (i in seq(length(apResult@clusters))) {
        apresLabels[as.numeric(apResult@clusters[[i]])] <- apclusnum
        apclusnum <- apclusnum + 1
    }
    if (bigSample) {
        resultsApres <- myExternalValidation(truthVec,apresLabels,"Leveraged_AP")
    }
    else {
        resultsApres <- myExternalValidation(truthVec,apresLabels,"Oracle_AP")
    }
    
    resultsApres <- rbind(resultsApres,length(unique(apresLabels)))
    rownames(resultsApres)[length(rownames(resultsApres))] <- "Num Clusters"
    resultsApres <- rbind(resultsApres,as.numeric(apTime[3]))
    rownames(resultsApres)[length(rownames(resultsApres))] <- "Elapsed Time"
    resultsApres <- rbind(resultsApres,length(which(table(mresLabels) > 24)))
    rownames(resultsApres)[length(rownames(resultsApres))] <- "Num Clusters > 25"

    print("affinity propagation complete")

    if (bigSample == 0) {
        omat <- cbind(resultsScamp,
                      resultsScamp20,
                      resultsKmeans,
                      resultsMclust,
                      resultsKmed,
                      resultsApres)
    }
    else {
        omat <- cbind(resultsScamp,
                      resultsScamp20,
                      resultsKmeans,
                      resultsMclust,
                      resultsClara,
                      resultsApres)
    }
    if (remainingIter == numIter) {
        resMat <- omat
    }
    else {
        resMat <- resMat + omat
    }
    remainingIter <- remainingIter-1
}

finalRes <- resMat/numIter

if (!dir.exists(paste0("./sepMix_bigSample_",bigSample,"_bigHet_",hetBigSample,"_twoMV_",twoMV))) {
    dir.create(paste0("./sepMix_bigSample_",bigSample,"_bigHet_",hetBigSample,"_twoMV_",twoMV))
}

if (!dir.exists(paste0("./sepMix_bigSample_",bigSample,"_bigHet_",hetBigSample,"_twoMV_",twoMV,"/",jobDesc))) {
    dir.create(paste0("./sepMix_bigSample_",bigSample,"_bigHet_",hetBigSample,"_twoMV_",twoMV,"/",jobDesc))
}

saveRDS(finalRes,paste0("./sepMix_bigSample_",bigSample,"_bigHet_",hetBigSample,"_twoMV_",twoMV,"/",jobDesc,"/jobOutput.rds"))
