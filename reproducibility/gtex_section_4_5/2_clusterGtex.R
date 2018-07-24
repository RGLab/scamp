library(irlba)
library(scamp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gmp)
suppressMessages(library(mclust))
suppressMessages(library(apcluster,lib="/home/egreene/local/R_libs"))
suppressMessages(library(cluster))
source("./helperScripts/clusterValidation.R")
geneData <- readRDS("./geneData.rds")
countData <- readRDS("./allCounts.rds")
dim(countData)
zeroCounts <- apply(countData,1,function(x){length(which(x==0))})
subCount <- countData[-which(zeroCounts>8554),]
dim(subCount)
l2SubCount <- apply(subCount,2,function(x){log2(1+x)})
l2NormCount <- apply(l2SubCount,2,function(x){x/sum(x)})
transposeData <- t(l2NormCount)
tsneData <- readRDS("./tsneData.rds")
scaleMatrix <- scale(transposeData,center=TRUE,scale=FALSE)
SVs <- irlba(scaleMatrix, nv=100)
svMatrixAll <- scaleMatrix %*% SVs$v
colnames(svMatrixAll) <- paste0("SV",seq(ncol(svMatrixAll)))
saveRDS(svMatrixAll,"./svMatrixAll.rds")
svMatrix <- svMatrixAll[,1:50]
saveRDS(svMatrix,"./svMatrix50.rds")
print(dim(svMatrix))
print("scamp 25 start")
scampLabels25 <- scamp(dataSet=svMatrix,
                       numberIterations=100,
                       pValueThreshold=0.25,
                       randomCandidateSearch=TRUE,
                       randomResidualCandidateSearch=TRUE,
                       maximumClusterNum=(100*ncol(svMatrix)),
                       finalAnnQs = c(0.495, 0.5, 0.505),
                       getDebugInfo=FALSE,
                       numberOfThreads=28,
                       randomSeed=72313)

maxLabels25 <- scampLabels25[[2]]
maxLabels25 <- as.numeric(as.factor(maxLabels25))
length(table(maxLabels25))

print("scamp 05 start")
scampLabels05 <- scamp(dataSet=svMatrix,
                       numberIterations=100,
                       pValueThreshold=0.05,
                       randomCandidateSearch=TRUE,
                       randomResidualCandidateSearch=TRUE,
                       maximumClusterNum=(100*ncol(svMatrix)),
                       finalAnnQs = c(0.495, 0.5, 0.505),
                       getDebugInfo=FALSE,
                       numberOfThreads=28,
                       randomSeed=72313)

maxLabels05 <- scampLabels05[[2]]
maxLabels05 <- as.numeric(as.factor(maxLabels05))
length(table(maxLabels05))


print("scamp 01 start")
scampLabels01 <- scamp(dataSet=svMatrix,
                       numberIterations=100,
                       pValueThreshold=0.01,
                       randomCandidateSearch=TRUE,
                       randomResidualCandidateSearch=TRUE,
                       maximumClusterNum=(100*ncol(svMatrix)),
                       finalAnnQs = c(0.495, 0.5, 0.505),
                       getDebugInfo=FALSE,
                       numberOfThreads=28,
                       randomSeed=72313)

maxLabels01 <- scampLabels01[[2]]
maxLabels01 <- as.numeric(as.factor(maxLabels01))
length(table(maxLabels01))

set.seed(28397)
print("apcluster start")
apResult <- apcluster(negDistMat(r=2),
                      svMatrix,
                      maxits=5000,
                      convits=5000,
                      lam=0.95,
                      q=0,
                      seed=71231)
                      
apresLabels <- rep(0,nrow(svMatrix))
apclusnum <- 1
for (i in seq(length(apResult@clusters))) {
    apresLabels[as.numeric(apResult@clusters[[i]])] <- apclusnum
    apclusnum <- apclusnum + 1
}

set.seed(28397)
print("mclust start")
mres <- Mclust(svMatrix,
               G=2:100)
mresLabels <- mres$classification

set.seed(28397)
print("kmeans start")
fit <- kmeans(svMatrix,length(table(mresLabels)),iter.max=10000,nstart=100)
kmLabels <- fit$cluster

set.seed(28397)
print("kmed start")
orcmedFig <- pam(svMatrix,k=length(table(mresLabels)))
orcmedLabels <- orcmedFig$clustering
sinfo <- colnames(geneData)[-c(1,2)]

scampDF <- data.frame(SAMPID=sinfo,
                      maxCluster01=maxLabels01,
                      maxCluster05=maxLabels05,
                      maxCluster25=maxLabels25,
                      apCluster=apresLabels,
                      mclustCluster=mresLabels,
                      kmeansCluster=kmLabels,
                      kmedoidCluster=orcmedLabels,
                      tsneX=tsneData$Y[,1],
                      tsneY=tsneData$Y[,2],
                      stringsAsFactors=FALSE)

sampleData <- readRDS("./samplesData.rds")
sampleData$SAMPID <- as.character(as.factor(sampleData$SAMPID))
mergeDF <- inner_join(scampDF,sampleData,by=c("SAMPID"))
saveRDS(mergeDF,"./sub_100_mergeDF_all_alpha_sel_top50.rds")
