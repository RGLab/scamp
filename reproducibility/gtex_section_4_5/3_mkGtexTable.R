library(xtable)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gmp)
library(viridis)
source("./helperScripts/clusterValidation.R")
mdf <- readRDS("./sub_100_mergeDF_all_alpha_sel_top50.rds")
clusterMethods <- c("maxCluster01","maxCluster05","maxCluster25","apCluster","mclustCluster","kmeansCluster","kmedoidCluster")
truthMethods <- c("primary.tissue","tissue")
firstTruth <- TRUE
for (truthMethod in truthMethods) {
    firstMethod <- TRUE
    for (method in clusterMethods) {
        clustering <- mdf[,method]
        truth <- mdf[,truthMethod]
        score <- myExternalValidation(as.numeric(as.factor(truth)),
                                      as.numeric(as.factor(clustering)),
                                      method)
        pt <- as.character(truth)
        rc <- as.numeric(as.factor(clustering))
        cn <- sort(unique(rc))
        vc <- rep("NONE",length(rc))
        for (num in cn) {
            lookup <- which(rc==num)
            numPart <- table(pt[lookup])
            maxLook <- which(numPart==max(numPart))
            if (length(maxLook) == 1) {
                vc[lookup] <- names(which(numPart==max(numPart)))
            }
            else if (length(maxLook) == 2) {
                lenLook <- length(lookup)
                subLen1 <- floor(lenLook/2)
                subLook <- sample(lookup,subLen1)
                subLook2 <- setdiff(lookup,subLook)
                firstName <- names(which(numPart==max(numPart))[1])
                secondName <- names(which(numPart==max(numPart))[2])
                if (firstName %in% vc) {
                    vc[lookup] <- secondName
                }
                else {
                    vc[lookup] <- firstName
                }
            }
            else if (length(maxLook) == 3) {
                lenLook <- length(lookup)
                subLen1 <- floor(lenLook/2)
                subLook <- sample(lookup,subLen1)
                subLook2 <- setdiff(lookup,subLook)
                                firstName <- names(which(numPart==max(numPart))[1])
                secondName <- names(which(numPart==max(numPart))[2])
                thirdName <- names(which(numPart==max(numPart))[3])
                if ((firstName %in% vc) && (secondName %in% vc)){
                    vc[lookup] <- thirdName
                }
                else if (firstName %in% vc) {
                    vc[lookup] <- secondName
                }
                else {
                    vc[lookup] <- firstName
                }
            }
            else {
                print("problem")
            }
        }
        score <- rbind(score,length(table(pt)))
        rownames(score)[nrow(score)] <- paste0("num",truthMethod)
        score <- rbind(score,length(table(vc)))
        rownames(score)[nrow(score)] <- paste0(truthMethod,"_ID")
        if (!firstTruth) {
            score <- rbind(score,length(table(clustering)))
        }
        if (firstMethod) {
            scoreMat <- score
            firstMethod <- FALSE
        }

        else scoreMat <- cbind(scoreMat,score)
    }
    scoreMat <- cbind(scoreMat,truthMethod)
    if (firstTruth) {
        outMat <- scoreMat
        firstTruth <- FALSE
    }
    else {
        outMat <- rbind(outMat,scoreMat)
    }
}

flookup <- sort(c(which(rownames(outMat)=="VarInfo"),
                  which(rownames(outMat)=="Adj Rand Index"),
                  which(rownames(outMat)=="numprimary.tissue"),
                  nrow(outMat),
                  which(rownames(outMat)=="numtissue")))
fmat <- outMat[flookup,]
fmat <- fmat[,-ncol(fmat)]
colnames(fmat) <- c("scamp01","scamp05","scamp25","ap","mclust","kmeans","kmedoid")
rownames(fmat) <- c("Primary Tissue: VI",
                    "Primary Tissue: ARI",
                    "# Primary Tissue Labels",
                    "Tissue: VI",
                    "Tissue: ARI",
                    "# Tissue Labels",
                    "# of Clusters by Method")     
xtable(fmat)
