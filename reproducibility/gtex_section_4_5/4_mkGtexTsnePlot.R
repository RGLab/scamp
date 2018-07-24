library(ggplot2)
library(ggalt)
library(cowplot)
library(dplyr)
library(gmp)
library(viridis)
source("./helperScripts/clusterValidation.R")
mdf <- readRDS("./sub_100_mergeDF_all_alpha_sel_top50.rds")
clusterMethods <- c("maxCluster01","maxCluster05","maxCluster25","apCluster","mclustCluster","kmeansCluster","kmedoidCluster")
truthMethod <- "tissue"
truth <- as.character(mdf[,truthMethod])
truthVals <- sort(unique(truth))
baseColors <- viridis(length(unique(truth)))
set.seed(13457)
baseColors <- baseColors[sample(seq(length(baseColors)))]
colorMap <- as.character(rep("#FFFFFFFF",nrow(mdf)))
for (colorNum in seq(length(truthVals))) {
    truthVal <- truthVals[colorNum]
    lookup <- which(truth==truthVal)
    colorMap[lookup] <- baseColors[colorNum]
}
residualColors <- magma(100)
colorMapList <- list()
for (method in clusterMethods) {
    clustering <- mdf[,method]
    pt <- truth
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
        else if (length(maxLook) == 4) {
            firstName <- names(which(numPart==max(numPart))[1])
            secondName <- names(which(numPart==max(numPart))[2])
            thirdName <- names(which(numPart==max(numPart))[3])
            fourthName <- names(which(numPart==max(numPart))[4])
            print(firstName)
            print(secondName)
            print(thirdName)
            print(fourthName)
            
            if ((firstName %in% vc) && (secondName %in% vc) && (thirdName %in% vc)) {
                vc[lookup] <- fourthName
            }
            else if ((firstName %in% vc) && (secondName %in% vc)){
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
            print(length(maxLook))
        }
    }
    uniqVals <- sort(unique(vc))
    uniqNums <- sort(unique(rc))
    coloredNums <- c()
    methodColorMap <- as.character(rep("#FFFFFFFF",nrow(mdf)))
    for (uval in uniqVals) {
        ulook <- which(vc == uval)
        usummary <- table(rc[ulook])
        maxName <- as.numeric(names(which(usummary==max(usummary)))[1])
        maxColor <- unique(colorMap[which(truth==uval)])
        if (length(maxColor) > 1) {
            print("Problem")
        }
        methodColorMap[which(rc==maxName)] <- maxColor
        if (maxName %in% coloredNums) {
            print(paste0(maxName,"colored twice."))
        }
        coloredNums <- append(coloredNums,maxName)
    }
    residNums <- sort(setdiff(uniqNums,coloredNums))
    if (length(residNums)) {
        colorStep <- floor(length(residualColors)/(length(residNums)+1))
        currentColor <- colorStep
        for (rcolnum in residNums) {
            methodColorMap[which(rc==rcolnum)] <- residualColors[currentColor]
            currentColor <- (currentColor + colorStep)
        }
    }
    colorMapList <- append(colorMapList,list(methodColorMap))
    names(colorMapList)[length(colorMapList)] <- method
}

p <- ggplot(mdf,aes(x=tsneX,y=tsneY)) +
    geom_point(color=colorMap) +
    theme_bw()+
   theme(legend.position="none") + ggtitle("Color by tissue label")
library(ggalt)
classNames <- c("Brain - Anterior cingulate cortex (BA24)", 
                "Brain - Caudate (basal ganglia)",          
                "Brain - Cerebellar Hemisphere",            
                "Cells - Transformed fibroblasts",
                "Skin - Not Sun Exposed (Suprapubic)",
                "Skin - Sun Exposed (Lower leg)")
vcp <- c("#0066ff","#0066ff",
         "#0066ff",
         "#ff66cc","#ff66cc","#ff66cc")
for (nameNum in seq(length(classNames))) {
    name <- classNames[nameNum]
    dfSub <- mdf[which(mdf$tissue==name),c("tsneX","tsneY")]
    p <- p + geom_encircle(data=dfSub,aes(x=tsneX,y=tsneY),colour=vcp[nameNum])
}

p1a <- ggplot(mdf,aes(x=tsneX,y=tsneY)) +
    geom_point(color=colorMapList[["maxCluster01"]]) +
    theme_bw()+
        theme(legend.position="none") +
            labs(title=expression("SCAMP "~alpha~" = 0.01, ARI = 0.649"))
p1b <- ggplot(mdf,aes(x=tsneX,y=tsneY)) +
    geom_point(color=colorMapList[["maxCluster05"]]) +
    theme_bw()+
        theme(legend.position="none") +
            labs(title=expression("SCAMP "~alpha~" = 0.05, ARI = 0.657"))
p1c <- ggplot(mdf,aes(x=tsneX,y=tsneY)) +
    geom_point(color=colorMapList[["maxCluster25"]]) +
    theme_bw()+
        theme(legend.position="none") +
                        labs(title=expression("SCAMP "~alpha~" = 0.25, ARI = 0.654"))
p2 <- ggplot(mdf,aes(x=tsneX,y=tsneY)) +
    geom_point(color=colorMapList[["apCluster"]]) +
    theme_bw()+
    theme(legend.position="none") + ggtitle("AP, ARI: 0.690")
p3 <- ggplot(mdf,aes(x=tsneX,y=tsneY)) +
    geom_point(color=colorMapList[["mclustCluster"]]) +
    theme_bw()+
    theme(legend.position="none") + ggtitle("Mclust, ARI: 0.606")
p4 <- ggplot(mdf,aes(x=tsneX,y=tsneY)) +
    geom_point(color=colorMapList[["kmeansCluster"]]) +
    theme_bw()+
    theme(legend.position="none") + ggtitle("K-means, ARI: 0.557")
p5 <- ggplot(mdf,aes(x=tsneX,y=tsneY)) +
    geom_point(color=colorMapList[["kmedoidCluster"]]) +
    theme_bw()+
    theme(legend.position="none") + ggtitle("K-medoid, ARI: 0.564")
pj <- plot_grid(p,p1a,
                p1b,p1c,
                p2,p3,
                p4,p5,
                nrow=3,ncol=3,labels=c("1","2","3","4","5","6","7","8"))
pj
