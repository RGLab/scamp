library(scamp)
library(ggalt)
library(classifly)
library(ClusterR)
library(Rtsne)
library(ggplot2)
library(cowplot)
library(viridis)
#helper function for coloring t-SNE maps
mapToUnitInterval <- function(xIn) {
    trimVals <- as.numeric(quantile(xIn,probs=c(0.025,0.975)))
    lowLookup <- which(xIn <= trimVals[1])
    highLookup <- which(xIn >= trimVals[2])
    x <- xIn
    x[lowLookup] <- trimVals[1]
    x[highLookup] <- trimVals[2]
    y <- (x-min(x))
    return((2*(y/max(y))-1))
}

#cluster data.
scampData <- olives[,-c(1,2)]
scampMatrix <- as.matrix(scampData)

#note: if number of threads changed, results may not match paper exactly.
#however, the results should be very similar.
#this is because the random candidate cluster search allows all launched
#threads to finish after the maximumClusterNum threshold has been crossed.
#each launched thread increments the random seed, and so
#modifying the threads parameter will affect all subsequent scamp iterations.
scampLabels <- scamp(dataSet=scampMatrix,
                     numberIterations=5000,
                     pValueThreshold=0.30,
                     clusterOutputString = c("./scampRnd5000"),
                     randomCandidateSearch=TRUE,
                     randomResidualCandidateSearch=TRUE,
                     maximumClusterNum=400,
                     numberOfThreads=28,
                     finalAnnQs = c(0.499, 0.5, 0.501),
                     gaussianScaleParameter=4,
                     getDebugInfo=FALSE,
                     randomSeed=8327)
#score clustering
maxClustering <- scampLabels[[2]]
length(table(maxClustering))
length(which(table(maxClustering) > 25))
length(table(maxClustering))
external_validation(as.numeric(as.factor(maxClustering)),
                    as.numeric(as.factor(olives$Area)),
                    summary_stats=TRUE)

#get t-SNE map
set.seed(239478)
tsneData <- Rtsne(scampMatrix,perplexity=30,theta=0,pca_scale=T)
unitData <- apply(as.matrix(scampMatrix),2,mapToUnitInterval)
tsneOut <- data.frame(tsneX=tsneData$Y[,1],
                      tsneY=tsneData$Y[,2],
                      Area=as.factor(olives$Area),
                      ScampLabels=as.factor(maxClustering),
                      Palmitic=unitData[,"palmitic"],
                      Palmitoleic=unitData[,"palmitoleic"],
                      Stearic=unitData[,"stearic"],
                      Oleic=unitData[,"oleic"],
                      Linoleic=unitData[,"linoleic"],
                      Linolenic=unitData[,"linolenic"],
                      Arachidic=unitData[,"arachidic"],
                      Eicosenoic=unitData[,"eicosenoic"])
#generate plots
set.seed(9498)
areaPlotColors <- viridis(length(unique(olives$Area)))
areaPlotColors <- areaPlotColors[sample(seq(length(areaPlotColors)))]
areaShapeNumbers <- seq(nlevels(tsneOut$Area))
scampShapeNumbers <- seq(nlevels(tsneOut$ScampLabels))
plotColors <- magma(length(unique(maxClustering)))
plotColors <- plotColors[sample(seq(length(plotColors)))]
p <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=Area,shape=Area))+
    geom_point()+
    scale_shape_manual(values=areaShapeNumbers) +
    scale_color_manual(values=areaPlotColors)+
    theme_bw()+
    theme(legend.position="none")+
    ggtitle("Italy Area")

areas <- c("South-Apulia","North-Apulia","Sicily","Calabria")
vcp <- c("#0066ff","#ff6600","#ff66ff","#00aaaa")
for (nameNum in seq(length(areas))) {
    name <- areas[nameNum]
    print(name)
    dfSub <- tsneOut[which(tsneOut$Area==name),c("tsneX","tsneY")]
    chSub <- dfSub[chull(dfSub),]
    colnames(chSub) <- c("subTX","subTY")
    colnames(dfSub) <- c("subTX","subTY")
    p <- p + geom_encircle(data=dfSub,
                           aes(x=subTX,subTY),
                           colour=vcp[nameNum],
                           linetype=(4*nameNum-1),
                           inherit.aes = FALSE)
}
p2 <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=ScampLabels,shape=ScampLabels))+
    geom_point()+
    scale_shape_manual(values=scampShapeNumbers) +
    scale_color_manual(values=plotColors)+
    theme_bw()+
    theme(legend.position="none")+
    ggtitle("SCAMP")

scamps <- c("palmitic_Highest_palmitoleic_Highest_stearic_Lowest_oleic_Lowest_linoleic_Highest_linolenic_Highest_arachidic_Highest_eicosenoic_Highest_",
            "palmitic_Highest_palmitoleic_Highest_stearic_MediumHigh_oleic_Lowest_linoleic_Lowest_linolenic_Highest_arachidic_Highest_eicosenoic_Highest_",
            "palmitic_Lowest_palmitoleic_Lowest_stearic_MediumHigh_oleic_Highest_linoleic_Lowest_linolenic_Highest_arachidic_Highest_eicosenoic_Highest_")
table(maxClustering)[scamps]
#palmitic   palmitoleic    stearic    oleic     linoleic   linolenic    arachidic   eicosenoic
#Highest    Highest        Lowest     Lowest    Highest    Highest      Highest     Highest
#Highest    Highest        MediumHigh Lowest    Lowest     Highest      Highest     Highest
#Lowest     Lowest         MediumHigh Highest   Lowest     Highest      Highest     Highest

#get summary stats for paper
crtb <- table(tsneOut$ScampLabels,tsneOut$Area)
vcp <- c("#0066ff","#ff66cc","#00aaaa")
apply(crtb[which(rownames(crtb) %in% scamps),],1,sum)
crtb[,which(colnames(crtb) =="North-Apulia")]
crtb[,which(colnames(crtb) %in% areas)]
crtb[which(rownames(crtb) %in% scamps),which(colnames(crtb) %in% areas)]
apply(crtb[which(rownames(crtb) %in% scamps),which(colnames(crtb) %in% areas)],2,sum)
apply(crtb[,which(colnames(crtb) %in% areas)],2,sum)
for (nameNum in seq(length(scamps))) {
    name <- scamps[nameNum]
    dfSub <- tsneOut[which(tsneOut$ScampLabels==name),c("tsneX","tsneY")]
    p2 <- p2 + geom_encircle(data=dfSub,aes(x=tsneX,y=tsneY),colour=vcp[nameNum],linetype=(4*nameNum-1),inherit.aes=FALSE)
}

p3 <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=Palmitic))+
    geom_point()+
    scale_colour_gradient(low = "white", high = "black")+
    theme_bw()+
    theme(legend.position="none")+
    ggtitle("Palmitic")

p4 <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=Palmitoleic))+
    geom_point()+
    scale_colour_gradient(low = "white", high = "black")+
    theme_bw()+
    ylab("")+
    theme(legend.position="none")+
    ggtitle("Palmitoleic")

p5 <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=Oleic))+
    geom_point()+
    scale_colour_gradient(low = "white", high = "black")+
    theme_bw()+
    ylab("")+
    theme(legend.position="none")+
    ggtitle("Stearic")

p6 <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=Oleic))+
    geom_point()+
    scale_colour_gradient(low = "white", high = "black")+
    theme_bw()+
    ylab("")+
    theme(legend.position="none")+
    ggtitle("Oleic")

p7 <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=Linoleic))+
    geom_point()+
    scale_colour_gradient(low = "white", high = "black")+
    theme_bw()+
    ylab("")+
    theme(legend.position="none")+
        ggtitle("Linoleic")

p8 <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=Linoleic))+
    geom_point()+
    scale_colour_gradient(low = "white", high = "black")+
    theme_bw()+
    ylab("")+
    theme(legend.position="none")+
    ggtitle("Linolenic")

p9 <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=Linoleic))+
    geom_point()+
    scale_colour_gradient(low = "white", high = "black")+
    theme_bw()+
    ylab("")+
    theme(legend.position="none")+
    ggtitle("Arachidic")

p10 <- ggplot(tsneOut,aes(x=tsneX,y=tsneY,color=Eicosenoic))+
    geom_point()+
    scale_colour_gradient(low = "white", high = "black")+
    theme_bw()+
    ylab("")+
    theme(legend.position="none")+
    ggtitle("Eicosenoic")

pjOne <- plot_grid(p,p2,nrow=1,ncol=2,labels=c("1","2"))
pjTwo <- plot_grid(p3,p4,p5,p6,p7,p8,p9,p10,nrow=2,ncol=4,labels=c("A","B","C","D","E","F","G","H"))
pj <- plot_grid(pjOne,pjTwo,nrow=2,ncol=1)
pj

      
