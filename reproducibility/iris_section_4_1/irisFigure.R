library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ClusterR)
library(scamp)
library(viridis)
library(tidyr)
clusterMatrix <- as.matrix(iris[,-5])
#
#Replicate Figure 4, Page 24
#
firstOb <- TRUE
for (i in seq(1000)) {
    pValues <- apply(clusterMatrix,2,function(x){singleDip(sort(addNoiseToDataVector(x,4,i)))})
    if (firstOb) {
        pValMat <- t(as.matrix(pValues))
        firstOb <- FALSE
    }
    else {
        pValMat <- rbind(pValMat,pValues)
    }
}
pValues <- apply(pValMat,2,mean)
plotDF <- as.data.frame(pValMat[,names(sort(pValues))])
colnames(plotDF) <- paste0("0",seq(ncol(plotDF)),"_",colnames(plotDF))
mdf <- gather(plotDF)
p <- ggplot(mdf,aes(x=key,y=value)) + geom_boxplot() + theme_bw() + ylab("Dip test p-value") +
    xlab("") + ggtitle("Iris data set, dip test p-value distribution across 1000 noise applications")+
    geom_hline(yintercept=0.01,color="blue",linetype="dashed")+
    geom_hline(yintercept=0.05,color="green",linetype="dotted")+
        geom_hline(yintercept=0.25,color="red")
p

#
#Replicate Figure 5, Page 25
#

#note: if number of threads changed, results may not match paper exactly.
#however, the results should be very similar.
#this is because the random candidate cluster search allows all launched
#threads to finish after the maximumClusterNum threshold has been crossed.
#each launched thread increments the random seed, and so
#modifying the threads parameter will affect all subsequent scamp iterations.
scampClustering <- scamp(dataSet=clusterMatrix,
                         numberIterations=500,
                         numberOfThreads=4,
                         getDebugInfo=FALSE,
                         randomSeed=8123)
scampMax <- scampClustering[[2]]
scampMax <- as.character(sapply(scampMax,function(x){gsub("_$","",x)}))
my.iris <- iris[-5]
my.iris$scampLabel <- as.factor(scampMax)
my.iris$Iris <- as.factor(iris[,5])
myColScale <- scale_color_manual(name="Iris",values= brewer.pal(n=3,name="Set2"))
myColScale2 <- scale_color_manual(name="scampLabel",values= brewer.pal(n=3,name="Set2"))
external_validation(as.numeric(as.factor(iris[,5])),as.numeric(as.factor(scampMax)),summary_stats=TRUE)
scampLabelPlot <- my.iris$scampLabel
scampLabelPlot <- sapply(scampLabelPlot,function(x){gsub("Petal.Length_Highest_Petal.Width_Highest","Petal.Length Highest\nPetal.Width Highest",x)})
scampLabelPlot <- sapply(scampLabelPlot,function(x){gsub("Petal.Length_Highest_Petal.Width_Medium","Petal.Length Highest\nPetal.Width Medium",x)})
scampLabelPlot <- sapply(scampLabelPlot,function(x){gsub("Petal.Length_Lowest_Petal.Width_Lowest","Petal.Length Lowest\nPetal.Width Lowest",x)})
my.iris$scampLabel <- scampLabelPlot
t1 <- ggplot(my.iris,aes(x=Petal.Length,y=Petal.Width,color=Iris))+
    geom_point()+
    theme_bw()+
    myColScale
p1 <- ggplot(my.iris,aes(x=Petal.Length,y=Petal.Width,color=scampLabel))+
    geom_point()+theme_bw()+
    scale_color_manual(values=c("#440154FF","#21908CFF","#FDE725FF"))
t2 <- ggplot(my.iris,aes(x=Sepal.Length,y=Sepal.Width,color=Iris))+
    geom_point()+
    theme_bw()+
    myColScale
p2 <- ggplot(my.iris,aes(x=Sepal.Length,y=Sepal.Width,color=scampLabel))+
    geom_point()+theme_bw()+
    scale_color_manual(values=c("#440154FF","#21908CFF","#FDE725FF"))
irisGridA <- plot_grid((t1+theme(legend.position="none")),
                      (t2+theme(legend.position="none")),nrow=2,ncol=1,labels=c("1","3"))
legendGA <- get_legend(t1 + theme(legend.position="bottom"))
pA <- plot_grid(irisGridA, legendGA, ncol = 1, rel_heights = c(1, .05))
irisGridB <- plot_grid((p1+theme(legend.position="none")),
                      (p2+theme(legend.position="none")),nrow=2,ncol=1,labels=c("2","4"))
legendGB <- get_legend(p1 + theme(legend.position="bottom"))
pB <- plot_grid(irisGridB, legendGB, ncol = 1, rel_heights = c(1, .05))
irisGrid <- plot_grid(pA,pB,nrow=1)
irisGrid

