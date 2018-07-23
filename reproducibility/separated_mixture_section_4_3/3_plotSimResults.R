library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

#
#Scenario 1
#
fs <- readRDS("./check25_firstSimOutMat.rds")
fsp <- fs
table(fsp$alg)
fsp$ClusterTruth <- rep(10,nrow(fsp))
fsp$ClusterTruth[which(fsp$cDim > 0)] <- 17
fsp$ClusterTruth <- as.factor(fsp$ClusterTruth)
fsp$myNVI <- fsp$"VarInfo"/log(3000)
getBP <- function(rVar,ylabStr,titleStr,transformVar=FALSE) {
    if ((rVar == "Num Clusters") || (rVar == "Num Clusters > 25")) {
        mdf <- as.data.frame(fsp[,c("alg",rVar,"ClusterTruth")])
        colnames(mdf) <- c("clustering","rv","trueClusterNum")
        mdf$trueClusterNum <- as.factor(mdf$trueClusterNum)
    }
    else {
        mdf <- as.data.frame(fsp[,c("alg",rVar)])
        colnames(mdf) <- c("clustering","rv")
    }
    if (transformVar) {
        mdf$rv <- log(mdf$rv)
    }
    p <- ggplot(mdf,aes(x=as.factor(clustering),y=rv,fill=as.factor(clustering)))+geom_boxplot()+theme_bw()+
        ggtitle(titleStr)+
        xlab("")+ylab(ylabStr)+theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))
    if ((rVar == "Num Clusters") || (rVar == "Num Clusters > 25")) {
        fDF <- data.frame(trueClusterNum=c(10,17),yiVal=c(log(10),log(17)))
        fDF$trueClusterNum <- as.factor(fDF$trueClusterNum)
        p <- p + 
            facet_grid(trueClusterNum ~ .) + 
            scale_fill_manual(values=c("#440154FF","#46337EFF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))+
            geom_hline(data=fDF,aes(yintercept=yiVal),color="red",linetype="dashed")
    }
    else {
        p <- p + scale_fill_manual(values=c("#440154FF","#46337EFF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))
    }
    return(p)
}


getDistFromBest <- function(rVar,ylabStr,titleStr,transformVar=FALSE) {
    am <- fsp[,c(rVar,"alg","cDim","meanVec","aDim","tType","sRegime")]
    amt <- am %>% spread("alg",rVar)   
    amtSub <- amt[,-which(colnames(amt) %in% c("cDim","meanVec","aDim","tType","sRegime"))]
    if (rVar == "myNVI") {
        mdf <- gather(as.data.frame(t(apply(amtSub,1,function(x){x-min(x)}))),key="clustering",value=rv)
    }
    else {
        mdf <- gather(as.data.frame(t(apply(amtSub,1,function(x){max(x)-x}))),key="clustering",value=rv)
    }
    p <- ggplot(mdf,aes(x=as.factor(clustering),y=rv,fill=as.factor(clustering)))+geom_boxplot()+theme_bw()+
        ggtitle(titleStr)+
       xlab("")+ylab(ylabStr)+theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
           scale_fill_manual(values=c("#440154FF","#46337EFF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))
    return(p)
}

p0a <- getBP("Adj Rand Index","Mean ARI","Mean Adj Rand Index")
p0b <- getDistFromBest("Adj Rand Index","Difference from max mean ARI","Max ARI - method ARI")
p3 <- getBP("Elapsed Time","Log Mean Run-time (in seconds)","Log(Mean Elapsed Time)",TRUE)
p4 <- getBP("Num Clusters","Log(Mean Num Clusters)","Log(Mean Num Clusters)",TRUE)
p4a <- getBP("Num Clusters > 25","Log(Mean Num Clusters > 25)","Log(Mean Num Clusters > 25)",TRUE)
p4c <- plot_grid(p4,p4a,nrow=2,ncol=1)
p <- plot_grid(p0a,p0b,p3,p4c,nrow=1,ncol=4)
save_plot("./sim01.png",p,base_width=16,base_height=8)

#
#Scenario 2
#
fs <- readRDS("./check25_secondSimOutMat.rds")
fsp <- fs
fsp$alg[which(fsp$alg=="AP")] <- "Leveraged_AP"
fsp$ClusterTruth <- rep(10,nrow(fsp))
fsp$ClusterTruth[which(fsp$cDim > 0)] <- 17
fsp$ClusterTruth <- as.factor(fsp$ClusterTruth)
fsp$myNVI <- fsp$"VarInfo"/log(30000)
getBP <- function(rVar,ylabStr,titleStr,transformVar=FALSE) {
    if ((rVar == "Num Clusters") || (rVar == "Num Clusters > 25")) {
        mdf <- as.data.frame(fsp[,c("alg",rVar,"ClusterTruth")])
        colnames(mdf) <- c("clustering","rv","trueClusterNum")
        mdf$trueClusterNum <- as.factor(mdf$trueClusterNum)
    }
    else {
        mdf <- as.data.frame(fsp[,c("alg",rVar)])
        colnames(mdf) <- c("clustering","rv")
    }
    if (transformVar) {
        mdf$rv <- log(mdf$rv)
    }
    p <- ggplot(mdf,aes(x=as.factor(clustering),y=rv,fill=as.factor(clustering)))+geom_boxplot()+theme_bw()+
        ggtitle(titleStr)+
        xlab("")+ylab(ylabStr)+theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))
    if ((rVar == "Num Clusters") || (rVar == "Num Clusters > 25")) {
        fDF <- data.frame(trueClusterNum=c(10,17),yiVal=c(log(10),log(17)))
        fDF$trueClusterNum <- as.factor(fDF$trueClusterNum)
        p <- p + 
            facet_grid(trueClusterNum ~ .) + 
            scale_fill_manual(values=c("#FECC8FFF","#E85362FF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))+
            geom_hline(data=fDF,aes(yintercept=yiVal),color="red",linetype="dashed")
    }
    else {
        p <- p + scale_fill_manual(values=c("#FECC8FFF","#E85362FF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))
    }
    return(p)
}

getDistFromBest <- function(rVar,ylabStr,titleStr,transformVar=FALSE) {
    am <- fsp[,c(rVar,"alg","cDim","meanVec","aDim","tType","sRegime")]
    amt <- am %>% spread("alg",rVar)   
    amtSub <- amt[,-which(colnames(amt) %in% c("cDim","meanVec","aDim","tType","sRegime"))]
    if (rVar == "myNVI") {
        mdf <- gather(as.data.frame(t(apply(amtSub,1,function(x){x-min(x)}))),key="clustering",value=rv)
    }
    else {
        mdf <- gather(as.data.frame(t(apply(amtSub,1,function(x){max(x)-x}))),key="clustering",value=rv)
    }
    p <- ggplot(mdf,aes(x=as.factor(clustering),y=rv,fill=as.factor(clustering)))+geom_boxplot()+theme_bw()+
        ggtitle(titleStr)+
        xlab("")+ylab(ylabStr)+theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
        scale_fill_manual(values=c("#FECC8FFF","#E85362FF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))
    return(p)
}

p0a <- getBP("Adj Rand Index","Mean ARI","Mean Adj Rand Index")
p0b <- getDistFromBest("Adj Rand Index","Difference from max mean ARI","Max ARI - method ARI")
p3 <- getBP("Elapsed Time","Log Mean Run-time (in seconds)","Log(Mean Elapsed Time)",TRUE)
p4 <- getBP("Num Clusters","Log(Mean Num Clusters)","Log(Mean Num Clusters)",TRUE)
p4a <- getBP("Num Clusters > 25","Log(Mean Num Clusters > 25)","Log(Mean Num Clusters > 25)",TRUE)
p4c <- plot_grid(p4,p4a,nrow=2,ncol=1)
p <- plot_grid(p0a,p0b,p3,p4c,nrow=1,ncol=4)
save_plot("./sim02.png",p,base_width=16,base_height=8)

#
#Scenario 3
#
fs <- readRDS("./check25_thirdSimOutMat.rds")
fsp <- fs
fsp$alg[which(fsp$alg=="AP")] <- "Leveraged_AP"
fsp$ClusterTruth <- rep(30,nrow(fsp))
fsp$ClusterTruth[which(fsp$cDim > 0)] <- 56
fsp$ClusterTruth <- as.factor(fsp$ClusterTruth)
fsp$myNVI <- fsp$"VarInfo"/log(30000)
getBP <- function(rVar,ylabStr,titleStr,transformVar=FALSE) {
    if ((rVar == "Num Clusters") || (rVar == "Num Clusters > 25")) {
        mdf <- as.data.frame(fsp[,c("alg",rVar,"ClusterTruth")])
        colnames(mdf) <- c("clustering","rv","trueClusterNum")
        mdf$trueClusterNum <- as.factor(mdf$trueClusterNum)
    }
    else {
        mdf <- as.data.frame(fsp[,c("alg",rVar)])
        colnames(mdf) <- c("clustering","rv")
    }
    if (transformVar) {
        mdf$rv <- log(mdf$rv)
    }
    p <- ggplot(mdf,aes(x=as.factor(clustering),y=rv,fill=as.factor(clustering)))+geom_boxplot()+theme_bw()+
        ggtitle(titleStr)+
        xlab("")+ylab(ylabStr)+theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))
    if ((rVar == "Num Clusters") || (rVar == "Num Clusters > 25")) {
        fDF <- data.frame(trueClusterNum=c(30,56),yiVal=c(log(30),log(56)))
        fDF$trueClusterNum <- as.factor(fDF$trueClusterNum)
        p <- p + 
            facet_grid(trueClusterNum ~ .) + 
            scale_fill_manual(values=c("#FECC8FFF","#E85362FF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))+
            geom_hline(data=fDF,aes(yintercept=yiVal),color="red",linetype="dashed")
    }
    else {
        p <- p + scale_fill_manual(values=c("#FECC8FFF","#E85362FF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))
    }
    return(p)
}

getDistFromBest <- function(rVar,ylabStr,titleStr,transformVar=FALSE) {
    am <- fsp[,c(rVar,"alg","cDim","meanVec","aDim","tType","sRegime")]
    amt <- am %>% spread("alg",rVar)   
    amtSub <- amt[,-which(colnames(amt) %in% c("cDim","meanVec","aDim","tType","sRegime"))]
    if (rVar == "myNVI") {
        mdf <- gather(as.data.frame(t(apply(amtSub,1,function(x){x-min(x)}))),key="clustering",value=rv)
    }
    else {
        mdf <- gather(as.data.frame(t(apply(amtSub,1,function(x){max(x)-x}))),key="clustering",value=rv)
    }
    p <- ggplot(mdf,aes(x=as.factor(clustering),y=rv,fill=as.factor(clustering)))+geom_boxplot()+theme_bw()+
        ggtitle(titleStr)+
        xlab("")+ylab(ylabStr)+theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
        scale_fill_manual(values=c("#FECC8FFF","#E85362FF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))
    return(p)
}

p0a <- getBP("Adj Rand Index","Mean ARI","Mean Adj Rand Index")
p0b <- getDistFromBest("Adj Rand Index","Difference from max mean ARI","Max ARI - method ARI")
p3 <- getBP("Elapsed Time","Log Mean Run-time (in seconds)","Log(Mean Elapsed Time)",TRUE)
p4 <- getBP("Num Clusters","Log(Mean Num Clusters)","Log(Mean Num Clusters)",TRUE)
p4a <- getBP("Num Clusters > 25","Log(Mean Num Clusters > 25)","Log(Mean Num Clusters > 25)",TRUE)
p4c <- plot_grid(p4,p4a,nrow=2,ncol=1)
p <- plot_grid(p0a,p0b,p3,p4c,nrow=1,ncol=4)
save_plot("./sim03.png",p,base_width=16,base_height=8)

