# preparing data
library(dplyr)
library(tidyr)
library(ggplot2)
library(Rtsne)
print("reading samples")
samples <- read.delim('../GTEx_Data_V6_Annotations_SampleAttributesDS.txt', sep='\t') %>%
    select(SAMPID, primary.tissue=SMTS, tissue=SMTSD)
print("reading gtex")
gtex <- read.table('../GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct',
                   skip=2,
                   colClasses=c('character', 'character', rep('numeric', 2921)),
                   stringsAsFactors=F,
                   header=T)  
print("renaming")
newColNames <- as.character(sapply(colnames(gtex),function(x){gsub('\\.','-',x)}))
colnames(gtex) <- newColNames
print("saving as rds")
saveRDS(gtex,"./geneData.rds")
saveRDS(samples,"./samplesData.rds")
geneData <- readRDS("./geneData.rds")
dim(geneData)
print("making count data")
countData <- geneData[,-c(1,2)]
countData <- matrix(as.numeric(unlist(countData)),nrow=nrow(countData))
print("saving count data")
saveRDS(countData,"./allCounts.rds")
countData <- readRDS("./allCounts.rds")
dim(countData)
zeroCounts <- apply(countData,1,function(x){length(which(x==0))})
print("removing 8555 zero counts")
subCount <- countData[-which(zeroCounts>8554),]
dim(subCount)
l2SubCount <- apply(subCount,2,function(x){log2(1+x)})
l2NormCount <- apply(l2SubCount,2,function(x){x/sum(x)})
print("transposing")
transposeData <- t(l2NormCount)
rm(subCount)
rm(l2SubCount)
rm(l2NormCount)
rm(countData)
rm(gtex)
rm(samples)
rm(geneData)
gc()
set.seed(78123)
print("making tsne")
tsneData <- Rtsne(transposeData)
saveRDS(tsneData,"./tsneData.rds")

