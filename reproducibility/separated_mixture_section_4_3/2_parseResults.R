library(stringr)
library(dplyr)
library(tidyr)

jn1 <- "./sepMix_bigSample_0_bigHet_0_twoMV_0"
jn2 <- "./sepMix_bigSample_0_bigHet_0_twoMV_1"
jn3 <- "./sepMix_bigSample_1_bigHet_0_twoMV_0"
jn4 <- "./sepMix_bigSample_1_bigHet_0_twoMV_1"
jn5 <- "./sepMix_bigSample_1_bigHet_1_twoMV_0"
jn6 <- "./sepMix_bigSample_1_bigHet_1_twoMV_1"

getJob <- function(jobName) {
    simulations <- list.files(jobName)
    firstSim <- TRUE
    for (sim in simulations) {
        jobc <- str_split_fixed(sim,"_",8)
        params <- jobc[seq(1,7,by=2)]
        settings <- jobc[seq(2,8,by=2)]
        jobDF <- data.frame(tmpcol=0)
        for (i in seq(length(params))) {
            jobDF[,params[i]] <- settings[i]
        }
        jobDF <- jobDF[,-1]
        res <- as.data.frame(t(readRDS(paste0(jobName,"/",sim,"/jobOutput.rds"))))
        res$alg <- rownames(res)
        resMat <- cbind(res,jobDF)
        if (firstSim) {
            outMat <- resMat
            firstSim <- FALSE
        }
        else  {
            outMat <- rbind(outMat,resMat)
        }
    }
    return(outMat)
}
om1 <- getJob(jn1)
om2 <- getJob(jn2)
om3 <- getJob(jn3)
om4 <- getJob(jn4)
om5 <- getJob(jn5)
om6 <- getJob(jn6)

om1$meanVec <- "One"
om2$meanVec <- "Two"
om3$meanVec <- "One"
om4$meanVec <- "Two"
om5$meanVec <- "One"
om6$meanVec <- "Two"
outMat <- rbind(om1,om2)
saveRDS(outMat,"./check25_firstSimOutMat.rds")
outMat2 <- rbind(om3,om4)
saveRDS(outMat2,"./check25_secondSimOutMat.rds")
outMat3 <- rbind(om5,om6)
saveRDS(outMat3,"./check25_thirdSimOutMat.rds")

