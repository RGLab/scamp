#' Selective Clustering Annotated using Modes of Projections
#'
#' This function clusters a dataSet using scamp, selective clustering
#' annotated using modes of projections. 
#' 
#' @param dataSet The input dataSet to be clustered by scamp.
#'   Must be an R matrix. The dataSet should only contain numerical data.
#' 
#' @param numberIterations The number of scamp iterations to perform.
#' 
#' @param clusterOutputString scamp generates a number of intermediate
#'   clusterings of the dataSet. When applied to large data sets,
#'   often it is useful to examine these intermediate clusterings.
#'   The clusterOutputString is a path to where this output will be stored.
#'
#' @param pValueThreshold The significance level for the DIP test of Hartigan
#'   and Hartigan (1985). scamp uses this level to choose to partition a projection
#'   in the dataSet with the taut-string of Davies and Kovac (2004). Setting it
#'   to smaller values (such as a conventional 0.05) usually produces faster
#'   exaustive clusterings of the dataSet.
#'
#' @param minimumClusterSize The size of the smallest admissible candidate cluster.
#'   The recursive search for candidate clusters will not consider subsets of observations
#'   below this value as admissible candidates -- it terminates along such branches without
#'   recording the subsets encountered. Increasing this parameter usually increases
#'   the speed of a single scamp iteration.
#'
#' @param maximumSearchDepth An upper bound on the depth of a search tree for candidate clusters.
#'   The recursive search terminates if it exceeds this depth.
#'
#' @param maximumClusterNum An upper bound on the number of candidate clusters scamp will find.
#'   If a scamp iteration is exhaustively searching the space of clustering trees, the search will
#'   terminate if scamp finds a collection of candidate clusters that exceed this bound.
#'
#'   If a scamp iteration performs a randomCandidateSearch (or if any values are restricted in the search
#'   process using the resValMatrix) \emph{\bold{the scamp search will continue until it
#'   finds this number of candidate clusters}}.
#'
#' @param maximumAntimodeNumber An upper bound on the number of possible antimodes found by the
#'   taut-string estimate of a projection. If the taut-string produces more than the
#'   maximumAntimodeNumber number of antimodes, the recursive search terminates along that branch.
#'
#' @param maximumSearchTime The maximum time in seconds a scamp iteration can search for candidate
#'   clusters. If a single iteration exceeds this value in either the intial search phase or
#'   and residual search phases, that scamp iteration is abandonded and the next iteration is 
#'   attempted.
#'
#'   On large dataSets (with either a large number of observations or variables), the
#'   initial addition of noise can produce configurations of the dataSet which take much longer
#'   to cluster than the average scamp iteration. This parameter is useful to prevent such iterations
#'   from extending a single scamp call from running too long.
#'
#' @param annotationVec A vector of strings that will be used to annotate the projections (columns) in the dataSet.
#'   For the scamp clustering to be interpetable, meaningful strings must be provided here.
#'
#' @param finalAnnQs A vector of three increasing values c(qLow,qMed,qHigh), with  0 < qLow < qMed < qHigh < 1.
#'   When scamp annotates selected clusters, it computes the sample quantiles qLow, qMed, and qHigh of each
#'   coordinate project for each selected cluster.
#'
#'   If scamp has found a single annotation boundary for a coordinate projection in the dataset,
#'   both qLow and qHigh are compared to that boundary. If qLow exceeds the annotation boundary, the cluster is labeled
#'   "Highest" along that coordinate. If qHigh falls below the annotation boundary, the cluster is labeled "Lowest"
#'   along that coordinate. Otherwise, the observations in the selected cluster are split into two sub-clusters:
#'   observations below the annotation boundary comprise one of the new sub-clusters, observation above the other.
#'
#'   If scamp has found multiple annotation boundaries for a coordinate, qMed is compared to them and used to annotate
#'   a selected cluster along that coordinate.
#' 
#' @param numberOfThreads The addition of noise to the dataSet, the random search for candidate clusters,
#'   and the annotation of selected clusters have been parallelized. For numberOfThreads greater than 2,
#'   numberOfThreads-1 worker threads will add noise to the dataSet, conduct the search for candidate clusters,
#'   and annotate the selected clusters. The value of 0 indicates scamp should detect the maximum number
#'   of threads supported by the hardware and then use that value.
#' 
#' @param getVotingHistory Boolean flag. If true, the complete voting history for 
#'   \emph{every} scamp iteration is written the clusterOutputString directory. Only the most recent
#'   voting history and the one preceding it are stored.
#'
#' @param getDebugInfo Boolean flag. If true, scamp reports job progress to the R interpreter/log file.
#' 
#' @param randomCandidateSearch Boolean flag. If true, the initial search for candidate clusters will
#'   sample from the space of clustering trees. At each node in the clustering tree, scamp will pick 
#'   the next clustering branch uniformly at random the set of admissible projections (according to the dip).
#'   Especially in small dataSets, this can cause the scamp search to find the same clustering tree multiple
#'   times.
#'
#' @param randomResidualCandidateSearch Boolean flag. If true, the search for candidate clusters of observations
#'   not clustered by the initial search and selection will also randomly sample from the space of clustering
#'   trees. Users can choose the combination of exhaustive search and random search for the two stages of 
#'   clustering that best suits their dataSet.
#'
#' @param anyValueRestricted Boolean flag. If true, indicates the user has some prior knowledge about bounds of
#'  admissible values. Such values must be identified by coordinate projection in the resValMatrix.
#'
#'  The motivation for this parameter is that some technologies produce dataSets which are: very zero-inflated;
#'  have a minimum limit of detection; have an upper limit of detection. Hence, one might wish to indicate that
#'  all values below some threshold be treated as low \emph{by default}, or above a threshold be treated as high
#'  \emph{by default}. Of critical importance: restricted values \emph{\bold{are not used by the taut-string density estimator
#'  in the search for candidate clusters}}.
#'  
#' @param resValMatrix R matrix indicating restricted values. Contains only three entries: {0,1,2}. If entry (i,j) is 0, it indicates
#'  that observation i is unrestricted along coordinate j (and so is used by the taut-string in the search for cnadidate clusters).
#'
#'  If entry (i,j) is 1, it indicates observation i is \emph{restricted from below} along coordinate j: this observation will 
#'  automatically be treated as "Lowest" along coordinate j, and will not be used by the taut-string in the search for candidate clusters.
#'
#'  If entry (i,j) is 2, it indicates observation i is \emph{restricted from above} along coordinate j: this observation will 
#'  automatically be treated as "Highest" along coordinate j, and will not be used by the taut-string in the search for candidate clusters.
#'
#' @param setAnnotationBoudaries Boolean flag. If true, indicates the user has supplied a list of numeric vectors which  
#'  will serve as the annotation boundaries \emph{for each} column in the data set.
#' 
#' @param annotationBoundaryValues R list of numeric vectors. Each entry in the list corresponds to the columns in the data set. Each vector
#'  contains at least one value that will serve as the annotation boundary for that column.
#'
#' @param gaussianScaleParameter Noise is added in a two-step procedure for SCAMP. First, uniform noise is added to
#' break ties observed in the columns of the dataSet. After uniform noise is added, Gaussian noise is added to perturb
#' the relative order of the observations in each column of the dataSet. The standard-devaition of the Gaussian noise
#' depends on the distance to the neighboring order statistics of an observation.
#'
#' The default value of 4.0 preserves the relative ordering of observations in a column, with high probability. As the
#' value decreases, observations are more and more likely to switch position in the order-statistics. 
#'
#'
#' @param randomSeed Set to a large integer value to reproduce scamp runs when numberOfThreads==1.
#' Default value of 0 leads to non-reproducible SCAMP runs regardless of numberOfThreads value.
#'
#' @param allowRepeatedSplitting Set to FALSE by default. If TRUE, indicates the search for candidate clusters can
#' induce modal clusters along the same column vector of dataSet multiple times throughout the search. No matter
#' the value of this parameter (TRUE OR FALSE), the search for candidate clusters is never allowed to induce modal clusters
#' along the same column vector twice in a row: from the point of view of SCAMP, this can only arise due to technical
#' error.
#'
#' @return scamp returns an R list with two entries. Each entry is a vector of strings.
#'  The first entry, named RunOffVote, is a clustering of the dataSet according to a run off vote across the scamp iterations.
#'  The second entry, named MaxVote, is a clustering of the dataSet according to the maximum vote across the scamp iterations.
#'
#' @examples
#' clusterMatrix <- as.matrix(iris[,-5])
#' scampClustering <- scamp(dataSet=clusterMatrix,
#'                          numberIterations=100,
#'                          clusterOutputString="./scampTest",
#'                          numberOfThreads=0) #0 automatically sets to
#'                                                          #hardware max.
#' table(scampClustering$RunOffVote)
#' table(scampClustering$MaxVote)
#' @export
#' @md

scamp <- function(dataSet,
                  numberIterations = 100,
                  clusterOutputString=c("."),
                  pValueThreshold = 0.25,
                  minimumClusterSize = 25,
                  maximumSearchDepth = 1e6,
                  maximumClusterNum = (50 * ncol(dataSet)),
                  maximumAntimodeNumber = 100,
                  maximumSearchTime=1e6,
                  annotationVec = colnames(dataSet),
                  finalAnnQs=c(0.499,0.5,0.501),
                  numberOfThreads=1,
                  getVotingHistory=FALSE,
                  getDebugInfo=FALSE,
                  randomCandidateSearch=TRUE,
                  randomResidualCandidateSearch=TRUE,
                  anyValueRestricted=FALSE,
                  resValMatrix = matrix(0,nrow=2,ncol=2),
                  setAnnotationBoudaries=FALSE,
                  annotationBoundaryValues=list(rep(0,ncol(dataSet))),
                  gaussianScaleParameter=4,
                  randomSeed=0,
                  allowRepeatedSplitting=FALSE){
    if (!is.matrix(dataSet))
        stop("Must provide a numeric R matrix")
    if (!is.double(pValueThreshold))
        stop("The p-value threshold must be a double precison number.")
    if ((pValueThreshold < 0.01) || (pValueThreshold > 0.99))
        stop("The p-value threshold for the dip test must be between  0.01 and 0.99")
    if ((!is.numeric(minimumClusterSize)) || (minimumClusterSize < 10))
        stop("A cluster cannot be smaller than 10 observations")
    if ((!is.numeric(maximumSearchDepth)) || (maximumSearchDepth < 1))
        stop("Search depth must be a natural number > 1")
    if ((!is.numeric(maximumClusterNum)) || (maximumClusterNum < 2))
        stop("Must search for more than 1 candidate")
    if (nrow(dataSet) < minimumClusterSize)
        stop("Too few observations for data")
    if (max(is.na(annotationVec)))
        stop('User must assign column names to input data set. colnames(dataSet) <- c("user", "entry",...)')
    if (length(finalAnnQs) != 3)
        stop("Must provide three quantile values for annotation boundaries.")
    if (length(unique(finalAnnQs)) != 3)
        stop("Must provide three distinct quantile values for  annotation boundaries.")
    if (min((sort(finalAnnQs) == finalAnnQs)) == 0)
        stop("finalAnnQs must be in increasing order.")
    if (gaussianScaleParameter < 0)
        stop("Gaussian scale parameter smaller than 0. Not possible.")
    if (!(identical(randomCandidateSearch,TRUE) || identical(randomCandidateSearch,FALSE))) 
        stop("User must set randomCandidateSearch either to TRUE or FALSE.")
    if (!(identical(randomResidualCandidateSearch,TRUE) || identical(randomResidualCandidateSearch,FALSE))) 
        stop("User must set randomResidualCandidateSearch either to TRUE or FALSE.")
    if (!(identical(anyValueRestricted,TRUE) || identical(anyValueRestricted,FALSE))) 
        stop("User must set anyValueRestricted either to TRUE or FALSE.")
    if (!(identical(setAnnotationBoudaries,TRUE) || identical(setAnnotationBoudaries,FALSE))) 
        stop("User must set setAnnotationBoudaries either to TRUE or FALSE.")
    if (!(identical(allowRepeatedSplitting,TRUE) || identical(allowRepeatedSplitting,FALSE))) 
        stop("User must set allowRepeatedSplitting either to TRUE or FALSE.")
    if (anyValueRestricted) {
        if ((ncol(resValMatrix) != ncol(dataSet)) || (nrow(resValMatrix) != nrow(dataSet))) {
            stop("The restricted value matrix must be the same dimension as the input dataSet.")
        }
        if (length(setdiff(unique(unlist(apply(resValMatrix,2,unique))),c(0,1,2)))) {
            stop("The restricted value matrix must only contain integer values in the set {0,1,2}")
        }
        resValMatOrderViolations <- function(x,dataV){
            zeroLookup <- which(x==0)
            oneLookup <- which(x==1)
            twoLookup <- which(x==2)
            zeroMin <- min(dataV[zeroLookup])
            zeroMax <- max(dataV[zeroLookup])
            if (length(oneLookup) && length(twoLookup))  {
                oneMax <- max(dataV[oneLookup])
                twoMin <- min(dataV[twoLookup])
                if ((oneMax > zeroMin) || (oneMax > twoMin)) return(TRUE)
                else if (twoMin < zeroMax) return(TRUE)
                else return(FALSE)
            }
            else if (length(oneLookup)) {
                oneMax <- max(dataV[oneLookup])
                if (oneMax > zeroMin) return(TRUE)
                else return(FALSE)
            }
            else if (length(twoLookup)) {
                twoMin <- min(dataV[twoLookup])
                if (twoMin < zeroMax) return(TRUE)
                else return(FALSE)
            }
            else {
                return(FALSE)
            }
        }
        scampOrderViolations <- rep(0,ncol(dataSet))
        for (i in seq(ncol(dataSet))) {
            scampOrderViolations[i] <- resValMatOrderViolations(resValMatrix[,i],dataSet[,i])
        }
        if (max(scampOrderViolations)) {
            print("In each column of the resValMatrix, all values marked 1 must correspond to values in the dataSet column")
            print("that fall below values in the same dataSet column with 0 or 2 in the resValMatrix.")
            print("Similarly, in each column of the resValMatrix, all values marked 2 must correspond to values in the dataSet column")
            print("that fall above values in the same dataSet column with 0 or 2 in the resValMatrix.")
            stop("Must initialize the restricted value matrix according to these constraints for a meaningful clustering of the dataSet.")
        }
    }
    if (setAnnotationBoudaries) {
        if (length(annotationBoundaryValues) != ncol(dataSet)) {
            stop("Must provide annotation boundaries values for every column in the dataSet.")
        }
        if (!min(unlist(lapply(annotationBoundaryValues,function(x){min(sapply(x,is.numeric))})))) {
            stop("Boundary values must be numeric values.")
        }
        if (!min(unlist(lapply(annotationBoundaryValues,function(x){length(x)==length(unique(x))})))) {
            stop("Boundary values must be unique numeric values.")
        }
        if (!min(unlist(lapply(annotationBoundaryValues,function(x){min(sort(x)==x)})))) {
            stop("boundary values must be in increasing order.")
        }
    }
    if (!dir.exists(file.path(clusterOutputString))) {
        dir.create(file.path(clusterOutputString))
    }
    else {
        print(paste0("Scamp results will be stored in ",file.path(clusterOutputString)))
    }
    if (!dir.exists(file.path(clusterOutputString,"scampResults"))) {
        dir.create(file.path(clusterOutputString,"scampResults"))
    }
    else {
        print(paste0("Writing all scamp output to ",file.path(clusterOutputString,"scampResults")))
    }

    aCols <- apply(dataSet,2,length)
    uCols <- apply(dataSet,2,function(x){length(unique(x))})
    rCols <- uCols/aCols
    mCol <- min(rCols)
    if (mCol < 0.05){
        print("Warning: there are a large number of duplicated values in the input data matrix.\n
               This can lead to large variation from fuzzing -- results will vary from run to run.\n
               Setting the 'numberIterations' parameter to a large value may be warranted.") 
    }

    #scale the dataSet to ensure L-moment scores are comparable across columns.
    scaledDataSet <- scale(dataSet)

    if (setAnnotationBoudaries) {
        print("Annotation boundary values before scaling.")
        scaledColNames <- colnames(scaledDataSet)
        for (j in seq(length(annotationBoundaryValues))) {
            print(paste0("Column ",scaledColNames[j]," values: "))
            print(annotationBoundaryValues[[j]])
        }
        for (i in seq(ncol(dataSet))) {
            cData <- dataSet[,i]
            vals <- annotationBoundaryValues[[i]]
            newBoundaryVals <- c()
            for (k in seq(length(vals))) {
                val <- vals[k]
                probVal <- length(which(cData <= val))/length(cData)
                newBoundaryVals <- append(newBoundaryVals,as.numeric(quantile(scaledDataSet[,i],probs=probVal)))
            }
            annotationBoundaryValues[[i]] <- newBoundaryVals
        }
        print("Annotation boundary values after scaling.")
        for (j in seq(length(annotationBoundaryValues))) {
            print(paste0("Column ",scaledColNames[j]," values: "))
            print(annotationBoundaryValues[[j]])
        }
    }

    
    startTime <- proc.time()
    timeDiff <- proc.time()-startTime
    elapsedTime <- as.numeric(timeDiff[3])
    scampLabels <- cppNoisyScamp(scaledDataSet,
                                 as.double(pValueThreshold),
                                 as.integer(minimumClusterSize),
                                 allowRepeatedSplitting,
                                 as.integer(maximumSearchDepth),
                                 as.integer(maximumClusterNum),
                                 annotationVec,
                                 maximumAntimodeNumber, 
                                 randomCandidateSearch,
                                 randomResidualCandidateSearch,
                                 finalAnnQs,
                                 numberOfThreads,
                                 anyValueRestricted,
                                 resValMatrix,
                                 setAnnotationBoudaries,
                                 annotationBoundaryValues,
                                 numberIterations,
                                 clusterOutputString,
                                 getVotingHistory,
                                 maximumSearchTime,
                                 getDebugInfo,
                                 gaussianScaleParameter,
                                 randomSeed)
    return(scampLabels)
}



