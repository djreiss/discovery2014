library(entropy)
library(randomForest)
library(doMC)


#' Generate estimate FWER for ranked features using the mProbes method outlined in 
#'  Van Anh Huynh-Thu, Yvan Saeys, Louis Wehenkel, and Pierre Geurts
#'  Statistical interpretation of machine learning-based feature importance scores for biomarker discovery
#'  Bioinformatics Advance Access published April 25, 2012
#'  See notebook #3, page 164 for algorithm
#'  NOTE:  mProbes is pretty strict, so it wouldn't be aweful to select ANY that have
#'    a score less than 1.
#'  Can I use CMIM?
#'
#' @param dataMat  The attribute matrix.  Columns are features, rows are examples
#' @param classIdx  The index for the class variable (DEFAULT: ncol(dataMat))
#' @param type  The feature selection:  Currently 'MI', 'randomForest', or 'Dynamic'.  Dynamic sets MI if there's more than 700 features, or randomForest otherwise. (DEFAULT: "Dynamic")
#' @param dynamicCutoff  The cutoff if 'type' is set to "Dynamic". (DEFAULT: 800)
#' @param ntree  The number of trees to include in a random forest if used. (DEFAULT: 1000)
#' @return a ranked list of features with their estimated FWER
#' @export
#' @usage features.FWER <- mProbes(dataMat, classIdx=ncol(dataMat), type="Dynamic", dynamicCutoff=800, ntree=1000)
mProbes <- function (dataMat, classData, type="Dynamic", dynamicCutoff=800, ntree=1000) {
    gc(reset=TRUE)

    classDF <- classData ##dataMat[, classIdx, drop=FALSE]
    ##className <- 'class' ##colnames(dataMat)[classIdx]
    dataPart <- dataMat ##[,1:ncol(dataMat) != classIdx]
    features.FWER <- rep(1, ncol(dataPart))
    names(features.FWER) <- colnames(dataPart)

    ##Dynamic select MI method.  Used MI only if there's too many features from randomForest
    if (type == "Dynamic") {
        if (ncol(dataPart) > dynamicCutoff) {
            type <- 'MI'
        } else {
            type <- 'randomForest'
        }
    }

    if (type == 'randomForest') {
        numFeat <- ncol(dataPart)
        if (ntree <= 2 * numFeat) {
            cat(numFeat, 'features but only', ntree, 'trees.  Consider using more trees.', '\n')
        }
    }

    cat('\nRunning mProbes with feature selector: ', type, ' on ', ncol(dataPart), ' candidates\n', sep="")

    ## I. Use randomForest to rank features 
    ## For Xi in all X features
    ##   1) Generate |X| permutations of Xi
    ##   2) Rank the 2*|X| features
    rankedFeat <- foreach (i = 1:ncol(dataPart)) %dopar% {
        if(i %% ceiling(ncol(dataPart)/100) == 0) { cat(round(i/ncol(dataPart),2)*100,'%,',sep="") }

        Xi <- dataPart[,i]
        randData <- as.matrix(sapply(1:ncol(dataPart), function(x) {Xi[order(runif(length(Xi)))]}))
        rownames(randData) <- rownames(dataPart)
        colnames(randData) <- paste('random.', 1:ncol(randData), sep="")
        curExp <- cbind(dataPart, randData)
        curExp <- apply(curExp, 2, as.numeric)
        rownames(curExp) <- rownames(dataPart)
        curExp <- cbind(as.data.frame(curExp), classDF)
        curForm <- as.formula(paste(colnames(curExp)[ncol(curExp)], "~", "."))
        if (type == "MI") {
            featureList <- xvalMIranker(curExp, classIdx=nrow(curExp), folds=2, verbose=FALSE)
            features.ranked <- featureList[,1]
            names(features.ranked) <- rownames(featureList)
            features.ranked <- sort(features.ranked, decreasing=TRUE)
        } else {
                                        #rf <- randomForest(curForm, data=curExp)
            rf <- randomForest(x=curExp[,1:(ncol(curExp)-1),drop=FALSE], y=curExp[,ncol(curExp)], ntree=ntree, importance=TRUE)
                                        #rf <- randomForest(curForm, data=curExp, ntree=20)
            features.ranked <-  importance(rf)[order(importance(rf)[,'MeanDecreaseAccuracy']),'MeanDecreaseAccuracy']
        }

        ##Enumerate the features that are better than random
        bestRandom <- grep('^random\\.[0-9]*$', names(features.ranked))
        bestRandom <- bestRandom[length(bestRandom)]
        if (bestRandom < length(features.ranked)) {
            features <- names(features.ranked)[(bestRandom+1):length(features.ranked)]
        } else {
            features <- NULL
        }

        return(features)
    } #rankedFeat <- foreach (i = 1:ncol(dataPart)) %dopar% {
    cat('\n')

                                        # II. Calculate the FWERs using the previous list
                                        # For Xi in all X features
                                        #   1) FWER = fraction of the |X| runs where a rand feature ranks above Xi
    cur.FWERs <- 1 - (table(unlist(rankedFeat)) / length(rankedFeat))
    features.FWER[match(names(cur.FWERs), names(features.FWER))] <- cur.FWERs

    return(features.FWER)
}

#' Rank features using mutual information
#'   Cross validate to pick consistently best features
#'
#' @param dataMat  The attribute matrix
#' @param classIdx  The index for the class variable (DEAFULT: nrow(dataMat))
#' @param folds  The number of folds to break up the data (DEFAULT: 10)
#' @param verbose  Set to TRUE to show output (DEFAULT: FALSE)
#' @export
#' @usage featureList <- xvalMIranker(dataMat, classIdx=nrow(dataMat), folds=10, verbose=FALSE)
xvalMIranker <- function(dataMat, classIdx=ncol(dataMat), folds=10, verbose=FALSE) {
    cuts <- cut(1:nrow(dataMat), folds)
    randIdxs <- order(runif(nrow(dataMat)))
    
    featureMatrix <- NULL
    for (curCut in unique(cuts)) {
        if (verbose == TRUE) { cat(curCut, ", ") }
        idxs <- randIdxs[cuts!=curCut]
        curTest <- dataMat[idxs,]
        ##featureList <- order(useMutualInformation(curTest), decreasing=FALSE)
        featureList <- useMutualInformation(curTest)
        names(featureList) <- colnames(curTest)
        featureMatrix <- cbind(featureMatrix, featureList)
                                        #rownames(featureMatrix) <- colnames(curTest)
    }
    colnames(featureMatrix) <- unique(cuts)

    rankMatrix <- apply(featureMatrix, 2, function(x) {rank(x, na.last = FALSE, ties.method='min')})
    rankMatrix <- nrow(featureMatrix) - rankMatrix

    featureMeans <- rowMeans(rankMatrix)
    featureSds <- apply(rankMatrix, 1, sd)
    rv <- data.frame(rankMean = featureMeans, rankSds = featureSds)
    rv <- rv[order(rv$rankMean), ]
    
    return(rv)
}

#' Rank features using mutual information
#'   This only seems to be the method in [R] that's reasonably quick
#'   The maximum score has the highest importance
#'
#' @param dataMat  The attribute matrix
#' @param classIdx  The index for the class variable (DEAFULT: nrow(dataMat))
#' @export
#' @usage featureList <- useMutualInformation(dataMat, classIdx=ncol(dataMat))
useMutualInformation <- function(dataMat, classVector) { ##=ncol(dataMat)) {
    classVector <- as.factor(classVector) ##dataMat[,classIdx])
    miVect <- rep(NA, ncol(dataMat))
    names(miVect) <- colnames(dataMat)
    for ( i in 1:ncol(dataMat)) {
        featureName <- colnames(dataMat)[i]
        if ( i != classIdx ) {
            curVector <- as.factor(dataMat[,i])
            miVect[[featureName]] <- mi.plugin(rbind(curVector, classVector))
        }
    }
                                        #miVect <- sort(miVect)
    return(miVect)
}
