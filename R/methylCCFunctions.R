#' @title pickTargetPositions
#'
#' @description Pick probes from target data using the indices in dmpRegions.
#'
#' @param targetgRanges add more here.
#' @param targetObject an optional argument which contains the meta-data for targetRanges. If
#' targetRanges already contains the meta-data, do not need to supply targetObject.
#' @param targetCVG add more here.
#' @param dmpRegions add more here.
#'
#' @import GenomicRanges
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise_all funs
#'
#' @name pickTargetPositions
#' @aliases pickTargetPositions
#'
pickTargetPositions <- function(targetgRanges, targetObject = NULL,
                                targetCVG = NULL, dmpRegions)
{
    if(!is.null(targetObject)){
        mcols(targetgRanges) <- targetObject
    }

    if(is(dmpRegions, "GRanges")){
        traingRanges = dmpRegions
    } else {
        traingRanges = makeGRangesFromDataFrame(dmpRegions, keep.extra.columns=TRUE)
    }

    traingRanges <- subsetByOverlaps(traingRanges, targetgRanges) # Zmat regions that overlap with target data
    subjectMtch <- subjectHits(findOverlaps(targetgRanges, traingRanges)) # regions from Zmat (CpG level)
    queryMtch <- queryHits(findOverlaps(targetgRanges, traingRanges)) # CpGs from target data (CpG level)
    intersectgRanges <- targetgRanges[queryHits(findOverlaps(targetgRanges, traingRanges)), ] # (CpG level)

    mcols(intersectgRanges)$subjectMtch <- subjectMtch
    mcols(intersectgRanges)$queryMtch <- queryMtch
    metaDat = as.data.frame(mcols(intersectgRanges)) # (CpG level, in order from chr1)

    if(!is.null(targetCVG)){
        targetgRangesCVG <- targetgRanges
        mcols(targetgRangesCVG) <- targetCVG

        intersectCVG <- targetgRangesCVG[queryMtch, ] # (CpG level, in order from chr1)
        mcols(intersectCVG)$subjectMtch <- subjectMtch
        mcols(intersectCVG)$queryMtch <- queryMtch
        metaDatCVG = as.data.frame(mcols(intersectCVG))

        # Sum M values and cvg values for all CpGs in each region (subject)
        out.M = metaDat %>% group_by(subjectMtch) %>% summarise_all(funs(sum(.,na.rm = TRUE)))
        out.CVG = metaDatCVG %>% group_by(subjectMtch) %>% summarise_all(funs(sum(.,na.rm = TRUE)))
        dmpGRanges <- traingRanges
        mcols(dmpGRanges) <- out.M[, !(colnames(out.M) %in% c("queryMtch", "subjectMtch"))] /
            out.CVG[, !(colnames(out.CVG) %in% c("queryMtch", "subjectMtch"))]
        mcols(intersectgRanges) <- metaDat[, !(colnames(metaDat) %in% c("queryMtch", "subjectMtch"))] /
            metaDatCVG[, !(colnames(metaDatCVG) %in% c("queryMtch", "subjectMtch"))]
    } else {
        # Average M values for all CpGs in each region (subject)
        out = metaDat %>% group_by(subjectMtch) %>% summarise_all(funs(mean(., na.rm = TRUE)))
        dmpGRanges <- traingRanges
        mcols(dmpGRanges) <- out[, !(colnames(out) %in% c("queryMtch", "subjectMtch"))]
        mcols(intersectgRanges) <- metaDat[, !(colnames(metaDat) %in% c("queryMtch", "subjectMtch"))]

    }

    list("probeGRanges" = intersectgRanges,
         "dmpGRanges" = dmpGRanges, "avgDMRDat" = mcols(dmpGRanges))
}


#' @title oneVec
#'
#' @description One dimensional vector of ones
#'
#' @param R add more here
#'
#' @name oneVec
#' @aliases oneVec
oneVec <- function(R){ matrix(rep(1, R), ncol = 1) }

generateSimFrame <- function(ids, pi = NULL, alpha0 = NULL, alpha1 = NULL,
                             sig0 = NULL, sig1 = NULL, likFun = NULL){
    mleNames <- c(ids, "alpha0", "alpha1",
                  "sig0", "sig1", "sig", "likFun")
    mleValues <- matrix(NA, nrow = 1, ncol = length(mleNames))
    colnames(mleValues) <- mleNames

    if(is.null(pi)){
        wts = runif(length(ids))
        pi <- wts / sum(wts) }
    if(is.null(alpha0)){ alpha0 <- runif(1,0.10, 0.3) }
    if(is.null(alpha1)){ alpha1 <- runif(1,0.70, 0.9) }
    if(is.null(sig0)){ sig0 <- runif(1, 0.01, 0.05)^2 }
    if(is.null(sig1)){ sig1 <- runif(1, 0.01, 0.05)^2 }
    sig <- runif(1, 0.01, 0.05)^2
    likFun <- 0
    mleValues[1,] <- c(pi, alpha0, alpha1, sig0, sig1, sig, likFun)
    mleValues <- as.data.frame(mleValues)
    return(mleValues)
}

WFun <- function(ZMat, est, n){
    W0 = c((1-ZMat) %*% t(est[seq_len(n)]))
    W1 = c(ZMat %*% t(est[seq_len(n)]))
    data.frame("W0" = W0, "W1" = W1)
}


deltaEstep <- function(ZMat, est, YMat, R, n, methStatus = 0)
{
    W <- WFun(ZMat=ZMat, est=est, n=n)
    mat11 = solve(diag(c((W$W0^2*est$sig0 + W$W1^2*est$sig1 + est$sig)), R))

    if(methStatus == 0){
        mat12.0 = diag(c(W$W0 * est$sig0), R)
        mat22.0 = diag(c(est$sig0), R)
        mom1 = est$alpha0*oneVec(R) +
            (t(mat12.0) %*% mat11 %*%
                 (YMat - est$alpha0*W$W0 - est$alpha1*W$W1))
        mom2 = mom1^2 + diag(mat22.0 - t(mat12.0) %*% mat11 %*% mat12.0)
    }
    if(methStatus == 1){
        mat12.1 = diag(c(W$W1 * est$sig1), R)
        mat22.1 = diag(c(est$sig1), R)
        mom1 = est$alpha1*oneVec(R) +
            (t(mat12.1) %*% mat11 %*%
                 (YMat - est$alpha0*W$W0 - est$alpha1*W$W1))
        mom2 = mom1^2 + diag(mat22.1 - t(mat12.1) %*% mat11 %*% mat12.1)
    }
    list("mom1" = mom1, "mom2" = mom2)
}

#' @title Maximization step
#'
#' @description Maximization step in EM Algorithm
#'
#' @param ZMat add more here
#' @param est add more here
#' @param YMat add more here
#' @param R add more here
#' @param n add more here
#' @param ids add more here
#' @param del0 add more here
#' @param del1 add more here
#'
#' @importFrom quadprog solve.QP
#' @name MStep
#' @aliases MStep
MStep <- function(ZMat, est, YMat, R, n, ids, del0, del1)
{
    W <- WFun(ZMat=ZMat, est=est, n=n)
    newMLEs <- generateSimFrame(ids = ids)
    newMLEs$alpha0 <- mean(del0$mom1)
    newMLEs$alpha1 <- mean(del1$mom1)
    newMLEs$sig0 <- mean(del0$mom2) - (newMLEs$alpha0)^2
    newMLEs$sig1 <- mean(del1$mom2) - (newMLEs$alpha1)^2
    newMLEs$sig <- est$sig

    xnew <- (sweep((1-ZMat), 1, del0$mom1, FUN = "*") +
                 sweep(ZMat, 1, del1$mom1, FUN = "*"))
    newMLEs[seq_len(n)] <-
        solve.QP(Dmat = (t(xnew)%*%xnew),
                 dvec = t(xnew)%*%YMat, Amat = cbind(rep(1,n), diag(n)),
                 bvec = c(1, rep(0, n)), meq = 1)$sol
    newMLEs$likFun <- -(R/2)*log(newMLEs$sig) -
        (1/(2*newMLEs$sig))*sum((YMat - xnew%*%t(newMLEs[seq_len(n)]))^2) -
        (R/2) * log(newMLEs$sig0) -
        (1/(2*newMLEs$sig0)) * sum((del0$mom1 - newMLEs$alpha0)^2) -
        (R/2) * log(newMLEs$sig1) -
        (1/(2*newMLEs$sig1)) * sum((del1$mom1 - newMLEs$alpha1)^2)

    return(as.data.frame(newMLEs))
}

#' @title preProcessEstCountsSeq
#'
#' @description add more here
#'
#' @param object an object can be an \code{RGChannelSet} or
#' \code{BSseq} object. The \code{object} can also be a
#' data frame or matrix with observations (e.g. CpGs)
#' on the rows and samples as the columns.
#' @param regionMat If the \code{object} is not an
#' \code{RGChannelSet} or \code{BSseq} object, then
#' the user must provide \code{regionMat}. 
#' This must have the same number of
#' rows as \code{object} (number of regions or CpGs)
#' and one column for each cell type. The values in
#' \code{regionMat} should be either 0 or 1 representing
#' regions (or CpGs) that are either methylated (1) or
#' not methylated (0).
#' @param verbose TRUE/FALSE argument specifying if verbose
#' messages should be returned or not. Default is TRUE.
#' @param initParamMethod add more here
#'
#' @import GenomicRanges
#' @importFrom Biobase ExpressionSet rowMin pData sampleNames
#' @importFrom minfi RGChannelSet preprocessIllumina mapToGenome getBeta
#' @importFrom bsseq BSseq getBSseq
#' @name preProcessEstCountsSeq
#' @aliases preProcessEstCountsSeq
preProcessEstCountsSeq <- function(object, regionMat, verbose=TRUE,
                                   initParamMethod)
{
    # load cell type specific DMRs
    suppressMessages(data(celltypeSpecificDMRegions))

    if(initParamMethod == "knownRegions"){
        if(!(is(object, "RGChannelSet") || is(object, "BSseq"))){
            stop("If initParamMethod is set to 'knownRegions', then 'object'\n
                 must either be an 'RGChannelSet' or be a 'BSseq' object.")
        } else {
            suppressMessages(data(offMethRegions))
            suppressMessages(data(onMethRegions))
            grd0 = makeGRangesFromDataFrame(offMethRegions)
            grd1 = makeGRangesFromDataFrame(onMethRegions)
        }
    }

    if(is(object, "RGChannelSet") || is(object, "BSseq")) {
        if(verbose){
            if(!is.null(regionMat)){
                message("[estimateCC] Ignoring regionMat.") }
        }

        if(is(object, "RGChannelSet")) {
            if(verbose){
                message("[estimateCC] Extracting RGChannelSet data.") }
            Mset <- preprocessIllumina(updateObject(object))
            Mset <- mapToGenome(Mset, mergeManifest = FALSE)
            pd = pData(Mset)
            grObject <- granges(Mset)
            betaValues = getBeta(Mset, type = "Illumina")
            rm(Mset)
            mcols(grObject) <- betaValues
            rm(betaValues)
            
            keepDMRs = pickTargetPositions(targetgRanges=grObject,
                                           dmpRegions=celltypeSpecificDMRegions)$dmpGRanges
            zmatSub <- subsetByOverlaps(celltypeSpecificDMRegions, keepDMRs, type = "equal")

            if(initParamMethod == "knownRegions"){
                y_d0 <- as.data.frame(pickTargetPositions(targetgRanges = grObject,
                                                dmpRegions = grd0)$avgDMRDat)
                a0init <- colMeans(y_d0)
                sig0init <- apply(y_d0, 2, var)
                y_d1 <- as.data.frame(pickTargetPositions(targetgRanges = grObject,
                                        dmpRegions = grd1)$avgDMRDat)
                a1init <- colMeans(y_d1)
                sig1init <- apply(y_d1, 2, var)
            }
        }

        if(is(object, "BSseq")){
            if(verbose){
                message("[estimateCC] Extracting BSSeq data.") }
            object <- updateObject(object)
            cvg_targetbs <- getBSseq(object, type = "Cov")
            M_targetbs <- getBSseq(object, type = "M")
            cvg_targetbs = getBSseq(object, type = "Cov")
            grObject <- granges(object)
            keepDMRs = pickTargetPositions(targetgRanges = grObject,
                            targetObject = M_targetbs, targetCVG = cvg_targetbs,
                            dmpRegions = celltypeSpecificDMRegions)$dmpGRanges
            grObject = pickTargetPositions(targetgRanges = grObject,
                                           targetObject = M_targetbs,
                                           targetCVG = cvg_targetbs,
                                           dmpRegions = grObject)$dmpGRanges
            zmatSub <- subsetByOverlaps(celltypeSpecificDMRegions, keepDMRs,
                                        type = "equal")
            if(verbose){
                mes <- "[estimateCC] BSseq object contained
                %s out of %s celltype-specific regions."
                message(sprintf(mes, length(keepDMRs), length(celltypeSpecificDMRegions)))
            }

            if(initParamMethod == "knownRegions"){
                y_d0 <- as.data.frame(pickTargetPositions(
                                            targetgRanges = grObject,
                                            targetObject = M_targetbs,
                                            targetCVG = cvg_targetbs,
                                            dmpRegions = grd0)$avgDMRDat)
                a0init <- colMeans(y_d0)
                sig0init <- apply(y_d0, 2, var)
                y_d1 <- as.data.frame(pickTargetPositions(
                                            targetgRanges = grObject,
                                            targetObject = M_targetbs,
                                            targetCVG = cvg_targetbs,
                                            dmpRegions = grd1)$avgDMRDat)
                a1init <- colMeans(y_d1)
                sig1init <- apply(y_d1, 2, var)
            }
        }

        ymat <- mcols(keepDMRs); colnames(ymat) <- sampleNames(object)
        nS <- ncol(mcols(keepDMRs))
        ids <- colnames(mcols(zmatSub))
        zmat <- as.matrix(as.data.frame(mcols(zmatSub)))

        } else {

            if(verbose){
                if(is.null(regionMat)){
                    stop("User must provide regionMat if object\n
                          is not an RGChannelSet or BSseq object.") }
            }
            ymat <- t(t(object))
            nS <- ncol(object)
            ids <- colnames(regionMat)
            zmat <- as.matrix(regionMat)
            grObject <- GRanges()
            keepDMRs <- GRanges()
        }
    if(initParamMethod == "knownRegions"){
        list(ymat = ymat, nS = nS, ids = ids, zmat = zmat,
             a0init = a0init, a1init = a1init,
             sig0init = sig0init, sig1init = sig1init,
             grObject = grObject, keepDMRs = keepDMRs)
    }

    list(ymat = ymat, nS = nS, ids = ids, zmat = zmat,
         grObject = grObject, keepDMRs = keepDMRs)
}



#' @title estimateCC
#'
#' @description Estimate cell composition from seq data
#' @name estimateCC
#' @aliases estimateCC
#'
#' @param object an object can be an \code{RGChannelSet} or
#' \code{BSseq} object. The \code{object} can also be a
#' data frame or matrix with observations (e.g. CpGs)
#' on the rows and samples as the columns.
#' @param regionMat If the \code{object} is not an
#' \code{RGChannelSet} or \code{BSseq} object, then
#' the user must provide \code{regionMat}. 
#' This must have the same number of
#' rows as \code{object} (number of regions or CpGs)
#' and one column for each cell type. The values in
#' \code{regionMat} should be either 0 or 1 representing
#' regions (or CpGs) that are either methylated (1) or
#' not methylated (0).
#' @param verbose TRUE/FALSE argument specifying if verbose
#' messages should be returned or not. Default is TRUE.
#' @param epsilon Threshold for Expectation-Maximization
#' (EM) algorithm
#' @param maxIter Maximum number of iterations for EM
#' algorithm
#' @param initParamMethod add more here
#' @param includeSE TRUE/FALSE argument specifying if
#' standard error estimates should be included. Default
#' is FALSE The parameter \code{B} must be set to a value
#' greater than 0 if \code{includeSE == TRUE}
#' @param B Number of bootstrap samples to simulate to
#' calculate SE estimates.
#' @param a0init add more here.
#' @param a1init add more here.
#' @param sig0init add more here.
#' @param sig1init add more here.
#'
#' @importFrom mclust Mclust
#' @importFrom quadprog solve.QP
#' @export
#'
estimateCC <- function(object, regionMat = NULL, verbose = TRUE,
                                  epsilon = 0.01, maxIter = 50,
                                  initParamMethod = "random", includeSE = FALSE,
                                  B = 0,
                                  a0init = NULL, a1init = NULL,
                                  sig0init = NULL, sig1init = NULL)
{

    if(initParamMethod == "fixed"){
        if(is.null(a1init) || is.null(a0init) ||
           is.null(sig0init) || is.null(sig1init)){
            stop("If initParamMethod is set to 'fixed', then user must\n
                 provide initial values for a0init, a1init, sig0init, sig1init.")
        }
    }

    dat = preProcessEstCountsSeq(object, regionMat=regionMat, verbose=verbose,
                                initParamMethod=initParamMethod)
    ymat <- dat$ymat
    nS <- dat$nS
    ids <- dat$ids
    zmat <- dat$zmat

    results <- new("estimateCC")
    results@YMat <- as.data.frame(ymat)
    results@ZMat <- as.data.frame(zmat)
    results@includeSE <- includeSE
    results@grObject <- dat$grObject
    results@keepDMRs <- dat$keepDMRs

    if(initParamMethod == "knownRegions"){
        a0init = dat$a0init; sig0init = dat$sig0init
        a1init = dat$a1init; sig1init = dat$sig1init
    }

    n = ncol(zmat)
    decoEngine <- function(YMat, ZMat){

        finalMLEs <- NULL
        count = 0
        repeat{
            # E-Step
            del0 <- deltaEstep(ZMat=ZMat, est=currentMLEs,
                               YMat=YMat, R=R, n=n, methStatus=0)
            del1 <- deltaEstep(ZMat=ZMat, est=currentMLEs,
                               YMat=YMat, R=R, n=n, methStatus=1)

            # M-Step
            updateMLEs <- MStep(ZMat=ZMat, est=currentMLEs,
                                YMat=YMat, R=R, n=n, ids=ids, del0, del1)

            if(abs(updateMLEs$likFun - currentMLEs$likFun) >= epsilon){
                finalMLEs <- rbind(finalMLEs, currentMLEs)
                currentMLEs <- updateMLEs
                count = count + 1
                if(count > maxIter){
                    break
                }
            } else {
                break
            }
        }
        return(finalMLEs)
    }
    initializeMLEs <- function(initParamMethod, Ys, Zs, ids, n,
                               a0init, a1init, sig0init, sig1init){

        if(initParamMethod == "random"){
            currentMLEs <- generateSimFrame(ids = ids)
        }

        if(initParamMethod=="fixed" || initParamMethod=="knownRegions"){
            currentMLEs <- generateSimFrame(ids = ids,
                                            alpha0 = a0init[i], alpha1 = a1init[i],
                                            sig0 = sig0init[i], sig1 = sig1init[i])
        }

        if(initParamMethod == "clust"){
            out = Mclust(Ys[,i], G = 2, modelNames = "V")
            mu <- out$parameters$mean
            currentMLEs <- generateSimFrame(ids = ids,
                                            alpha0 = min(mu[1], 0.49),
                                            alpha1 = max(0.5, mu[2]))
        }

        # Initial guess for cell composition weights
        x0 = (1-Zs) * currentMLEs$alpha0 + (Zs) * currentMLEs$alpha1
        currentMLEs[1,seq_len(n)] <-
            round(solve.QP(Dmat = (t(x0)%*%x0),
                           dvec = t(x0)%*%Ys,
                           Amat = cbind(rep(1,n), diag(n)),
                           bvec = c(1, rep(0, n)), meq = 1)$sol,3)
        currentMLEs$sig <-
            mean((Ys - (x0 %*% t(currentMLEs[1, seq_len(n)])))^2)

        return(currentMLEs)
    }

    cellcounts = data.frame(array(NA, dim=c(nS,n)))
    mles = data.frame(array(NA, dim=c(nS, (n+6))))
    mleList <- vector("list", nS)
    if(includeSE){
        countsMean = array(NA, dim=c(nS, n))
        countsSE = array(NA, dim=c(nS, n))
        mleSE = array(NA, dim=c(nS, (n+5)))
        mleSEboot = array(NA, dim=c(nS, (n+6), B))
    }

    nRegions = array(NA, dim = nS)
    for(i in seq_len(nS)){
        keepme <- c(!is.na(ymat[,i]))
        YMatSub <- t(t(ymat[keepme,i]))
        ZMatSub <- zmat[keepme,]
        R = nrow(ZMatSub)
        nRegions[i] <- sum(keepme)

        ## Include verbose messages about parameter estimation
        if(verbose){
            mes <- "[estimateCC] Starting parameter estimation for sample %s
                    using %s regions."
            message(sprintf(mes, i, nRegions[i]))
        }

        # Initialize MLEs
        currentMLEs <- initializeMLEs(initParamMethod=initParamMethod,
                                      Ys=YMatSub, Zs=ZMatSub, ids=ids, n=n,
                                      a0init=a0init, a1init=a1init,
                                      sig0init=sig0init, sig1init=sig1init)

        # Run EM algorithm
        finalMLEs <- decoEngine(YMat = YMatSub, ZMat = ZMatSub)

        if(verbose){
            mes <- "[estimateCC] Parameter estimates obtained in %s iterations for sample %s."
            message(sprintf(mes, nrow(finalMLEs) - 1, i))
        }

        cellcounts[i,] = finalMLEs[nrow(finalMLEs), seq_len(n)]
        mles[i,] = finalMLEs[nrow(finalMLEs), ]
        mleList[[i]] = finalMLEs


        if(includeSE){
            if(B<=0){
                stop("B must be greater than 0 for boostrap SE estimation.")
            }

            if(verbose){
                mes <- "[estimateCC] Performing bootstrap SE estimation using %s random samples with replacement.\n
                        This may take a few minutes."
                message(sprintf(mes, B))
            }

            b = 0
            repeat{
                shuffID = sample(seq_len(length(YMatSub)), replace = TRUE)
                Yboot <- YMatSub[shuffID,]
                Zboot <- ZMatSub[shuffID,]

                tmp <- try(initializeMLEs(initParamMethod=initParamMethod,
                                          Ys=Yboot, Zs=Zboot, ids=ids, n=n,
                                          a0init=a0init, a1init=a1init,
                                          sig0init=sig0init, sig1init=sig1init), silent = TRUE)
                if('try-error' %in% class(tmp)){
                    next
                } else {
                    currentMLEs <- tmp
                }

                # Run EM algorithm
                finalMLEs <- decoEngine(YMat = Yboot, ZMat = Zboot)

                mleSEboot[i,,b] = as.numeric(finalMLEs[nrow(finalMLEs), ])
                if(b==B){ break } else { b = b + 1}
            }
            countsMean[i,] <- apply(mleSEboot[i, seq_len(n),], 1, mean)
            countsSE[i,] <- apply(mleSEboot[i, seq_len(n),], 1, sd)
            mleSE[i,] <- apply(mleSEboot[i,-(n+6),], 1, sd)
        }
    }

    colnames(cellcounts) <- ids
    colnames(mles) <- colnames(finalMLEs)
    results@cellcounts <- cellcounts
    results@mles <- mles
    results@mleList <- mleList
    results@summary <- list("class" = class(object),
                            "nSamples" = nS, "cellTypes" = ids,
                            "sampleNames" = colnames(ymat),
                            "initParamMethod" = initParamMethod,
                            "nRegionsPerSample" = nRegions)

    if(includeSE){
        colnames(countsMean) <- ids
        colnames(countsSE) <- ids
        colnames(mleSE) <- colnames(finalMLEs)[-(n+6)]
        results@countsMean <- countsMean
        results@countsSE <- countsSE
        results@mleSE <- mleSE
    } else {
        results@countsMean <- array(NA)
        results@countsSE <- array(NA)
        results@mleSE <- array(NA)
    }

    return(results)
}


