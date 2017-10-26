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
    if(is(dmpRegions, "GRanges")){
        traingRanges = dmpRegions
    } else {
        traingRanges = makeGRangesFromDataFrame(dmpRegions, keep.extra.columns=TRUE)
    }

    traingRanges <- subsetByOverlaps(traingRanges, targetgRanges) # Zmat regions that overlap with target data
    subjectMtch <- subjectHits(findOverlaps(targetgRanges, traingRanges)) # regions from Zmat (CpG level)
    queryMtch <- queryHits(findOverlaps(targetgRanges, traingRanges)) # CpGs from target data (CpG level)
    intersectgRanges <- targetgRanges[queryMtch, ] # (CpG level)
    if(!is.null(targetObject)){
      mcols(intersectgRanges) <- targetObject[queryMtch, ]
    }

    mcols(intersectgRanges)$subjectMtch <- subjectMtch
    mcols(intersectgRanges)$queryMtch <- queryMtch
    metaDat = as.data.frame(mcols(intersectgRanges)) # (CpG level, in order from chr1)

    if(!is.null(targetCVG)){
        intersectCVG <- targetgRanges[queryMtch, ] # (CpG level, in order from chr1)
        mcols(intersectCVG) <- targetCVG[queryMtch, ]
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


#' @title initialize_theta
#'
#' @description Creates a container with initial theta
#' parameter estimates
#'
#' @param n Number of samples
#' @param K Number of cell types
#' @param alpha0 Default NULL. Initial mean methylation level in unmethylated regions
#' @param alpha1 Default NULL. Initial mean methylation level in methylated regions
#' @param sig0 Default NULL. Initial var methylation level in unmethylated regions
#' @param sig1 Default NULL. Initial var methylation level in methylated regions
#' @param tau Default NULL. Initial var for measurement error
#'
#' @name initialize_theta
#' @aliases initialize_theta
initialize_theta <- function(n, K, alpha0 = NULL, alpha1 = NULL,
                       sig0 = NULL, sig1 = NULL, tau = NULL){
    if(is.null(alpha0)){ alpha0 <- runif(1,0.10, 0.3) }
    if(is.null(alpha1)){ alpha1 <- runif(1,0.70, 0.9) }
    if(is.null(sig0)){ sig0 <- runif(1, 0.01, 0.05)^2 }
    if(is.null(sig1)){ sig1 <- runif(1, 0.01, 0.05)^2 }
    if(is.null(tau)){ tau <- runif(1, 0.001, 0.01)^2 }
    return(data.frame("alpha0" = alpha0, "alpha1" = alpha1,
                      "sig0" = sig0, "sig1" = sig1, "tau" = tau,
                      "likFun" = 0))
}

#' @title Helper function to take the product of Z and
#' cell composition estimates
#'
#' @description Helper function which is the product of
#' Z and pi.mle
#'
#' @param ZMat Cell type specific regions of dimension R x K
#' @param pi.mle cell composition MLE estimates
#'
#' @name WFun
#' @aliases WFun
WFun <- function(ZMat, pi.mle){
    W0 = (1-ZMat) %*% t(pi.mle)
    W1 = ZMat %*% t(pi.mle)
    list("W0" = W0, "W1" = W1)
}

#' @title Expectation step
#'
#' @description Expectation step in EM algorithm for methylCC
#'
#' @param ZMat Cell type specific regions of dimension R x K
#' @param pi.mle cell composition MLE estimates of dimension K x n
#' @param theta other parameter estimates in EM algorithm
#' @param YMat observed methylation levels in samples provided by user
#' of dimension R x n
#' @param methStatus Indicator function corresponding to regions that
#' are unmethylated (methStatus=0) or methylated (methStatus=1)
#'
#' @name methylCC_EStep
#' @aliases methylCC_EStep
#'
methylCC_Estep <- function(ZMat, pi.mle, theta, YMat, methStatus = 0)
{
    R <- nrow(ZMat)
    n <- ncol(YMat)
    W <- WFun(ZMat=ZMat, pi.mle=pi.mle)

    # dim(mat11.inv) = (R*n x R*n) in block diagonal form
      mat11.inv <- lapply(seq_len(nrow(ZMat)), function(r){
              # take inverse of individual blocks
            solve( t(t(W$W0[r,])) %*% t(W$W0[r,]) * theta$sig0 +
            t(t(W$W1[r,])) %*% t(W$W1[r,]) * theta$sig1 +
            diag(theta$tau, n) )    })

    if(methStatus == 0){
        # dim(mat12.0) = (R*n x R) in block diagonal form
        mat12.0 <- lapply(seq_len(nrow(ZMat)), function(r){
          t(t(W$W0[r,])) * theta$sig0 })

        # dim(mat22.0) = (R x R)
        mat22.0 <- diag(theta$sig0, R)

        # dim(mom1) = dim(mom2) = (R x 1)
        fast.mat.mom1 <- t(t(sapply(seq_len(nrow(ZMat)), function(r){
          t(mat12.0[[r]]) %*% mat11.inv[[r]] %*%
            (t(t(YMat[r,])) - theta$alpha0*t(t(W$W0[r,])) - theta$alpha1*t(t(W$W1[r,])) )
        })))
        fast.mat.mom2 <- diag(sapply(seq_len(nrow(ZMat)), function(r){
          t(mat12.0[[r]]) %*% mat11.inv[[r]] %*% mat12.0[[r]]
        }), R)

        mom1 = theta$alpha0 * matrix(rep(1, R), ncol=1) + fast.mat.mom1
        mom2 = mom1^2 + t(t(diag(mat22.0 - fast.mat.mom2 )))
    }
    if(methStatus == 1){
      # dim(mat12.1) = (R*n x R) in block diagonal form
      mat12.1 <- lapply(seq_len(nrow(ZMat)), function(r){
        t(t(W$W1[r,])) * theta$sig1 })

      # dim(mat22.1) = (R x R)
      mat22.1 <- diag(theta$sig1, R)

      # dim(mom1) = dim(mom2) = (R x 1)
      fast.mat.mom1 <- t(t(sapply(seq_len(nrow(ZMat)), function(r){
        t(mat12.1[[r]]) %*% mat11.inv[[r]] %*%
          (t(t(YMat[r,])) - theta$alpha0*t(t(W$W0[r,])) - theta$alpha1*t(t(W$W1[r,])) )
      })))
      fast.mat.mom2 <- diag(sapply(seq_len(nrow(ZMat)), function(r){
        t(mat12.1[[r]]) %*% mat11.inv[[r]] %*% mat12.1[[r]]
      }), R)

      mom1 = theta$alpha1 * matrix(rep(1, R), ncol=1) + fast.mat.mom1
      mom2 = mom1^2 + t(t(diag(mat22.1 - fast.mat.mom2 )))
    }
    list("mom1" = mom1, "mom2" = mom2)
}

#' @title Maximization step
#'
#' @description Maximization step in EM Algorithm for methylCC
#'
#' @param ZMat Cell type specific regions of dimension R x K
#' @param pi.mle cell composition MLE estimates of dimension K x n
#' @param theta other parameter estimates in EM algorithm
#' @param YMat observed methylation levels in samples provided by user
#' of dimension R x n
#' @param EStep0 Results from expectation step for unmethylated regions
#' @param EStep1 Results from expectation step for methylated regions
#'
#' @importFrom quadprog solve.QP
#' @name methylCC_MStep
#' @aliases methylCC_MStep
methylCC_MStep <- function(ZMat, pi.mle, theta, YMat, EStep0, EStep1)
{
  R <- dim(ZMat)[1]
  K <- dim(ZMat)[2]
  n <- ncol(YMat)
  new.theta <- data.frame("alpha0" = mean(EStep0$mom1),
                          "alpha1" = mean(EStep1$mom1),
                          "sig0" = mean(EStep0$mom2) - mean(EStep0$mom1)^2,
                          "sig1" = mean(EStep1$mom2) - mean(EStep1$mom1)^2)
  x0 <- (sweep((1-ZMat), 1, EStep0$mom1, FUN = "*") +
         sweep(ZMat,     1, EStep1$mom1, FUN = "*"))
  new.pi.mle <- abs(t(apply(YMat, 2, function(ys){
    solve.QP(Dmat = (t(x0)%*%x0),
             dvec = t(x0)%*%t(t(ys)), Amat = cbind(rep(1,K), diag(K)),
             bvec = c(1, rep(0, K)), meq = 1)$sol })))
  tau.est <- (1/(R*n)) * sum((YMat - x0 %*% t(new.pi.mle))^2)
  new.theta <- cbind(new.theta, "tau" = tau.est,
                     "likFun" = -((R*n)/2)*log(2*pi*tau.est) -
                       (1/(2*tau.est)) * sum((YMat - x0 %*% t(new.pi.mle))^2) -
                       ((R*n)/2) * log(2*pi*new.theta$sig0) -
                       (n/(2*new.theta$sig0)) * sum((EStep0$mom1 - new.theta$alpha0)^2) -
                       ((R*n)/2) * log(2*pi*new.theta$sig1) -
                       (n/(2*new.theta$sig1)) * sum((EStep1$mom1 - new.theta$alpha1)^2))

  return(list("theta" = new.theta, "pi.mle" = new.pi.mle))
}

#' @title preProcessEstCountsSeq
#'
#' @description add more here
#'
#' @param object an object can be an \code{RGChannelSet} or
#' \code{BSseq} object. The \code{object} can also be a
#' data frame or matrix with observations (e.g. CpGs)
#' on the rows and whole blood samples as the columns.
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
#' @param initParamMethod method to initialize parameter estimates.
#' Choose between "random" (randomly sample) or "knownRegions"
#' (uses unmethyalted and methylated regions that were identified
#' based on Reinus et al. (2012) cell sorted data.).
#' Defaults to "random".
#'
#' @import GenomicRanges
#' @importFrom Biobase ExpressionSet rowMin pData sampleNames
#' @importFrom minfi RGChannelSet preprocessIllumina mapToGenome getBeta
#' @importFrom bsseq BSseq getBSseq
#' @name preProcessEstCountsSeq
#' @aliases preProcessEstCountsSeq
# dat <- preProcessEstCountsSeq(BS, regionMat = NULL, initParamMethod = "random")
preProcessEstCountsSeq <- function(object, regionMat, verbose=TRUE,
                                   initParamMethod = "random")
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
                message("[estimateCC] Warning: ignoring regionMat argument.") }
        }

        if(is(object, "RGChannelSet")) {
            if(verbose){
                message("[estimateCC] Preprocessing RGChannelSet object.") }
            Mset <- preprocessIllumina(object)
            Mset <- mapToGenome(Mset, mergeManifest = FALSE)
            pd = pData(Mset)
            grObject <- granges(Mset)
            betaValues = getBeta(Mset, type = "Illumina")
            rm(Mset)

            keepDMRs = pickTargetPositions(targetgRanges=grObject, targetObject = betaValues,
                                           dmpRegions=celltypeSpecificDMRegions)$dmpGRanges
            zmatSub <- subsetByOverlaps(celltypeSpecificDMRegions, keepDMRs, type = "equal")

            if(initParamMethod == "knownRegions"){
                y_d0 <- as.data.frame(pickTargetPositions(targetgRanges = grObject,
                                                          targetObject = betaValues,
                                                dmpRegions = grd0)$avgDMRDat)
                a0init <- mean(colMeans(y_d0))
                sig0init <- mean(apply(y_d0, 2, var))
                y_d1 <- as.data.frame(pickTargetPositions(targetgRanges = grObject,
                                                          targetObject = betaValues,
                                        dmpRegions = grd1)$avgDMRDat)
                a1init <- mean(colMeans(y_d1))
                sig1init <- mean(apply(y_d1, 2, var))
            }
            rm(betaValues)
        }

        if(is(object, "BSseq")){
            if(verbose){
                message("[estimateCC] Extracting BSSeq data.") }
            cvg_targetbs <- getBSseq(object, type = "Cov")
            M_targetbs <- getBSseq(object, type = "M")
            grObject <- granges(object)
            keepDMRs = pickTargetPositions(targetgRanges = grObject,
                            targetObject = M_targetbs, targetCVG = cvg_targetbs,
                            dmpRegions = celltypeSpecificDMRegions)$dmpGRanges
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
                a0init <- mean(colMeans(y_d0))
                sig0init <- mean(apply(y_d0, 2, var))
                y_d1 <- as.data.frame(pickTargetPositions(
                                            targetgRanges = grObject,
                                            targetObject = M_targetbs,
                                            targetCVG = cvg_targetbs,
                                            dmpRegions = grd1)$avgDMRDat)
                a1init <- mean(colMeans(y_d1))
                sig1init <- mean(apply(y_d1, 2, var))
            }
        }

        ymat <- mcols(keepDMRs); colnames(ymat) <- sampleNames(object)
        n <- ncol(mcols(keepDMRs))
        ids <- colnames(mcols(zmatSub))
        zmat <- as.matrix(as.data.frame(mcols(zmatSub)))

        } else {

            if(verbose){
                if(is.null(regionMat)){
                    stop("User must provide regionMat if object\n
                          is not an RGChannelSet or BSseq object.") }
            }
            ymat <- t(t(object))
            n <- ncol(object)
            ids <- colnames(regionMat)
            zmat <- as.matrix(regionMat)
            grObject <- GenomicRanges::GRanges()
            keepDMRs <- GenomicRanges::GRanges()
        }

  if(initParamMethod == "knownRegions"){
    list(ymat = ymat, n = n, ids = ids, zmat = zmat,
             a0init = a0init, a1init = a1init,
             sig0init = sig0init, sig1init = sig1init,
             grObject = grObject, keepDMRs = keepDMRs)
    } else {

    list(ymat = ymat, n = n, ids = ids, zmat = zmat,
         grObject = grObject, keepDMRs = keepDMRs)
    }
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
#' @param epsilon Threshold for EM algorithm to check
#' for convergence. Default is 0.01.
#' @param maxIter Maximum number of iterations for EM
#' algorithm. Default is 100 iterations.
#' @param initParamMethod method to initialize parameter estimates.
#' Choose between "random" (randomly sample) or "knownRegions"
#' (uses unmethyalted and methylated regions that were identified
#' based on Reinus et al. (2012) cell sorted data.).
#' Defaults to "random".
#' @param a0init Default NULL. Initial mean methylation level in unmethylated regions
#' @param a1init Default NULL. Initial mean methylation level in methylated regions
#' @param sig0init Default NULL. Initial var methylation level in unmethylated regions
#' @param sig1init Default NULL. Initial var methylation level in methylated regions
#' @param tauinit Default NULL. Initial var for measurement error
#'
#' @importFrom quadprog solve.QP
#' @export
#'
estimateCC <- function(object, regionMat = NULL, verbose = TRUE,
                       epsilon = 0.01, maxIter = 100, initParamMethod = "random",
                       a0init = NULL, a1init = NULL, sig0init = NULL,
                       sig1init = NULL, tauinit = NULL)
{

    if(!(initParamMethod %in% c("random", "knownRegions")) ){
      stop("The initParamMethod must be set to 'random' or 'knownRegions'.")
    }
    dat = preProcessEstCountsSeq(object, regionMat=regionMat, verbose=verbose,
                                 initParamMethod=initParamMethod)
    ymat <- as.matrix(dat$ymat[apply(dat$ymat, 1, function(x){all(!is.na(x))}), ])
    # ymat <- apply(dat$ymat, 2, function(x){
    #  x[is.na(x)] <- mean(x, na.rm = TRUE)
    #  x })
    n <- dat$n
    ids <- dat$ids
    zmat <- dat$zmat[apply(dat$ymat, 1, function(x){all(!is.na(x))}), ]
    K = ncol(zmat)
    R = nrow(zmat)

    results <- new("estimateCC")
    results@YMat <- as.data.frame(dat$ymat)
    results@ZMat <- as.data.frame(zmat)
    results@grObject <- dat$grObject
    results@keepDMRs <- dat$keepDMRs

    if(initParamMethod == "knownRegions"){
        a0init = dat$a0init; sig0init = dat$sig0init
        a1init = dat$a1init; sig1init = dat$sig1init
    }

    # Helper functions to initialize MLEs and run
    initializeMLEs <- function(initParamMethod, Ys, Zs, a0init, a1init,
                               sig0init, sig1init, tauinit){

      n <- ncol(Ys)
      K <- ncol(Zs)
      current.theta <- initialize_theta(n=n, K=K, alpha0=a0init,
                                        alpha1=a1init, sig0=sig0init,
                                        sig1=sig1init, tau=tauinit)

      # Initial guess for cell composition estimates
      x0 = (1-Zs) * current.theta$alpha0 + (Zs) * current.theta$alpha1
      current.pi.mle <- abs(t(apply(Ys, 2, function(ys){
        solve.QP(Dmat = (t(x0)%*%x0),
                 dvec = t(x0)%*%t(t(ys)), Amat = cbind(rep(1,K), diag(K)),
                 bvec = c(1, rep(0, K)), meq = 1)$sol })))

      return(list("pi.mle" = current.pi.mle, "theta" = current.theta))
    }

    methylCCEngine <- function(YMat, ZMat, current.pi.mle, current.theta,
                               epsilon, maxIter){

        final.theta <- NULL
        final.pi.mle <- NULL
        count = 0
        repeat{
            # E-Step
            EStep0 <- methylCC_Estep(ZMat=ZMat, pi.mle=current.pi.mle, theta=current.theta,
                               YMat=YMat, methStatus=0)
            EStep1 <- methylCC_Estep(ZMat=ZMat, pi.mle=current.pi.mle, theta=current.theta,
                               YMat=YMat, methStatus=1)

            # M-Step
            update.step <- methylCC_MStep(ZMat=ZMat, pi.mle=current.pi.mle, theta = current.theta,
                                YMat=YMat, EStep0=EStep0, EStep1=EStep1)

           if(abs(update.step$theta$likFun - current.theta$likFun) >= epsilon){
                final.pi.mle <- update.step$pi.mle
                final.theta <- rbind(final.theta, update.step$theta)
                current.pi.mle <- update.step$pi.mle
                current.theta <- update.step$theta
                count = count + 1

                if(count > maxIter){
                    break
                }
            } else {
                break
            }
        }
        return(list("pi.mle" = final.pi.mle, "theta" = final.theta))
    }

    ## Include verbose messages about parameter estimation
    if(verbose){
      mes <- "[estimateCC] Starting parameter estimation using %s regions."
      message(sprintf(mes, R))
    }

    # Initialize MLEs
    init.step <- initializeMLEs(initParamMethod=initParamMethod,
                                Ys=ymat, Zs=zmat, a0init=a0init,
                                a1init=a1init, sig0init=sig0init,
                                sig1init=sig1init, tauinit=tauinit)

    # Run EM algorithm
    finalMLEs <- methylCCEngine(YMat = ymat, ZMat = zmat,
                                current.pi.mle = init.step$pi.mle,
                                current.theta = init.step$theta,
                                epsilon=epsilon, maxIter = maxIter)

    if(verbose){
          mes <- "[estimateCC] Parameter estimates obtained in %s iterations."
          message(sprintf(mes, nrow(finalMLEs$theta) - 1))
    }

    cellcounts = as.data.frame(finalMLEs$pi.mle)
    colnames(cellcounts) <- ids
    results@cellcounts <- cellcounts
    results@theta <- finalMLEs$theta
    results@summary <- list("class" = class(object),
                            "nSamples" = n, "cellTypes" = ids,
                            "sampleNames" = colnames(ymat),
                            "initParamMethod" = initParamMethod,
                            "nRegions" = R)
    return(results)
}


