#' @title .splitit
#' 
#' @description helper function to split along a variable
#' 
#' @param x A vector
#' 
.splitit <- function(x) {
  split(seq(along = x), x)
}

#' @title .pick_target_positions
#'
#' @description Pick probes from target data using the indices in \code{dmp_regions}.
#'
#' @param target_granges add more here.
#' @param target_object an optional argument which contains the meta-data for targetRanges. If
#' targetRanges already contains the meta-data, do not need to supply \code{target_object}.
#' @param target_cvg add more here.
#' @param dmp_regions add more here.
#'
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps 
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise_all funs
#'
.pick_target_positions <- function(target_granges, target_object = NULL,
                                   target_cvg = NULL, dmp_regions)
{
  if(is(dmp_regions, "GRanges")){
    train_granges = dmp_regions
  } else {
    train_granges = makeGRangesFromDataFrame(dmp_regions, keep.extra.columns=TRUE)
  }
  
  train_granges <- subsetByOverlaps(train_granges, target_granges) # Zmat regions that overlap with target data
  subject_mtch <- subjectHits(findOverlaps(target_granges, train_granges)) # regions from Zmat (CpG level)
  query_mtch <- queryHits(findOverlaps(target_granges, train_granges)) # CpGs from target data (CpG level)
  intersect_granges <- target_granges[query_mtch, ] # (CpG level)
  if(!is.null(target_object)){
    mcols(intersect_granges) <- target_object[query_mtch, ]
  }
  
  mcols(intersect_granges)$subject_mtch <- subject_mtch
  mcols(intersect_granges)$query_mtch <- query_mtch
  metadata = as.data.frame(mcols(intersect_granges)) # (CpG level, in order from chr1)
  
  if(!is.null(target_cvg)){
    intersect_cvg <- target_granges[query_mtch, ] # (CpG level, in order from chr1)
    mcols(intersect_cvg) <- target_cvg[query_mtch, ]
    mcols(intersect_cvg)$subject_mtch <- subject_mtch
    mcols(intersect_cvg)$query_mtch <- query_mtch
    metadata_cvg = as.data.frame(mcols(intersect_cvg))
    
    # Sum M values and cvg values for all CpGs in each region (subject)
    out_m = metadata %>% group_by(subject_mtch) %>% summarise_all(funs(sum(.,na.rm = TRUE)))
    out_cvg = metadata_cvg %>% group_by(subject_mtch) %>% summarise_all(funs(sum(.,na.rm = TRUE)))
    dmp_granges <- train_granges
    mcols(dmp_granges) <- out_m[, !(colnames(out_m) %in% c("query_mtch", "subject_mtch"))] /
      out_cvg[, !(colnames(out_cvg) %in% c("query_mtch", "subject_mtch"))]
    mcols(intersect_granges) <- metadata[, !(colnames(metadata) %in% c("query_mtch", "subject_mtch"))] /
      metadata_cvg[, !(colnames(metadata_cvg) %in% c("query_mtch", "subject_mtch"))]
  } else {
    # Average M values for all CpGs in each region (subject)
    out = metadata %>% group_by(subject_mtch) %>% summarise_all(funs(mean(., na.rm = TRUE)))
    dmp_granges <- train_granges
    mcols(dmp_granges) <- out[, !(colnames(out) %in% c("query_mtch", "subject_mtch"))]
    mcols(intersect_granges) <- metadata[, !(colnames(metadata) %in% c("query_mtch", "subject_mtch"))]
    
  }
  
  list("probeGRanges" = intersect_granges,
       "dmp_granges" = dmp_granges, "avg_dmr_dat" = mcols(dmp_granges))
}


#' @title .preprocess_estimatecc
#'
#' @description add more here
#'
#' @param object an object can be an \code{RGChannelSet} or
#' \code{BSseq} object. The \code{object} can also be a
#' data frame or matrix with observations (e.g. CpGs)
#' on the rows and whole blood samples as the columns.
#' @param region_mat If the \code{object} is not an
#' \code{RGChannelSet} or \code{BSseq} object, then
#' the user must provide \code{region_mat}.
#' This must have the same number of
#' rows as \code{object} (number of regions or CpGs)
#' and one column for each cell type. The values in
#' \code{region_mat} should be either 0 or 1 representing
#' regions (or CpGs) that are either methylated (1) or
#' not methylated (0).
#' @param verbose TRUE/FALSE argument specifying if verbose
#' messages should be returned or not. Default is TRUE.
#' @param init_param_method method to initialize parameter estimates.
#' Choose between "random" (randomly sample) or "known_regions"
#' (uses unmethyalted and methylated regions that were identified
#' based on Reinus et al. (2012) cell sorted data.).
#' Defaults to "random".
#' @param celltype_specific_dmrs Add here.
#'
#' @import minfi
#' @import IlluminaHumanMethylation450kmanifest
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps 
#' @importFrom stats var
#' @importFrom utils data
#' @importFrom Biobase pData sampleNames
#' @importFrom bsseq BSseq getBSseq
#' 
# dat <- preprocess_estimatecc(BS, region_mat = NULL, init_param_method = "random")
.preprocess_estimatecc <- function(object, region_mat, verbose=TRUE,
                                    init_param_method = "random", 
                                    celltype_specific_dmrs = celltype_specific_dmrs)
{
  
  if(init_param_method == "known_regions"){
    if(!(is(object, "RGChannelSet") || 
         is(object, "GenomicMethylSet") || 
         is(object, "BSseq"))){
      stop("If init_param_method is set to 'known_regions', then 'object'\n
                 must be a 'RGChannelSet', 'GenomicMethylSet' or be a 'BSseq' object.")
    } else {
      suppressMessages(data(offMethRegions))
      suppressMessages(data(onMethRegions))
      grd0 = makeGRangesFromDataFrame(offMethRegions)
      grd1 = makeGRangesFromDataFrame(onMethRegions)
    }
  }
  
  if(is(object, "RGChannelSet") || 
     is(object, "GenomicMethylSet") || 
     is(object, "BSseq")) {
    
    if(verbose){
      if(!is.null(region_mat)){
        message("[estimatecc] Warning: ignoring region_mat argument.") }
    }
    
    if(is(object, "RGChannelSet")) {
      if(verbose){
        message("[estimatecc] Preprocessing RGChannelSet object.") }
      Mset <- preprocessIllumina(object)
      Mset <- mapToGenome(Mset, mergeManifest = FALSE)
      pd = pData(Mset)
      gr_object <- granges(Mset)
      beta_values = getBeta(Mset, type = "Illumina")
      rm(Mset)
    } 
    
    if(is(object, "GenomicMethylSet")) {
      pd = pData(object)
      gr_object <- granges(object)
      beta_values = getBeta(object, type = "Illumina")
    }
    
    if(is(object, "RGChannelSet") || 
       is(object, "GenomicMethylSet")) {       
      
      keep_dmrs = .pick_target_positions(target_granges=gr_object, target_object = beta_values,
                                         dmp_regions=celltype_specific_dmrs)$dmp_granges
      zmat_sub <- subsetByOverlaps(celltype_specific_dmrs, keep_dmrs, type = "equal")
      
      if(init_param_method == "known_regions"){
        y_d0 <- as.data.frame(.pick_target_positions(target_granges = gr_object,
                                                     target_object = beta_values,
                                                     dmp_regions = grd0)$avg_dmr_dat)
        a0init <- mean(colMeans(y_d0))
        sig0init <- mean(apply(y_d0, 2, var))
        y_d1 <- as.data.frame(.pick_target_positions(target_granges = gr_object,
                                                     target_object = beta_values,
                                                     dmp_regions = grd1)$avg_dmr_dat)
        a1init <- mean(colMeans(y_d1))
        sig1init <- mean(apply(y_d1, 2, var))
      }
      rm(beta_values)
    }
    
    if(is(object, "BSseq")){
      if(verbose){
        message("[estimatecc] Extracting BSSeq data.") }
      cvg_targetbs <- getBSseq(object, type = "Cov")
      M_targetbs <- getBSseq(object, type = "M")
      gr_object <- granges(object)
      keep_dmrs = .pick_target_positions(target_granges = gr_object,
                                         target_object = M_targetbs, target_cvg = cvg_targetbs,
                                         dmp_regions = celltype_specific_dmrs)$dmp_granges
      zmat_sub <- subsetByOverlaps(celltype_specific_dmrs, keep_dmrs,
                                   type = "equal")
      if(verbose){
        mes <- "[estimatecc] BSseq object contained
                %s out of %s celltype-specific regions."
        message(sprintf(mes, length(keep_dmrs), length(celltype_specific_dmrs)))
      }
      
      if(init_param_method == "known_regions"){
        y_d0 <- as.data.frame(
          .pick_target_positions(
            target_granges = gr_object,
            target_object = M_targetbs,
            target_cvg = cvg_targetbs,
            dmp_regions = grd0)$avg_dmr_dat)
        a0init <- mean(colMeans(y_d0))
        sig0init <- mean(apply(y_d0, 2, var))
        y_d1 <- as.data.frame(
          .pick_target_positions(
            target_granges = gr_object,
            target_object = M_targetbs,
            target_cvg = cvg_targetbs,
            dmp_regions = grd1)$avg_dmr_dat)
        a1init <- mean(colMeans(y_d1))
        sig1init <- mean(apply(y_d1, 2, var))
      }
    }
    
    ymat <- mcols(keep_dmrs); 
    if(!is.null(sampleNames(object))){ 
      colnames(ymat) <- sampleNames(object)
    }
    n <- ncol(mcols(keep_dmrs))
    ids <- colnames(mcols(zmat_sub))
    zmat <- as.matrix(as.data.frame(mcols(zmat_sub)))
    
  } else {
    
    if(verbose){
      if(is.null(region_mat)){
        stop("User must provide region_mat if object\n
                    is not an RGChannelSet or BSseq object.") }
    }
    ymat <- t(t(object))
    n <- ncol(object)
    ids <- colnames(region_mat)
    zmat <- as.matrix(region_mat)
    gr_object <- GRanges()
    keep_dmrs <- GRanges()
  }
  
  if(init_param_method == "known_regions"){
    list(ymat = ymat, n = n, ids = ids, zmat = zmat,
         a0init = a0init, a1init = a1init,
         sig0init = sig0init, sig1init = sig1init,
         gr_object = gr_object, keep_dmrs = keep_dmrs)
  } else {
    
    list(ymat = ymat, n = n, ids = ids, zmat = zmat,
         gr_object = gr_object, keep_dmrs = keep_dmrs)
  }
}


#' @title .initialize_theta
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
#' @importFrom stats runif
#'
.initialize_theta <- function(n, K, alpha0 = NULL, alpha1 = NULL,
                              sig0 = NULL, sig1 = NULL, tau = NULL){
  if(is.null(alpha0)){ alpha0 <- runif(1, 0.10, 0.3) }
  if(is.null(alpha1)){ alpha1 <- runif(1, 0.70, 0.9) }
  if(is.null(sig0)){ sig0 <- runif(1, 0.01, 0.05)^2 }
  if(is.null(sig1)){ sig1 <- runif(1, 0.01, 0.05)^2 }
  if(is.null(tau)){ tau <- runif(1, 0.001, 0.01)^2 }
  return(data.frame("alpha0" = alpha0, "alpha1" = alpha1,
                    "sig0" = sig0, "sig1" = sig1, "tau" = tau,
                    "lik_fun" = 0))
}

#' @title .initializeMLEs 
#' 
#' @description Helper functions to initialize 
#' MLEs in estimatecc. 
#' 
#' @param init_param_method Add here.
#' @param n Add here. 
#' @param K Add here.
#' @param Ys Add here. 
#' @param Zs Add here. 
#' @param a0init Add here. 
#' @param a1init Add here. 
#' @param sig0init Add here. 
#' @param sig1init Add here. 
#' @param tauinit Add here. 
#' 
.initializeMLEs <- function(init_param_method, n, K, Ys, Zs, a0init, a1init,
                            sig0init, sig1init, tauinit)
{
  
  init_theta <- .initialize_theta(n=n, K=K, alpha0=a0init,
                                  alpha1=a1init, sig0=sig0init,
                                  sig1=sig1init, tau=tauinit)
  
  # Initial guess for cell composition estimates
  x0 = (1-Zs) * init_theta$alpha0 + (Zs) * init_theta$alpha1
  init_pi_mle <- abs(t(apply(Ys, 2, function(ys){
    solve.QP(Dmat = (t(x0)%*%x0),
             dvec = t(x0)%*%t(t(ys)), 
             Amat = cbind(rep(1,K), diag(K)),
             bvec = c(1, rep(0, K)), 
             meq = 1)$sol })))
  
  return(list("init_pi_mle" = init_pi_mle, "init_theta" = init_theta))
}

#' @title Helper function to take the product of Z and
#' cell composition estimates
#'
#' @description Helper function which is the product of
#' Z and pi_mle
#'
#' @param Zs Cell type specific regions of dimension R x K
#' @param pi_mle cell composition MLE estimates
#' 
.WFun <- function(Zs, pi_mle){
  W0 = (1-Zs) %*% t(pi_mle)
  W1 = Zs %*% t(pi_mle)
  list("W0" = W0, "W1" = W1)
}

#' @title Expectation step
#'
#' @description Expectation step in EM algorithm for methylCC
#'
#' @param Ys observed methylation levels in samples provided by user
#' of dimension R x n
#' @param Zs Cell type specific regions of dimension R x K
#' @param current_pi_mle cell composition MLE estimates of dimension K x n
#' @param current_theta other parameter estimates in EM algorithm
#' @param meth_status Indicator function corresponding to regions that
#' are unmethylated (meth_status=0) or methylated (meth_status=1)
#'
.methylcc_estep <- function(Ys, Zs, current_pi_mle, 
                             current_theta,
                             meth_status = 0)
{
  R <- nrow(Zs)
  n <- ncol(Ys)
  W <- .WFun(Zs=Zs, pi_mle=current_pi_mle)
  
  # dim(mat11.inv) = (R*n x R*n) in block diagonal form
  mat11.inv <- lapply(seq_len(nrow(Zs)), function(r){
    # take inverse of individual blocks
    solve( t(t(W$W0[r,])) %*% t(W$W0[r,]) * current_theta$sig0 +
             t(t(W$W1[r,])) %*% t(W$W1[r,]) * current_theta$sig1 +
             diag(current_theta$tau, n) )    })
  
  if(meth_status == 0){
    # dim(mat12.0) = (R*n x R) in block diagonal form
    mat12.0 <- lapply(seq_len(nrow(Zs)), function(r){
      t(t(W$W0[r,])) * current_theta$sig0 })
    
    # dim(mat22.0) = (R x R)
    mat22.0 <- diag(current_theta$sig0, R)
    
    # dim(mom1) = dim(mom2) = (R x 1)
    fast.mat.mom1 <- t(t(sapply(seq_len(nrow(Zs)), function(r){
      t(mat12.0[[r]]) %*% mat11.inv[[r]] %*%
        (t(t(Ys[r,])) - current_theta$alpha0*t(t(W$W0[r,])) - current_theta$alpha1*t(t(W$W1[r,])) )
    })))
    fast.mat.mom2 <- diag(sapply(seq_len(nrow(Zs)), function(r){
      t(mat12.0[[r]]) %*% mat11.inv[[r]] %*% mat12.0[[r]]
    }), R)
    
    mom1 = current_theta$alpha0 * matrix(rep(1, R), ncol=1) + fast.mat.mom1
    mom2 = mom1^2 + t(t(diag(mat22.0 - fast.mat.mom2 )))
  }
  if(meth_status == 1){
    # dim(mat12.1) = (R*n x R) in block diagonal form
    mat12.1 <- lapply(seq_len(nrow(Zs)), function(r){
      t(t(W$W1[r,])) * current_theta$sig1 })
    
    # dim(mat22.1) = (R x R)
    mat22.1 <- diag(current_theta$sig1, R)
    
    # dim(mom1) = dim(mom2) = (R x 1)
    fast.mat.mom1 <- t(t(sapply(seq_len(nrow(Zs)), function(r){
      t(mat12.1[[r]]) %*% mat11.inv[[r]] %*%
        (t(t(Ys[r,])) - current_theta$alpha0*t(t(W$W0[r,])) - current_theta$alpha1*t(t(W$W1[r,])) )
    })))
    fast.mat.mom2 <- diag(sapply(seq_len(nrow(Zs)), function(r){
      t(mat12.1[[r]]) %*% mat11.inv[[r]] %*% mat12.1[[r]]
    }), R)
    
    mom1 = current_theta$alpha1 * matrix(rep(1, R), ncol=1) + fast.mat.mom1
    mom2 = mom1^2 + t(t(diag(mat22.1 - fast.mat.mom2 )))
  }
  list("mom1" = mom1, "mom2" = mom2)
}

#' @title Maximization step
#'
#' @description Maximization step in EM Algorithm for methylCC
#'
#' @param Ys observed methylation levels in samples provided by user
#' of dimension R x n
#' @param Zs Cell type specific regions of dimension R x K
#' @param current_pi_mle cell composition MLE estimates of dimension K x n
#' @param current_theta other parameter estimates in EM algorithm
#' @param estep0 Results from expectation step for unmethylated regions
#' @param estep1 Results from expectation step for methylated regions
#'
#' @importFrom quadprog solve.QP
#' 
.methylcc_mstep <- function(Ys, Zs, current_pi_mle, 
                             current_theta, 
                             estep0, estep1)
{
  R <- dim(Zs)[1]
  K <- dim(Zs)[2]
  n <- ncol(Ys)
  new_theta <- data.frame("alpha0" = mean(estep0$mom1),
                          "alpha1" = mean(estep1$mom1),
                          "sig0" = mean(estep0$mom2) - mean(estep0$mom1)^2,
                          "sig1" = mean(estep1$mom2) - mean(estep1$mom1)^2)
  x0 <- (sweep((1-Zs), 1, estep0$mom1, FUN = "*") +
           sweep(Zs,     1, estep1$mom1, FUN = "*"))
  new_pi_mle <- abs(t(apply(Ys, 2, function(ys){
    solve.QP(Dmat = (t(x0)%*%x0),
             dvec = t(x0)%*%t(t(ys)), Amat = cbind(rep(1,K), diag(K)),
             bvec = c(1, rep(0, K)), meq = 1)$sol })))
  tau_est <- (1/(R*n)) * sum((Ys - x0 %*% t(new_pi_mle))^2)
  new_theta <- cbind(new_theta, "tau" = tau_est,
                     "lik_fun" = -((R*n)/2)*log(2*pi*tau_est) -
                       (1/(2*tau_est)) * sum((Ys - x0 %*% t(new_pi_mle))^2) -
                       ((R*n)/2) * log(2*pi*new_theta$sig0) -
                       (n/(2*new_theta$sig0)) * sum((estep0$mom1 - new_theta$alpha0)^2) -
                       ((R*n)/2) * log(2*pi*new_theta$sig1) -
                       (n/(2*new_theta$sig1)) * sum((estep1$mom1 - new_theta$alpha1)^2))
  
  return(list("new_theta" = new_theta, "new_pi_mle" = new_pi_mle))
}


#' @title .methylcc_engine
#' 
#' @description Helper function for estimatecc
#' 
#' @param Ys observed methylation levels in samples provided by user
#' of dimension R x n
#' @param Zs Cell type specific regions of dimension R x K
#' @param current_pi_mle cell composition MLE estimates of dimension K x n
#' @param current_theta other parameter estimates in EM algorithm
#' @param epsilon Add here.
#' @param max_iter Add here.
#' 
.methylcc_engine <- function(Ys, Zs, current_pi_mle, 
                              current_theta, epsilon, max_iter)
{
  final_theta <- NULL
  final_pi_mle <- NULL
  count = 0
  
  repeat{
    # E-Step
    EStep0 <- 
      .methylcc_estep(Ys = Ys, Zs = Zs, 
                       current_pi_mle = current_pi_mle, 
                       current_theta = current_theta, 
                       meth_status = 0)
    EStep1 <- 
      .methylcc_estep(Ys = Ys, Zs = Zs, 
                       current_pi_mle = current_pi_mle, 
                       current_theta = current_theta, 
                       meth_status = 1)
    
    # M-Step
    update_step <- 
      .methylcc_mstep(Ys = Ys, Zs = Zs, 
                       current_pi_mle = current_pi_mle, 
                       current_theta = current_theta, 
                       estep0=EStep0,
                       estep1=EStep1)
    
    if(abs(update_step$new_theta$lik_fun - current_theta$lik_fun) >= epsilon){
      final_pi_mle <- update_step$new_pi_mle
      final_theta <- rbind(final_theta, update_step$new_theta)
      current_pi_mle <- update_step$new_pi_mle
      current_theta <- update_step$new_theta
      count = count + 1
      
      if(count > max_iter){
        break
      }
    } else {
      break
    }
  }
  return(list("pi_mle" = final_pi_mle, "theta" = final_theta))
}

