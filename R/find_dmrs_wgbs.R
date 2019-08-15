#' @title find_dmrs_wgbs
#' 
#' @description Add here. 
#'
#' @param verbose Add here. 
#' @param gr_target Add here. 
#' @param include_cpgs Add here. 
#' @param include_dmrs Add here. 
#' @param num_cpgs Add here. 
#' @param num_regions Add here. 
#' @param bumphunter_beta_cutoff Add here. 
#' @param dmr_up_cutoff Add here. 
#' @param dmr_down_cutoff Add here. 
#' @param dmr_pval_cutoff Add here. 
#' @param cpg_pval_cutoff Add here. 
#' @param cpg_up_dm_cutoff Add here. 
#' @param cpg_down_dm_cutoff Add here. 
#' @param pairwise_comparison Add here. 
#' @param mset_train_flow_sort Add here. 
#' 
#' @export 
#' 
#' @return Add here. 
#' 
#' @details 
#' Add here. 
#' 
#' @aliases find_dmrs_wgbs
#' 
#' @docType methods
#' @name find_dmrs_wgbs
#' 
#' @import minfi 
#' @import FlowSorted.Blood.450k
#' @importFrom plyranges arrange 
#' @importFrom S4Vectors queryHits
#' @importFrom utils head
#' @importFrom bumphunter clusterMaker
#' @importFrom stats model.matrix
#' @importFrom genefilter rowttests
#' 
#' @examples
#' found_regions <- find_dmrs_wgbs()
#' 
#' @rdname find_dmrs_wgbs
#' @export
find_dmrs_wgbs <- function(verbose=TRUE, gr_target=NULL,
                      num_regions=50, bumphunter_beta_cutoff = 0.2, 
                      dmr_up_cutoff = 0.5, dmr_down_cutoff = 0.4,
                      dmr_pval_cutoff = 1e-11, 
                      dmrseq_search = FALSE, workers = 1) {
  
  
  library(rhdf5) 
  library(HDF5Array)
  library(bsseq)
  
  dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"
  hdf5_bs_se_path <- file.path(dataPath, "files_bsseq_hdf5_se")
  bs <- loadHDF5SummarizedExperiment(hdf5_bs_se_path)
  
  IDs = c("Gran", "CD4T", "CD8T", "Bcell","Mono", "NK")

  # extract beta values, phenotypic information and GRanges objects
  cell <- factor(bs$cell_type, levels = IDs)
  cell_levels <- levels(cell)
  
  # define design matrix to search for DMRs
  xmat = model.matrix(~ cell - 1)
  colnames(xmat) = cell_levels
  pData(bs) <- cbind(pData(bs), xmat)
  
  all_poss = diag(length(cell_levels))
  all_poss <- (all_poss == TRUE)
  colnames(all_poss) <- cell_levels
  
  # find celltype specific regions using only overlap CpGs in target object
  if(!is.null(gr_target)){ 
    # which of the WGBS CpGs overlap with the target CpGs
    zz <- findOverlaps(granges(bs), gr_target) 
    bs <- bs[queryHits(zz), ]
    if(verbose){
      mes <- "[estimatecc] gr_target is not null. Using %s overlapping CpGs."
      message(sprintf(mes, nrow(bs)))
    }
  }
  
  regs <- vector("list", length(IDs))
  regions_all <- GRanges() 
  zmat <- c() # regions_all, will contain all celltype-specific DMRs
  
  for(ind in seq_len(nrow(all_poss))){
    
    if(!dmrseq_search){
      if(verbose){
        mes <- "[estimatecc] Loading %s cell type-specific regions."
        message(sprintf(mes, cell_levels[ind]))
      }
      regs[[ind]] <- readRDS(file = file.path("data",
                                              paste0("blueprint_blood_regions_dmrseq_", 
                                                     tolower(IDs[ind]),".RDS"))) 
    } else {
      library(BiocParallel)
      register(MulticoreParam(workers))
    
      if(verbose){
        mes <- "[estimatecc] Searching for %s cell type-specific regions."
        message(sprintf(mes, cell_levels[ind]))
      }
      regs[[ind]] <- dmrseq(bs=bs, testCovariate=IDs[ind], minNumRegion = 3, 
                            cutoff = bumphunter_beta_cutoff, verbose = TRUE)
    } 
  
    keep_dmrs = .pick_target_positions(target_granges = granges(bs),
                    target_object = getCoverage(bs, type = "M"), 
                    target_cvg = getCoverage(bs, type = "Cov"), 
                    dmp_regions = regs[[ind]])$dmp_granges
    
    # y_regions are the beta values collapsed (CpGs averaged) by regions from bumphunter
    y_regions <- as.matrix(mcols(keep_dmrs))
      
    tmp <- genefilter::rowttests(x = y_regions,fac = factor(pData(bs)[, IDs[ind]]))
    regs[[ind]]$ttest_pval <- tmp$p.value
    regs[[ind]]$ttest_dm   <- tmp$dm 
    regs[[ind]]$dmr_up_max_diff <- apply(abs(sweep(y_regions, 2, (pData(bs)[, IDs[ind]]), FUN = "-")), 1, max)
    regs[[ind]]$dmr_down_max_diff <- apply(abs(sweep(y_regions, 2, (1 - pData(bs)[, IDs[ind]]), FUN = "-")), 1, max)
     
    # Apply quality control 
    keep_ind_regions <- (regs[[ind]]$qval < 0.05) & # (regs[[ind]]$area > 1) # & 
                        (regs[[ind]]$ttest_pval < dmr_pval_cutoff) # ideally less than 1e-11
    regs[[ind]] <- regs[[ind]][keep_ind_regions, ]
    gr_regions_up <- regs[[ind]][(regs[[ind]]$stat > 0) & 
                                   (regs[[ind]]$dmr_up_max_diff < dmr_up_cutoff), ]
    gr_regions_down <- regs[[ind]][(regs[[ind]]$stat < 0) & 
                                    (regs[[ind]]$dmr_down_max_diff < dmr_down_cutoff), ]
      
    if(length(gr_regions_up) > 0){
        gr_regions_up <- gr_regions_up[, names(mcols(gr_regions_up)) %in% 
                                         c("L", "area", "stat", "ttest_dm", 
                                           "ttest_pval", "dmr_up_max_diff")]
        names(mcols(gr_regions_up))[
          names(mcols(gr_regions_up)) == "dmr_up_max_diff"] <- "dmr_max_diff"
        gr_regions_up <- gr_regions_up[order(-gr_regions_up$stat), ]
        gr_regions_up <- gr_regions_up[seq_len(min(length(gr_regions_up), num_regions)), ]
      } else {
        gr_regions_up <- GRanges()
      }

    if(length(gr_regions_down) > 0){
      gr_regions_down <- gr_regions_down[, names(mcols(gr_regions_down)) %in% 
                                       c("L", "area", "stat", "ttest_dm", 
                                         "ttest_pval", "dmr_down_max_diff")]
      names(mcols(gr_regions_down))[
        names(mcols(gr_regions_down)) == "dmr_down_max_diff"] <- "dmr_max_diff"
      gr_regions_down <- gr_regions_down[order(gr_regions_down$stat), ]
      gr_regions_down <- gr_regions_down[seq_len(min(length(gr_regions_down), num_regions)), ]
    } else {
      gr_regions_down <- GRanges()
    }
    
    mcols(gr_regions_up)$status <- rep("Up", length(gr_regions_up))
    mcols(gr_regions_down)$status <- rep("Down", length(gr_regions_down))
    bump_mat_all <- c(gr_regions_up, gr_regions_down)
    mcols(bump_mat_all)$cellType <- rep(cell_levels[ind], length(bump_mat_all)) 
    
    if(verbose){
        mes <- "[estimatecc] Found %s %s cell type-specific regions."
        message(sprintf(mes, length(bump_mat_all), cell_levels[ind]))
    }
    if(length(bump_mat_all) > 0){ 
      regions_all <- c(regions_all, bump_mat_all)
    }
    if(length(gr_regions_up) > 0){
      zmat <- rbind(zmat, t(replicate(min(length(gr_regions_up), num_regions), 
                                      as.numeric(all_poss[ind,]))))
    }
    if(length(gr_regions_down) > 0){
      zmat <- rbind(zmat, t(replicate(min(length(gr_regions_down), num_regions), 
                                      as.numeric(!all_poss[ind,]))))
    }
  }

  colnames(zmat) <- cell_levels
  
  keep_dmrs = .pick_target_positions(target_granges = granges(bs),
                                                target_object = getCoverage(bs, type = "M"), 
                                                target_cvg = getCoverage(bs, type = "Cov"), 
                                                dmp_regions = regions_all)$dmp_granges
  y_regions <- as.matrix(mcols(keep_dmrs))
  
  profiles <- sapply(.splitit(cell),function(ind){ rowMeans(y_regions[,ind]) } )
  
  removeMe <- duplicated(regions_all)
  list(regions_all = regions_all[!removeMe,], 
       zmat = zmat[!removeMe,], 
       y_regions = y_regions[!removeMe,], 
       profiles = profiles[!removeMe,], 
       cell = cell, cell_mat = all_poss, 
       cell_levels = cell_levels, pd = pData(bs))
} 
