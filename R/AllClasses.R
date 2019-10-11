#' @title the estimatecc class
#'
#' @description Objects of this class store all 
#' the values needed information to work with a
#' estimatecc object
#' 
#' @slot summary information about the samples and 
#' regions used to estimate cell composition
#' @slot cell_counts cell composition estimates
#' 
#' @return \code{summary} returns the summary information 
#' about the cell composition estimate procedure and 
#' \code{cell_counts} returns the cell composition estimates
#' 
#' @name estimatecc-class
#' @import methods
#' @exportClass estimatecc
#' @aliases estimatecc-class
#'
#' @examples
#' # This is a reduced version of the FlowSorted.Blood.450k 
#' # dataset available by using BiocManager::install("FlowSorted.Blood.450k),
#' # but for purposes of the example, we use the smaller version 
#' # and we set \code{demo=TRUE}. For any case outside of this example for 
#' # the package, you should set \code{demo=FALSE} (the default). 
#' 
#' dir <- system.file("data", package="methylCC")
#' files <- file.path(dir, "FlowSorted.Blood.450k.sub.RData") 
#' if(file.exists(files)){
#'     load(file = files)
#' 
#'     set.seed(12345)
#'     est <- estimatecc(object = FlowSorted.Blood.450k.sub, demo = TRUE) 
#'     cell_counts(est)
#'  }   
#' 
setClass(
    Class = "estimatecc", 
    slot = list(
        summary = "list",
        cell_counts = "data.frame")
)

setMethod("show", "estimatecc",
          function(object){
            cat("estimatecc: Estimate Cell Composition of Whole Blood 
                Samples using DNA methylation\n")
            cat("   Input object class: ", object@summary$class, "\n")
            cat("   Reference cell types: ", object@summary$celltypes, "\n")
            cat("   Number of Whole Blood Samples: ", 
                paste(object@summary$n_samples), "\n")
            cat("   Name of Whole Blood Samples: ", 
                paste(object@summary$sample_names), "\n")
          }
)