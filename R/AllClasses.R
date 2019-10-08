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
#' \dontrun{
#' library(FlowSorted.Blood.450k)
#' data(FlowSorted.Blood.450k)
#' rgset <- FlowSorted.Blood.450k[,
#'      pData(FlowSorted.Blood.450k)$CellTypeLong %in% "Whole blood"]
#' est <- estimatecc(object = rgset) 
#' }
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