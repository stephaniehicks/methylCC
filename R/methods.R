#' @title estimatecc
#'
#' @description This is the S4 class estimatecc container
#'
#' @exportClass estimatecc
#'
setClass(Class = "estimatecc",
         representation = representation(
           summary = "list",
           cell_counts = "data.frame",
           theta = "data.frame",
           ymat = "data.frame",
           zmat = "data.frame",
           gr_object = "GRanges",
           keep_dmrs = "GRanges")
)

setMethod("show", "estimatecc",
          function(object){
            cat("estimatecc: Estimate Cell Composition of Whole Blood Samples using DNA methylation\n")
            cat("   Input object class: ", object@summary$class, "\n")
            cat("   Reference cell types: ", object@summary$celltypes, "\n")
            cat("   Number of Whole Blood Samples: ", paste(object@summary$n_samples), "\n")
          }
)

#' @title Accessors for the 'summary' slot of a estimatecc object.
#'
#' @description Accessors for the 'summary' slot of a estimatecc object.
#'
#' @usage
#' \S4method{summary}{estimatecc}(object)
#'
#' @docType methods
#' @name summary
#' @rdname summary
#' @aliases summary summary,estimatecc-method
#' @param object a \code{estimatecc} object
summary.estimatecc <- function(object) object@summary

#' @title summary
#' @rdname summary
#' @export
setMethod("summary", signature(object="estimatecc"), summary.estimatecc)


#' @title Accessors for the 'cell_counts' slot of a estimatecc object.
#'
#' @description Accessors for the 'cell_counts' slot of a estimatecc object.
#'
#' @usage
#' \S4method{cell_counts}{estimatecc}(object)
#'
#' @docType methods
#' @name cell_counts
#' @rdname cell_counts
#' @aliases cell_counts cell_counts,estimatecc-method
#' @param object a \code{estimatecc} object
cell_counts.estimatecc <- function(object) object@cell_counts

#' @title cell_counts
#' @rdname cell_counts
#' @export
setMethod("cell_counts", signature(object="estimatecc"), cell_counts.estimatecc)

