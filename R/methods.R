#' @title estimateCC
#'
#' @description This is the S4 class estimateCC container
#'
#' @exportClass estimateCC
#'
setClass(Class = "estimateCC",
         representation = representation(
           summary = "list",
           cellcounts = "data.frame",
           theta = "data.frame",
           YMat = "data.frame",
           ZMat = "data.frame",
           grObject = "GRanges",
           keepDMRs = "GRanges")
)

setMethod("show", "estimateCC",
          function(object){
            cat("estimateCC: Estimate Cell Composition of Whole Blood Samples using DNA methylation\n")
            cat("   Input object class: ", object@summary$class, "\n")
            cat("   Reference cell types: ", object@summary$cellTypes, "\n")
            cat("   Number of Whole Blood Samples: ", paste(object@summary$nSamples), "\n")
          }
)

#' @title Accessors for the 'summary' slot of a estimateCC object.
#'
#' @description Accessors for the 'summary' slot of a estimateCC object.
#'
#' @usage
#' \S4method{summary}{estimateCC}(object)
#'
#' @docType methods
#' @name summary
#' @rdname summary
#' @aliases summary summary,estimateCC-method
#' @param object a \code{estimateCC} object
summary.estimateCC <- function(object) object@summary

#' @title summary
#' @rdname summary
#' @export
setMethod("summary", signature(object="estimateCC"), summary.estimateCC)


#' @title Accessors for the 'cellcounts' slot of a estimateCC object.
#'
#' @description Accessors for the 'cellcounts' slot of a estimateCC object.
#'
#' @usage
#' \S4method{cellcounts}{estimateCC}(object)
#'
#' @docType methods
#' @name cellcounts
#' @rdname cellcounts
#' @aliases cellcounts cellcounts,estimateCC-method
#' @param object a \code{estimateCC} object
cellcounts.estimateCC <- function(object) object@cellcounts

#' @title cellcounts
#' @rdname cellcounts
#' @export
setMethod("cellcounts", signature(object="estimateCC"), cellcounts.estimateCC)

