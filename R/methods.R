#' @title estimateCC Class
#'
#' @description This is the estimateCC class
#' @aliases estimateCC
#'
#' @exportClass estimateCC
#'
setClass(Class = "estimateCC",
         representation = representation(
           summary = "list",
           includeSE = "logical",
           cellcounts = "data.frame",
           mles = "data.frame",
           mleList = "list",
           YMat = "data.frame",
           ZMat = "data.frame",
           grObject = "GRanges",
           keepDMRs = "GRanges",
           countsMean = "array",
           countsSE = "array",
           mleSE = "array")
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
#' @param ... other
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
#' @param ... other
cellcounts.estimateCC <- function(object) object@cellcounts

#' @title cellcounts
#' @rdname cellcounts
#' @export
setMethod("cellcounts", signature(object="estimateCC"), cellcounts.estimateCC)

#' @title Accessors for the 'mles' slot of a estimateCC object.
#'
#' @description Accessors for the 'mles' slot of a estimateCC object.
#'
#' @usage
#' \S4method{mles}{estimateCC}(object)
#'
#' @docType methods
#' @name mles
#' @rdname mles
#' @aliases mles mles,estimateCellCountSeq-method
#' @param object a \code{estimateCC} object
#' @param ... other
mles.estimateCC <- function(object) object@mles

#' @title mles
#' @rdname mles
#' @export
setMethod("mles", signature(object="estimateCC"), mles.estimateCC)


#' @title Accessors for the 'mleList' slot of a estimateCC object.
#'
#' @description Accessors for the 'mleList' slot of a estimateCC object.
#'
#' @usage
#' \S4method{mleList}{estimateCC}(object)
#'
#' @docType methods
#' @name mleList
#' @rdname mleList
#' @aliases mleList mleList,estimateCellCountSeq-method
#' @param object a \code{estimateCC} object
#' @param ... other
mleList.estimateCC <- function(object) object@mleList

#' @title mleList
#' @rdname mleList
#' @export
setMethod("mleList", signature(object="estimateCC"), mleList.estimateCC)
