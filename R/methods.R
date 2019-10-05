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
#' @param object an object of class \code{estimatecc}.
#' 
#' @return Returns the cell composition estimates
#' 
#' @export
#' 
#' @examples
#' library(FlowSorted.Blood.450k)
#' data(FlowSorted.Blood.450k)
#' rgset <- FlowSorted.Blood.450k[,
#'      pData(FlowSorted.Blood.450k)$CellTypeLong %in% "Whole blood"]
#' est <- estimatecc(object = rgset) 
#' cell_counts(est)
#' 
setMethod(
  f = "cell_counts", 
  signature = "estimatecc",
  definition = function(object) {
    return(object@cell_counts)
  }
)
