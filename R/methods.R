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
setMethod(
  f = "cell_counts", 
  signature = "estimatecc",
  definition = function(object) {
    return(object@cell_counts)
  }
)
