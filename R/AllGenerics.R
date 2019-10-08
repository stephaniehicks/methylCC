#' Generic function that returns the cell composition 
#' estimates
#'
#' Given a estimatecc object, this function returns the 
#' cell composition estimates 
#' @rdname cell_counts
setGeneric(
    name = "cell_counts", 
    def = function(object) { 
        standardGeneric("cell_counts") 
    }
)
