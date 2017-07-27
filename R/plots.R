#' @title Plot results from \code{estimateCC} function.
#'
#' @description This function plots the parameter estimates
#' and standard error estimates from the \code{estimateCC}
#' function.
#'
#' @param objectEst an estimateCC object from \code{estimateCC}
#' @param savePlot a TRUE/FALSE object argument
#' determining if the plot will be saved for further use or
#' immediately displayed on the screen.
#'
#' @return The cell composition estimates for each sample.
#' If \code{includeSE} argument was set to TRUE in
#' \code{estimateCC()}, then the standard errors will
#' also be plotted.
#'
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom dplyr inner_join
#'
#' @export
estimateCCPlot <- function(objectEst, savePlot = FALSE){

    counts <- cellcounts(objectEst)
    sampIDs <- objectEst@summary$sampleNames

    df = gather(cbind("samples" = sampIDs, counts), celltype, est, -samples)

    gmat <- ggplot(df, aes(y=est, x = celltype)) + facet_wrap(~samples) +
        geom_point() + xlab("Celltypes") + ylab("Cell composition estimates")

    if(objectEst@includeSE){
        se = gather(cbind("samples" = sampIDs, as.data.frame(objectEst@countsSE)),
                celltype, se, -samples)
        df <- inner_join(df, se, by = c("samples", "celltype"))
        limits <- aes(ymax = est + se, ymin=est - se)

        gmat <- gmat + geom_errorbar(limits, data = df, width=0.2)
    }

    if (savePlot) {
        suppressMessages(gmat)
    } else {
        suppressMessages(print(gmat))
    }
}
