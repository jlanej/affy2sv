#' Create an object for CytoScan visualization
#' 
#' This function creates an object of class \code{CytoQCViewer} that allows
#' to perform a visual quality control. Three plots are avaialble: snp, 
#' calling and sc.
#' 
#' @param path Location of the .cychp files.
#' @param visualization "snp", "int" or "sc".
#' @param individual name of the invividual for calling and sc visualizations.
#' @param verbose If set to \code{TRUE} usefull information is shown.
#'
#' @export CytoQCView
CytoQCView <- function(path, visualization, individual, verbose = FALSE) {
    visualization <- value(visualization, c("snp", "int", "sc"))
    if (is.na(visualization)) {
        stop("Wrong type of 'visualization'. Available are: 'snp', 'calling' and 'sc'.")
    }

    new("CytoQCViewer", location = path, type = visualization,
        individual = individual, verbose = verbose)
}