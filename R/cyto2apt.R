#' Function to preprocess the .CEL files from Affymetrix CytoScan 750K/HD
#' 
#' This function run as a command line tool the \code{apt-copynumber-cyto}
#' included into the Affymetrix Port Tools to process the raw .CEL files and 
#' generate the \code{.cychp.txt} files that are requested by \code{Cyto2Mad}
#' and \code{Cyto2SnpMatrix}.
#'
#' @param arguments Object of class \code{APTCopyNumberParam} with the 
#' 'standard' or 'custom' arguments set.
#' @param verbose By default FALSE. If TRUE the function will shows messages 
#' indicating the process.
#' @examples 
#' \dontrun{
#' aptParam <- APTparam(
#'  type="cytoscan", 
#'  level="standard",
#'  cel.list="C:/Users/brge/Desktop/TEA_1de2", 
#'  output.path="C:/Users/brge/Desktop/APT_OUT", 
#'  analysis.path="C:/Users/brge/Desktop/lib_cytoHD",
#'  cdf="CytoScanHD_Array.cdf", 
#'  chrX="CytoScanHD_Array.chrXprobes", 
#'  chrY="CytoScanHD_Array.chrYprobes", 
#'  qca="CytoScanHD_Array.r1.qca", 
#'  qcc="CytoScanHD_Array.r1.qcc", 
#'  snp="CytoScanHD_Array.snplist.txt", 
#'  annot.db="CytoScanHD_Array.na32.3.annot.db", 
#'  refmodel="CytoScanHD_Array.na32.3.v1.REF_MODEL"
#'  )
#'  Cyto2APT(aptParam)
#' }
#' @export Cyto2APT
Cyto2APT <- function(arguments, verbose=FALSE) {

    # SETUPT
    OptionsAffy.Init("Cyto2APT", verbose)
    # /SETUP

    # SYSTEM
    CheckSystem()
    # /SYSTEM

    ScreenAffy("Creating input file for APT.")
    cel.df <- data.frame(c("cel_files", 
        list.celfiles(arguments@cel.list, full.names=TRUE)))
    
    cel.f <- RandomName(prefix="tmp.cyto.apt.")
    arguments@cel.path <- paste0(cel.f, ".1")
    
    write.table(cel.df, file=arguments@cel.path, quote=FALSE, 
      row.names=FALSE, col.names=FALSE)

    

    ScreenAffy("Creating call to APT.")
    # call <- NULL
    # if(Sys.info()["sysname"] == "Linux") {
    #     call <- CytoCopynumber.Linux(cel.f, output.path, analysis.path,
    #         cdf.file, chrX.probes, chrY.probes, analysis.qca, 
    #         analysis.qcc, analysis.snp, analysis.annotdb, analysis.refmodel)
    # } else {
    #     if(Sys.info()[ "sysname" ] == "Windows") {
    #         call <- CytoCopynumber.Windows(cel.f, output.path, analysis.path,
    #             cdf.file, chrX.probes, chrY.probes, analysis.qca, 
    #             analysis.qcc, analysis.snp, analysis.annotdb, analysis.refmodel)
    #     } else {
    #         call <- CytoCopynumber.Mac(cel.f, output.path, analysis.path,
    #         cdf.file, chrX.probes, chrY.probes, analysis.qca, 
    #         analysis.qcc, analysis.snp, analysis.annotdb, analysis.refmodel)
    #     }
    # }


    # if (!is.null(call)) {
    #     ScreenAffy("Calling APT")
    #     system(call)
    # }

    run(arguments)

    unlink(cel.f, recursive=TRUE, force=TRUE)
}


