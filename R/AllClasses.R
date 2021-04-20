#' Internal representation of a call to \code{apt-copynumber-cyto}
#' 
#' This class encapulate all the argumnets requested by the command-line
#' applitaction \code{apt-copynumber-cyto} of the suit Affymetrix Power
#' Tools.
#'
#' @exportClass APTCopyNumberParam
setClass( "APTCopyNumberParam",
    representation =
        representation( 
            tool          = "character",
            level         = "character",
            os            = "character",
            cel.list      = "character",
            cel.path      = "character",
            output.path   = "character",
            analysis.path = "character",
            cdf           = "character",
            chrX          = "character", 
            chrY          = "character",
            qca           = "character",
            qcc           = "character",
            snp           = "character",
            annot.db      = "character",
            refmodel      = "character",
            params        = "character",
            call          = "character"
        )
)

#' Internal representation of a call to \code{apt-geno-qc}
#' 
#' This class encapulate all the argumnets requested by the command-line
#' applitaction \code{apt-geno-qc} of the suit Affymetrix Power Tools.
#' 
#' @exportClass APTGenoQcParam
setClass( "APTGenoQcParam",
    representation =
        representation(
            tool          = "character",
            level         = "character",
            os            = "character",
            analysis.path = "character",
            xml.config    = "character",
            cel.list      = "character",
            output.path   = "character",
            call          = "character"
        )
)

#' Internal representation of a call to \code{apt-geno-qc}
#' 
#' This class encapulate all the argumnets requested by the command-line
#' applitaction \code{apt-geno-qc} of the suit Affymetrix Power Tools.
#' 
#' @exportClass APTProbeSetGenotypeParam
setClass( "APTProbeSetGenotypeParam",
    representation = 
        representation(
            tool          = "character",
            level         = "character",
            os            = "character",
            analysis.path = "character",
            xml.config    = "character",
            cel.list      = "character",
            output.path   = "character",
            params        = "character",
            call          = "character"
        )
)

#' Class result of the \code{CytoQCView}.
#' 
#' This class encapulate all information to draw one of the three 
#' visualizations ('snp', 'calling' and 'sc') of Affymetrix CytoScan.
#' 
#' @exportClass CytoQCViewer
setClass( "CytoQCViewer",
    representation = 
        representation(
            location      = "character",
            type          = "character",
            individual    = "character",
            verbose       = "logical"
        )
)