#' Function to create compatible files with MAD or PennCNV from Affymetrix 
#' GenomeWide SNP (5.0/6.0) .CEL files.
#' 
#' This function uses the raw .CEL files from both GenomeWide SNP 5.0 and 6.0 
#' arrays to create a file for each individual compatible with the R packages 
#' \code{MAD} and \code{R-GADA}. It can also be used to create a compatible 
#' file with \code{PennCNV}.
#'
#' @param cel.files Location for the raw .CEL files.
#' @param output.name Filename in case \code{output.type} is set to 
#' \code{"penncnv"}.
#' @param output.type By default \code{"mad"} can be set to \code{"penncnv"} to
#' generate a file compatible with \code{PennCNV}.
#' @param cel.platform By default \code{"affy6"}, it expects raw .CEL files 
#' from GenomeWide SNP 6.0. Can be set to \code{"affy5"} to use raw .CEL files 
#' from GenomeWide SNP 5.0.
#' @param cel.validate By default \code{TRUE}. Can be set to \code{FALSE}
#' to avoid the step 'validation of .CEL files' from \code{clrmm}.
#' @param cel.genome By default set to \code{"hg19"}, can be set to 
#' \code{"hg18"}. 
#' @param cel.attime By default '4'. If the package is run over an OS 
#' compatible with the R package \code{parallel}, number of .CEL files that 
#' will be processed at the same time.
#' @param markers.attime Number of markers processed at a time.
#' @param fill.24 Default \code{TRUE}. It \code{TRUE} a false chromosome
#' 24 will be added to complete the output files.
#' @param verbose If set to \code{TRUE}, usefull information is show.
#' @examples 
#' \donotrun{
#' # Generate MAD files
#' Gw2Mad(
#'  cel.files="gw6.raw", 
#'  output.name="gw6.mad", 
#'  output.type="mad"
#' )
#' }
#' \donotrun{
#' # Generate PennCNV compatible file
#' Gw2Mad(
#'  cel.files="gw6.raw", 
#'  output.name="gw6.mad", 
#'  output.type="penncnv"
#' )
#' }
#' @export Gw2Mad
Gw2Mad <- function(cel.files, output.name, output.type = "mad",
    cel.platform = "affy6", cel.validate = TRUE, cel.genome = "hg19",
    cel.attime = 4, markers.attime = 50e3, fill.24 = TRUE, verbose = FALSE) {

    # SETUP
    OptionsAffy.Init("Gw2Mad", verbose)
    OptionsAffy.Cores(cel.attime)
    # /SETUP

    # FUN CHECKS
    cel.platform <- value(cel.platform, c("affy6", "affy5"))
    if (is.na(cel.platform)) {
        stop("Argument 'cel.platform' must be 'affy6' or 'affy5' for Genomewide 6.0 or 5.0, respectively.")
    }
    if (cel.files == "") { # TODO: Must be checked on windows
        cel.files <- "./"
    }
    cel.genome <- value(cel.genome, c("hg18", "hg19"))
    if (is.na(cel.genome)) {
        stop("Argument 'cel.genome' must be 'hg18' or 'hg19'.")
    }
    output.type <- value(output.type, c("mad", "penncnv"))
    if (is.na(output.type)) {
        stop("Argument 'output.type' must be 'mad' or 'penncnv'.")
    }
    if (cel.files == "") { # TODO: Must be checked on windows
        cel.files <- "./"
    }
    if(is.na(output.type)) {
        stop("Argument 'output.type' must be 'mad' or 'penncnv'.")
    }
    if(is.na(output.name) | output.name == '') {
        stop("Argument 'output.name' must be filled.")
    }
    # /FUN CHEKS

    crlmm.lim <- 30;
    crlmm.oc <- c(2, 1e3)
    cdfname <- ifelse(cel.platform=="affy6",
        "genomewidesnp6","genomewidesnp5")
    output.temp <- RandomName(prefix="tmp.gw.mad.")

    # ARGS
    m <- paste0( "\n\n",
        " - Input:          ", cel.files, "\n",
        " - Output type:    ", output.type, "\n",
        " - Output name:    ", output.name, "\n",
        " - Temp path:      ", output.temp, "\n",
        " - Platform:       ", cel.platform, "\n",
        " - Genome:         ", cel.genome, "\n",
        " - CEL x time:     ", cel.attime, "\n",
        " - markers x time: ", markers.attime, "\n",
        " - CEL Validate:   ", cel.validate, "\n",
        " - fill.24:        ", fill.24, "\n" )
    ScreenAffy(m)
    # /ARGS

    ocProbesets(markers.attime)
    ocSamples(cel.attime)

    # Locating CEL files
    ScreenAffy("Locating CEL files.")
    cel.list <- list.celfiles(cel.files, full.names = TRUE)
    if (length(cel.list) == 0) {
        stop("No CEL files in directory '", cel.files, "'.")
    }
    # /Locating CEL files

    # Validating CEL files
    if(cel.validate) {
        ScreenAffy("Validating CEL files.")
        tryCatch( {
                validCEL( cel.list )
            }, warnning=function(w) {
                stop("Bad (w) CEL files.\n", w)
            }, error=function(e) {
                stop("Bad (e) CEL files.\n", e)
            }
        )
    }
    # /Validating CEL files

    unlink(output.temp, recursive = TRUE, force = TRUE)
    unlink(output.name, recursive = TRUE, force = TRUE)
    dir.create(output.temp, showWarnings=FALSE)

    ldPath(output.temp)
    ScreenAffy("Loading CEL files (tmp: '", output.temp, "').")
    cnSet <- tryCatch({
            crlmm::constructAffyCNSet(cel.list, batch=rep("1",
                length(cel.list)), cdfName=cdfname, genome=cel.genome)
        }, warnning=function(w) {
            stop("Warning during 'cnSet' construction.\n", w)
        }, error=function(e) {
            stop("Error during 'cnSet' construction.\n", e)
        }
    )

    ScreenAffy("Performing 'CNRMA'.")
    tryCatch({
            crlmm::cnrmaAffy(cnSet)
        }, warnning=function(w) {
            stop("Warning performing 'CNRMA'.\n", w)
        }, error=function(e) {
            stop("Error performing 'CNRMA'.\n", e )
        }
    )

    ScreenAffy("Performing 'SNPRMA'.")
    tryCatch({
            crlmm::snprmaAffy(cnSet)
        },  warnning=function(w) {
            stop("Warning performing 'SNPRMA'.\n", w)
        }, error=function(e) {
            stop("Error performing 'SNPRMA'.\n", e)
        }
    )

    ScreenAffy("Genotyping.")
    tryCatch({
            crlmm::genotypeAffy(cnSet, gender=NULL)
        }, warnning=function(w) {
            stop("Warning on genotyping.\n", w)
        }, error=function(e) {
            stop("Error on genotyping.\n", e)
        }
    )

    ScreenAffy("Performing 'CRLMM'.")
    tryCatch({
            crlmm::crlmmCopynumber(cnSet, fit.linearModel=FALSE)
        }, warnning=function(w) {
            stop("Warning performing 'CRLMM'.\n", w)
        }, error=function(e) {
            stop("Error performing 'CRLMM'.\n", e)
        }
    )

    ScreenAffy("Coercing to 'BafLrrSetList'.")
    oligoList <- tryCatch( {
            BafLrrSetList(cnSet)
        }, warnning=function(w) {
            stop("Warning on coercing to 'BafLrrSetList'.\n", w)
        }, error=function(e) {
            stop("Error on coercing to 'BafLrrSetList'.\n", e)
        }
    )

    data <- CreateBafLrrContainer(cnSet, oligoList)

    if(output.type == "mad") {
        dir.create(output.name, showWarnings = FALSE)
        GwProcess.Mad(data, output.name, fill.24, verbose)
    }
    else {
       GwProcess.PennCNV(data, output.name, verbose)
    }

    unlink(output.temp, recursive=TRUE, force=TRUE)
}