#' Function to create a SnpMatrix or a file compatible with PLINK from 
#' Affymetrix GenomeWide SNP (6.0/5.0) .CEL files.
#' 
#' This function uses the raw .CEL files from both GenomeWide SNP 5.0 and 6.0 
#' arrays to create a list containing a map file and a SnpMatrix with the 
#' genotype of each probe in the array. It can also be used to create a 
#' \code{.tped} file, compatible with \code{PLINK} (is in charge of the user 
#' to create the partner \code{.tfam} file).
#' 
#' @param cel.files Location where the \code{.CEL} files.
#' @param cel.platform By default \code{"affy6"}, it expects raw .CEL files 
#' from GenomeWide SNP 6.0. Can be set to \code{"affy5"} to use raw .CEL files 
#' from GenomeWide SNP 5.0.
#' @param cel.attime By default '4'. If the package is run over an OS 
#' compatible with the R package \code{parallel}, number of .CEL files that 
#' will be processed at the same time.
#' @param cel.validate By default \code{TRUE}. Can be set to \code{FALSE}
#' to avoid the step 'validation of .CEL files' from \code{clrmm}.
#' @param output.type By default \code{"snpmatrix"}. It can be set to 
#' \code{"plink"}. Use to determiner the type of object/file the function will 
#' generate.
#' @param output.name If \code{"output.type"} is set to \code{"plink"}, it 
#' must contain the name of the output \code{.tped} file (with no extension).
#' @param vervose By default \code{FALSE}. If TRUE the function will shows 
#' messages indicating the process.
#' @examples 
#' \donotrun{
#' # To create a SnpMatrix Container (a list with a map and a SnpMatrix)
#' smc <- Gw2SnpMatrix(
#'  cel.files="C:/Users/brge/Desktop/GW6",
#'  output.type="snpmatrix"
#' )
#' }
#' \donotrun{
#' # To create the .tped file compatible with PLINK
#' Gw2SnpMatrix(
#'  cel.files="C:/Users/brge/Desktop/GW6",
#'  output.type="plink",
#'  output.name="C:/Users/brge/Desktop/PRAD"
#' )
#' }
#' @export Gw2SnpMatrix
Gw2SnpMatrix <- function (cel.files, cel.platform = "affy6", cel.attime = 4, 
    cel.validate = TRUE, output.type = "snpmatrix", output.name = NA, 
    verbose = FALSE) {

    # SETUP
    OptionsAffy.Init("Gw2SnpMatrix", verbose)
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
    output.type <- value(output.type, c("snpmatrix", "plink"))
    if (is.na(output.type)) {
        stop("Argument 'output.type' must be 'snpmatrix' or 'plink'.")
    }
    if (output.type == "plink" & is.na(output.name)) {
        stop("If 'output.type' is 'plink', 'output.name' must be filled.")
    }
    # /FUN CHECKS

    # ARGS
    m <- paste0("\n\n",
        " - Input:      ", cel.files, "\n",
        " - Platform:   ", cel.platform, "\n",
        " - CEL x time: ", cel.attime, "\n",
        " - Output t.:  ", output.type)
    if(output.type == "plink") {
        m <- paste0(m, "\n - Output f.:  ", output.name)
    }
    ScreenAffy(paste0(m, "\n"))
    # /ARGS

    # Locating CEL files
    ScreenAffy("Locating CEL files.")
    cel.list <- list.celfiles(cel.files, full.names=TRUE)
    if (length(cel.list) == 0) { 
        stop("No CEL files in directory '", cel.files, "'.")
    }
    # /Locating CEL files

    # Validating CEL files
    if (cel.validate) {
        ScreenAffy("Validating CEL files.")
        tryCatch( { 
                validCEL(cel.list)
            }, warnning=function(w){ 
                stop("Bad (w) CEL files.\n", w)
            }, error=function(e){ 
                stop("Bad (e) CEL files.\n", e)
            }
        )
    }
    # /Validating CEL files

    # Creating the MAP content
    if (cel.platform == "affy6") {
        ScreenAffy("Getting annotation from 'pd.genomewidesnp.6'.")
        con <- pd.genomewidesnp.6@getdb();
    } else {
        ScreenAffy("Getting annotation from 'pd.genomewidesnp.5'.")
        con <- pd.genomewidesnp.5@getdb();
    }
    snp.t <-  dbGetQuery(con, "select * from featureSet")
    # /Creating the MAP content


    if (output.type == "snpmatrix") {
        ScreenAffy(
            "Processing individuals to generate a 'SnpMatrix Container'."
        )

        object = list()
        object$genotype <- GwProcess.SnpMatrix(cel.list)

        ScreenAffy("Generating MAP content.")
        object$map <- snp.t[ , c("dbsnp_rs_id", "chrom", "physical_pos", "allele_a",
            "allele_b", "strand")]
        rownames(object$map) <- snp.t[ , "man_fsetid"]
        names(object$map)[2:3] <- c("chromosome", "position")

        nmap <- cbind(object$map[ , c("chromosome", "dbsnp_rs_id")], rep(NA, dim(object$map)[1]),
            object$map[ , c("position", "allele_a", "allele_b")])
        colnames(nmap) <- c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")

        object$map <- nmap[colnames(object$genotype), ]
        object$map <- object$map[order(object$map$chromosome, object$map$position), ]

        return(object)
    } else {
        ScreenAffy(
            "Processing individuals to generate a 'TPED' file."
        )

        ScreenAffy("Selecting SNP data from annotation file.")
        common <- cbind( 
            snp.t[ , c("chrom", "dbsnp_rs_id")],
            rep(0, dim(snp.t)[1]),
            snp.t[ , c("physical_pos")])
        rownames(common) <- snp.t[ , "man_fsetid"]

        ScreenAffy("Creating the TPED content.")
        output.name <- paste0(output.name, ".tped")

        GwProcess.Plink(cel.list, common, snp.t$allele_a, 
            snp.t$allele_b, snp.t$strand, output.name)

        object <- NA
    }

    return(object)
}
