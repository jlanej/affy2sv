#' Function to create compatible txt files for \code{MAD} or for \code{PennCNV}
#' 
#' This function uses the files created with \code{Cyto2APT} to create a file 
#' for each individual compatible with the R packages \code{MAD} and 
#' \code{R-GADA}. It can also be used to create a compatible file with 
#' \code{PennCNV}.
#' 
#' @param cychp.files Location where the .cychp.txt files (obtained from Cyto2APT) are placed.
#' @param output.name If 'output.type' is set to 'mad', location where the files compatible with MAD/R-GADa will be saved. If 'output.type' is set to 'penncnv', name of the generated file.
#' @param annotation.file NetAffx Annotation database file in CSV format.
#' @param output.type By default 'mad'. It can be set to 'penncnv'. Use to determiner the type of file(s) the function will generate.
#' @param cychp.attime By default '4'. If the package is run over an OS compatible with the R package 'parallel', number of cychp.txt files that will be processed at the same time.
#' @param bottom.quantile By default '0.25'. Bottom threshold on BAF quantile normalization.
#' @param top.quantile By default '0.75'. Top threshold on BAF quantile normalization.
#' @param verbose By default FALSE. If TRUE the function will shows messages indicating the process.
#' @examples 
#' \dontrun{
#' # To create files compatible with MAD/R-GADA
#'  Cyto2Mad(
#'  cychp.files="C:/Users/brge/Desktop/APT_OUT",
#'  output.name="C:/Users/brge/Desktop/MAD",
#'  annotation.file="C:/Users/brge/Desktop/lib_cytoHD/CytoScanHD_Array.na32.3.annot.csv", 
#'  output.type="mad",
#'  )
#' }
#' \dontrun{
#' # To create files compatible with PennCNV
#'  Cyto2Mad(
#'  cychp.files="C:/Users/brge/Desktop/APT_OUT",
#'  output.name="C:/Users/brge/Desktop/ASD",
#'  annotation.file="C:/Users/brge/Desktop/lib_cytoHD/CytoScanHD_Array.na32.3.annot.csv", 
#'  output.type="penncnv",
#'  )
#' }
#' @export Cyto2Mad
Cyto2Mad <- function(cychp.files, output.name, annotation.file, output.type = "mad",
    cychp.attime = 4, bottom.quantile=0.25, top.quantile=0.75, verbose = FALSE) {

    # SETUP
    OptionsAffy.Init("Cyto2Mad", verbose)
    OptionsAffy.Cores(cychp.attime)
    # /SETUP

    # SYSTEM
    CheckSystem(bits=FALSE)
    # /SYSTEM

    # FUN CHECKS
    if(cychp.files == "") {
		stop("Argument 'cychp.files' must be filled.")
    }
  	if(missing(annotation.file)) {
  		stop("Argument 'annotation.file' must be filled.")
  	}
    output.type <- value(output.type, c("mad", "penncnv"))
    if (is.na(output.type)) {
        stop("Argument 'output.type' must be 'mad' or 'penncnv'.")
    }
    if (is.na(output.name) | output.name == '') {
        stop("Argument 'output.name' must be filled.")
    }
    files.complete <- list.files(cychp.files, pattern=".cyhd.cychp.txt",
        full.names=TRUE)
    if (length(files.complete) == 0) {
        stop("No 'cychp' files in folder '", cychp.files, "'.")
    }
    # /FUN CHECKS

    output.temp <- RandomName(prefix="tmp.cyto.mad.")

    # ARGS
    m <- paste0( "\n\n",
        " - Input:        ", cychp.files, "\n",
		    " - Annotation:   ", annotation.file, "\n",
        " - Temp path:    ", output.temp, "\n",
        " - cychp x time: ", cychp.attime, "\n",
        " - Output.t:     ", output.type, "\n",
        " - Output.f:     ", output.name, "\n");
    ScreenAffy(m)
    # /ARGS

    # CLEAN
    unlink(output.name, recursive=TRUE, force=TRUE)
    unlink(output.temp, recursive=TRUE, force=TRUE)
    dir.create(output.name, showWarnings=FALSE)
    dir.create(output.temp, showWarnings=FALSE)
    # /CLEAN

	# ANNOTATION
    ScreenAffy("Loading annotation file.")
    ann <- tryCatch({
            read.csv(annotation.file, header=TRUE, comment.char="#")
        }, warnning=function(w) {
            stop("Invalid file provided in \"annotation.file\".\n",
            w, "\n")
        }, error=function(e) {
            stop("Invalid file provided in \"annotation.file\".\n",
            e, "\n")
        })
	ann <- ann[ , c("Probe.Set.ID", "Chromosome", "Physical.Position")]
  rownames(ann) <- ann$Probe.Set.ID
    # /ANNOTATION

    data <- plapply(X=files.complete, FUN=CytoProcess.Mad,
  		annotation=ann,
  		output.name=output.name,
      output.temp=output.temp,
  		bt=bottom.quantile,
  		tt=top.quantile,
      save=(output.type == "mad")
	)

    ScreenAffy("Number of individuals: ", length(data))
	if (output.type == "penncnv") {
        unlink(output.name, recursive=TRUE, force=TRUE)

		y <- do.call(cbind, data)
		result <- cbind(ann, data)

		nn <- sapply(basename(files.complete), function(x) strsplit(x, "\\.")[[1]][1])
		nn <-  unlist(lapply(nn, function(x) paste(x, c(".Log R Ratio", ".GType", ".B Allele Freq"), sep="")))
		colnames(result) <- c("Name", "Chr", "Position", nn)

        write.table(result, file=paste0(output.name, ".penncnv"),
            quote=FALSE, row.names=FALSE, sep="\t", na="NA");
    }

    unlink(output.temp, recursive=TRUE, force=TRUE)
    ScreenAffy("Remove temp folder '", output.temp, "'.")
}
