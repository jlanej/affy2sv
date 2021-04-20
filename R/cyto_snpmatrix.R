#' Function to create a \code{SnpMatrix} or a file compatible with \code{PLINK}
#' 
#' This function uses the files created with \code{Cyto2APT} to create a list 
#' containing a map file and a SnpMatrix with the genotype of each probe in 
#' the array. It can also be used to create a \code{.tped} file, compatible 
#' with \code{PLINK} (is in charge of the user to create the partner 
#' \code{.tfam} file).
#' 
#' @param cychp.files Location where the .cychp.txt files (obtained from Cyto2APT) are placed.
#' @param annotation.file NetAffx Annotation database file in CSV format.
#' @param cychp.attime By default 4, number of \code{.cychp} fieles processed at a time.
#' @param output.type By default 'snpmatrix'. It can be set to 'plink'. Use to determiner the type of object/file the function will generate.
#' @param output.name If 'output.type' is set to 'plink', it must contain the name of the output .tped file (with no extension).
#' @param verbose By default FALSE. If TRUE the function will shows messages indicating the process.
#' @examples 
#' \dontrun{
#' # To create a SnpMatrix Container (a list with a map and a SnpMatrix)
#'  smc <- Cyto2SnpMatrix(
#'  cychp.files="C:/Users/brge/Desktop/APT_OUT", 
#'  annotation.file="C:/Users/brge/Desktop/lib_cytoHD/CytoScanHD_Array.na32.3.annot.csv", 
#'  output.type="snpmatrix"
#'  )
#' }
#' \dontrun{
#' # To create the .tped file compatible with PLINK
#'  Cyto2SnpMatrix(
#'  cychp.files="C:/Users/brge/Desktop/APT_OUT", 
#'  annotation.file="C:/Users/brge/Desktop/lib_cytoHD/CytoScanHD_Array.na32.3.annot.csv", 
#'  output.type="plink",
#'  output.name="ASD"
#'  )
#' }
#' @export Cyto2SnpMatrix
Cyto2SnpMatrix <- function(cychp.files, annotation.file, cychp.attime = 4,
    output.type = "snpmatrix", output.name = NA, verbose = FALSE) {

    # SETUP
    OptionsAffy.Init("Cyto2SnpMatrix", verbose)
    OptionsAffy.Cores(cychp.attime)
    # /SETUP

    # SYSTEM
    CheckSystem(bits=FALSE)
    # /SYSTEM

    # FUN CHECKS
    if( cychp.files == "" ) {
        cychp.files = "./";
    }
    output.type <- value(output.type, c("snpmatrix", "plink"))
    if (is.na(output.type)) {
        stop("Argument 'output.type' must be 'snpmatrix' or 'plink'.")
    }
    if (output.type == "plink" & is.na(output.name)) {
        stop("If 'output.type' is 'plink', 'output.name' must be filled.")
    }
    files.complete <- list.files(cychp.files, pattern=".cyhd.cychp.txt",
        full.names=TRUE)
    if (length(files.complete) == 0) {
        stop("No 'cychp' files in folder '", cychp.files, "'.")
    }
    # /FUN CHECKS

    output.temp <- RandomName(prefix="tmp.cyto.snpmatrix.")

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
    dir.create(output.temp, showWarnings=FALSE)
    # /CLEAN

	  # ANNOTATION
    ScreenAffy("Loading annotation file.")
    annot <- tryCatch({
            read.csv(annotation.file, comment.char="#", header=TRUE)
        }, warnning=function(w) {
            stop("Invalid file provided in \"annotation.file\".\n",
            w, "\n")
        }, error=function(e) {
            stop("Invalid file provided in \"annotation.file\".\n",
            e, "\n")
        })
     annot <- annot[ , c("Probe.Set.ID", "Chromosome", "Physical.Position", "dbSNP.RS.ID", "Allele.A", "Allele.B", "Strand", "ChrX.pseudo.autosomal.region.1")]
    # /ANNOTATION


    # CALLING
    ScreenAffy("Calling individuals.")
    calls.ind <<- mclapply(X=files.complete, FUN=CytoProcess.SnpMatrix,
        output.temp = output.temp,
        ann.p=as.character(annot$Probe.Set.ID),
        mc.cores=cychp.attime
    )
    OptionsAffy.CheckErrors()
    # /CALLING



    # FILTERING CALLS, GET COMMON
    # ScreenAffy("Filtering calls.")
    # common <- tryCatch({
            # Reduce(intersect,
                # Map(function(x) { x$names }, calls.ind))
        # }, warnning=function(w) {
            # stop(w, "\n", "class(calls.ind): ", class(calls.ind), "\n",
                # "class(calls.ind[1]): ", class(calls.ind[1]))
        # }, error=function(e) {
            # stop(e, "\n", "class(calls.ind): ", class(calls.ind), "\n",
                # "class(calls.ind[1]): ", class(calls.ind[1]))
        # })
    # calls <- Map(function(x){ x$genotype[common] }, calls.ind)
    # /FILTERING CALLS, GET COMMON

    # FLATTING CALLS TO SINGLE OBJECT
    # ScreenAffy("Flatting data.")
    # calls.ind <- data.frame(calls[[1]])
    # if (length(calls) > 1) {
        # for(ii in 2:length(calls)) {
            # calls.ind <- cbind(calls.ind, calls[[ii]])
        # }
    # }
    # rm("calls")
    # FLATTING CALLS TO SINGLE OBJECT

    ScreenAffy("Calling complet sampleset.")
  	calls <- do.call(cbind, calls.ind)
  	# rm(calls.ind)
  	colnames(calls) <- gsub(".cyhd.cychp.txt", "", basename(files.complete))
  	rownames(calls) <- as.character(annot$Probe.Set.ID)

    if (output.type == "snpmatrix") {
        ScreenAffy(
            "Processing individuals to generate a 'SnpMatrix Container'."
        )
    		map <- annot[ , c("dbSNP.RS.ID", "Chromosome", "Physical.Position", "Allele.A", "Allele.B")]
    		map$cM <- NA
    		colnames(map) <- c("snp.name", "chromosome", "position", "allele.1", "allele.2", "cM")
    		rownames(map) <- annot$Probe.Set.ID

        # CREATING MAP
        # map <- data.frame(chromosome=annot$Chromosome,
            # snp.name=annot$dbSNP.RS.ID,
            # cM=rep(NA, length(annot$dbSNP.RS.ID)),
            # position=as.numeric(annot$Physical.Position),
            # allele.1=annot$Allele.A,
            # allele.2=annot$Allele.B )
        # rownames(map) <- rownames(annot)
        # rm("annot")
        # /CREATING MAP

        # NAMING CALLS
        # colnames(calls.ind) <- gsub(".cyhd.cychp.txt", "",
            # basename(files.complete))
        # rownames(calls.ind) <- common
        # rm("common")
        # /NAMING CALLS

        # COERCION
        ScreenAffy("Coercing to a SnpMatrix.")
        calls.sm <- new("SnpMatrix", t(as(calls, "matrix")))
        # rm("calls.ind")
        # /COERCION

        ScreenAffy("Creating container.")
        object <- list()
        object$map <- map
        object$genotype <- calls.sm
        object$map <- object$map[colnames(object$genotype), ]
        object$map <- object$map[order(object$map$chromosome, object$map$position), ]

    } else {
        ScreenAffy(
            "Processing individuals to generate a 'TPED' file."
        )

        strand <- as.character(annot$Strand)
        allele.A <- as.character(annot$Allele.A)
        allele.B <- as.character(annot$Allele.B)

        ScreenAffy("Getting A and B alleles.")
        allComplA <- c("A", "C", "G", "T")[
            match(allele.A, c( "T", "G", "C", "A"))
        ]
        forwA <- ifelse(strand=="+", allele.A, allComplA)

        allComplB <- c("A", "C", "G", "T")[
            match( allele.B, c("T", "G", "C", "A"))
        ]
        forwB <- ifelse(strand == "+", allele.B, allComplB)
        rm("allele.A", "allele.B", "strand")

        ScreenAffy("Creating CALLS dataset." )
        alleleA <- alleleB <- matrix("", nrow=nrow(calls),
            ncol=ncol(calls), dimnames=dimnames(calls))
        for( ii in 1:ncol(calls)) {
            alleleA[ , ii] <- ifelse(as.numeric(calls[ , ii]) < 3,
                forwA, forwB)
            alleleB[ , ii] <- ifelse(as.numeric(calls[ , ii]) < 2,
                forwA, forwB)
        }

        ScreenAffy("Coercing TPED to matrix.")
        p1 <- seq(1, ncol(calls)*2, 2)
        p2 <- seq(2, ncol(calls)*2, 2)

        tped <- matrix("", ncol=ncol(calls) * 2,
            nrow=nrow(calls))

        for(ii in 1:ncol(calls)) {
            tped[ , p1[ii]] <- alleleA[ , ii]
            tped[ , p2[ii]] <- alleleB[ , ii]
        }

        rownames(tped) <- annot$Probe.Set.ID
        rm("alleleA", "alleleB", "forwA", "forwB")

        ScreenAffy("Selecting SNP data from annotation file.")
        #annot <- annot[rownames(tped), ]
        tped.an <- data.frame(annot$Chromosome, annot$dbSNP.RS.ID,
            as.numeric(annot$ChrX.pseudo.autosomal.region.1),
            as.numeric(annot$Physical.Position)
        )

        tped <- cbind(tped.an, tped)
        ScreenAffy("Writing TPED file to disk '", output.name, ".tped'.")
        write.table(tped, file=paste0(output.name, ".tped"), quote=FALSE,
            row.names=FALSE, col.names=FALSE)

        rm("tped")
        object <- NA
    }

    unlink(output.temp, recursive=TRUE, force=TRUE)
    return(object)
}
