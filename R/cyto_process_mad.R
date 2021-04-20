
# CytoPreprocess
# -----------------------------------------------------------------------------
# Created: 24/12/2013

CytoProcess.Mad <- function(cychp.file, output.name, annotation, output.temp,
    save, bt=0.25, tt=0.75) {

    name <- gsub(".cyhd.cychp.txt", "", basename(cychp.file))
    ScreenAffy("Calling individual '", name ,"'.")

    # Call the python parser's to clean the data
    CytoParser(cychp.file, output.temp)
    # /Call the python parser's to clean the data

    # Create files' names
    genotype.file <- file.path(output.temp,
        paste0(name, ".cyhd.cychp.genotype.f.txt"), fsep=.Platform$file.sep)
    copynumber.file <- file.path(output.temp,
        paste0(name, ".cyhd.cychp.cnds.f.txt"), fsep=.Platform$file.sep)
    allelepeaks.file <- file.path(output.temp,
        paste0(name, ".cyhd.cychp.apds.f.txt"), fsep=.Platform$file.sep )
    # Create files' names


    # Load data
    ScreenAffy("Loading calls ('", name ,"').")
    continue <- TRUE
    gnds <- tryCatch( {
          fread( genotype.file, header=TRUE, colClasses="character")#, as.is=TRUE )
        }, warnning=function(w) {
            OptionsAffy.AddError(error.name = genotype.file,
                "Invalid file provided in 'genotype.file' (",
                    genotype.file ,") some error from 'python'.")
            continue <- FALSE
        }, error=function(e) {
            OptionsAffy.AddError(error.name = genotype.file,
                "Invalid file provided in 'genotype.file' (",
                    genotype.file ,") some error from 'python'.")
            continue <- FALSE
        }
    )
    if (!continue) {
        ScreenAffy("Failed calling (gnds) for individual '", name ,"'.")
        return(NA)
    }

    ScreenAffy("Loading allele peaks ('", name ,"').")
    apds <- tryCatch( {
          fread(allelepeaks.file, sep="\t", header=TRUE, colClasses="character")#, as.is=TRUE)
        }, warnning=function(w) {
            OptionsAffy.AddError(error.name = allelepeaks.file,
                "Invalid file provided in 'allelepeaks.file' (",
                    allelepeaks.file ,") some error from 'python'.")
        }, error=function(e) {
            OptionsAffy.AddError(error.name = allelepeaks.file,
                "Invalid file provided in 'allelepeaks.file' (",
                    allelepeaks.file ,") some error from 'python'.")
        }
    )
    if (!continue) {
        ScreenAffy("Failed calling (apds) for individual '", name ,"'.")
        return(NA)
    }

    ScreenAffy("Loading copynumber ('", name ,"').")
    cnds <- tryCatch( {
            fread(copynumber.file, header=TRUE, colClasses="character")#, as.is=TRUE)
        }, warnning=function(w) {
            OptionsAffy.AddError(error.name = copynumber.file,
                "Invalid file provided in 'copynumber.file' (",
                    copynumber.file ,") some error from 'python'.")
            continue <- FALSE
        }, error=function(e) {
            OptionsAffy.AddError(error.name = copynumber.file,
                "Invalid file provided in 'copynumber.file' (",
                    copynumber.file ,") some error from 'python'.")
        }
    )
    if (!continue) {
        ScreenAffy("Failed calling (cnds) for individual '", name ,"'.")
        return(NA)
    }
    # /Load data


    # Process
    setkey(gnds, ProbeSetName)
    setkey(apds, ProbeSetName)
    setkey(cnds, ProbeSetName)


    ScreenAffy("Quantile filtering for individual '", name ,"'.")
    apds$AllelePeaks0 <- as.numeric(apds$AllelePeaks0)

    n <- quantile(as.numeric(apds$AllelePeaks0), c(bt, tt), rm.na=TRUE)
    apds$AllelePeaks0[apds$AllelePeaks0 < n[1]] <- n[1]
    apds$AllelePeaks0[apds$AllelePeaks0 > n[2]] <- n[2]

    ap <- (apds$AllelePeaks0 - n[1] ) / (n[2] - n[1])
    ap <- 1 - ap;
	  names(ap) <- apds$ProbeSetName
    ScreenAffy("Creating data.frame for individual '", name ,"'.")

    name <- file.path(output.name, name, fsep=.Platform$file.sep)

  	annotation$Log2Ratio <- NA
  	annotation$GType <- NA
  	annotation$B.Allele.Freq <- NA

  	annotation$Log2Ratio <- cnds[rownames(annotation), "Log2Ratio", with=FALSE]$Log2Ratio
  	annotation$GType <- gnds[rownames(annotation), "Call", with=FALSE]$Call
    annotation[names(ap), "B.Allele.Freq" ] <- ap

    # madfile <- data.frame(cnds[rownames(apds), "ProbeSetName"],
        # cnds[rownames(apds), "Chromosome"],
        # cnds[rownames(apds), "Position"],
        # cnds[rownames(apds), "Log2Ratio"],
        # gnds[rownames(apds), "Call"], ap)

    colnames(annotation) <- c("Name", "Chr", "Position", "Log.R.Ratio", "GType",
        "B.Allele.Freq")

    rownames(annotation) <- annotation$Name
    rm(apds, gnds, cnds)

    ScreenAffy("Converting 'NC' to 'NA' ('", name ,"').")
    annotation$GType[annotation$GType=="NC"] <- NA

    rownames(annotation) <- annotation$Name

    if (save) {
        ScreenAffy("Writing file for individual '", name ,"'.")
        write.table(annotation, file=name, quote=FALSE, row.names=FALSE,
            sep="\t", na="NA")

        #probes <- as.character(annotation$Name)
        rm(annotation)
    } else {
        return(annotation[, c("Log.R.Ratio", "GType", "B.Allele.Freq")])
    }
    # /Process
}




## CytoPostprocess.Mad
## -----------------------------------------------------------------------------
## Created: 26/12/2013

# CytoPostprocess.Mad <- function(file.name, probes.common) {

    # ScreenAffy("Post-process for individual '", file.name ,"'.")

    # madfile <- read.table(file.name, header=TRUE, sep="\t")
    # rownames(madfile) <- madfile$Name
    # madfile <- madfile[probes.common, ]

    # write.table(madfile, file=file.name, quote=FALSE, row.names=FALSE,
        # sep="\t", na="NA")

    # return(NA)
# }
