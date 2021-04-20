
# CY2APT
# -----------------------------------------------------------------------------
# 14/12/2013

Cyto2APT.old <- function(cel.files, output.path, analysis.path, 
    cdf.file, chrX.probes, chrY.probes, analysis.qca,
    analysis.qcc, analysis.snp, analysis.annotdb, analysis.refmodel, 
    verbose = FALSE ) {

    # SETUPT
    OptionsAffy.Init("Cyto2APT", verbose)
    # /SETUP

    # SYSTEM
    CheckSystem()
    # /SYSTEM

    ScreenAffy("Creating input file for APT.")
    cel.df <- data.frame(c("cel_files", 
        list.celfiles(cel.files, full.names=TRUE)))
    
    cel.f <- RandomName(prefix="tmp.cyto.apt.")
    write.table(cel.df, file=paste0(cel.f, ".1"), quote=FALSE, 
      row.names=FALSE, col.names=FALSE)

    ScreenAffy("Creating call to APT.")
    call <- NULL
    if(Sys.info()["sysname"] == "Linux") {
        call <- CytoCopynumber.Linux(cel.f, output.path, analysis.path,
            cdf.file, chrX.probes, chrY.probes, analysis.qca, 
            analysis.qcc, analysis.snp, analysis.annotdb, analysis.refmodel)
    } else {
        if(Sys.info()[ "sysname" ] == "Windows") {
            call <- CytoCopynumber.Windows(cel.f, output.path, analysis.path,
                cdf.file, chrX.probes, chrY.probes, analysis.qca, 
                analysis.qcc, analysis.snp, analysis.annotdb, analysis.refmodel)
        } else {
            call <- CytoCopynumber.Mac(cel.f, output.path, analysis.path,
            cdf.file, chrX.probes, chrY.probes, analysis.qca, 
            analysis.qcc, analysis.snp, analysis.annotdb, analysis.refmodel)
        }
    }


    if (!is.null(call)) {
        ScreenAffy("Calling APT")
        system(call)
    }
}


