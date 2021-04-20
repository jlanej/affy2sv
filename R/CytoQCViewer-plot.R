setMethod( "plot", 
  signature("CytoQCViewer"),
  function(x, cychp.attime = 4, ...) {
    plot <- NA
    if (x@verbose) {
        cat("[plot] Type:", x@type, "\n")
    }

    OptionsAffy.Init("CytoQCViewer-plot", x@verbose)
    OptionsAffy.Cores(cychp.attime)

    output.path <- RandomName(prefix="tmp.cyto.plotQc.")
    dir.create(output.path, showWarnings=FALSE)

    if (x@type == "snp") {
        if (x@verbose) {
            cat("[plot] Location: ", x@location, "\n")
            cat("[plot] SNP: ", x@individual, "\n")
        }

        if (x@verbose) {
            cat("[plot] Parsing individuals.\n")
        }
        
        # PARSING ALL DATA-FILES
        files <- list.files(x@location, full.name=T, pattern=".cyhd.cychp.txt")
        
        
        parse <- function(x, verbose) {
            pypa <- system.file("exec/parser.py", package="affy2sv")
            pypa <- paste0("python ", pypa, " -i ", x, " -o ", output.path)
            if(verbose) {
                message("[plot][parser]", pypa)
            }
            system(pypa)
        }
        plapply(files, parse, verbose=x@verbose)
        # /

        plot <- plotSnpIntensities(output.path, x@individual, x@verbose)
    }
    if (x@type == "int") {
        gnds.path <- paste0(x@location, .Platform$file.sep, x@individual, collapse='', sep='')

        if (x@verbose) {
            cat("[plot] Parsing file '", gnds.path, "' to '", output.path, "'\n")
        }

        # PARSING SINGLE FILE
        pypa <- system.file("exec/parser.py", package="affy2sv")
        pypa <- paste0("python ", pypa, " -i ", gnds.path,
                 " -o ", output.path)
        system(pypa)
        # /

        f.f <- strsplit(x@individual, "\\.")[[1]]
        f.f <- paste0(paste0(f.f[-length(f.f)], collapse="", sep="."), "genotype.f.txt")
        f.f <- paste0(output.path, .Platform$file.sep, f.f, collapse="", sep="")

        if (x@verbose) {
            cat("[plot] Location: ", gnds.path, "\n")
            cat("[plot] Individual: ", x@individual, "\n")
            cat("[plot] Genotype: ", f.f, "\n")
        }

        gnds <- fread(f.f, header=TRUE)
        plot <- plotIndIntensities(gnds, x@individual, x@verbose)
    }
    if (x@type == "sc") {
        gnds.path <- paste0(x@location, .Platform$file.sep, x@individual, collapse='', sep='')

        if (x@verbose) {
            cat("[plot] Parsing file '", gnds.path, "' to '", output.path, "'\n")
        }

        # PARSING SINGLE FILE
        pypa <- system.file("exec/parser.py", package="affy2sv")
        pypa <- paste0("python ", pypa, " -i ", gnds.path,
                 " -o ", output.path)
        system(pypa)
        # /

        f.f <- strsplit(x@individual, "\\.")[[1]]
        f.f <- paste0(paste0(f.f[-length(f.f)], collapse="", sep="."), "genotype.f.txt")
        f.f <- paste0(output.path, .Platform$file.sep, f.f, collapse="", sep="")

        if (x@verbose) {
            cat("[plot] Location: ", gnds.path, "\n")
            cat("[plot] Individual: ", x@individual, "\n")
            cat("[plot] Genotype: ", f.f, "\n")
        }
        gnds <- fread(f.f, header=TRUE)
        plot <- plotIndStrength(gnds, x@individual, x@verbose)
    }

    if (x@verbose) {
        cat("[plot] Cleaning.\n")
    }
    unlink(output.path, recursive=TRUE, force=TRUE)

    return(plot)
  }
)

plotIndIntensities <- function(gnds, individual, verbose) {
    # Individual Intensities
    plot1 <- ggplot(data = gnds, aes(x=log2(ASignal), y=log2(BSignal), color=Call)) + geom_point(shape = 20, size=1)
    plot1 <- plot1 + scale_color_manual(name = "GType", labels = c("AA", "AB", "BB", "NC"), values = c("#008000", "#FF6347", "#87CEEB", "black"))
    plot1 <- plot1 + ggtitle(paste0("ASignal vs. BSignal for ", individual))
    #plot1 <- plot1 + expand_limits(x = 0, y = 0)
    return(plot1)
}

plotIndStrength <- function(gnds, individual, verbose) {
    # Individual Strength vs. Contrast
    plot2 <- ggplot(data = gnds, aes(x=Contrast, y=SignalStrength, color=Call)) + geom_point(shape = 20, size=1)
    plot2 <- plot2 + scale_color_manual(name = "GType", labels = c("AA", "AB", "BB", "NC"), values = c("#008000", "#FF6347", "#87CEEB", "black"))
    plot2 <- plot2 + ggtitle(paste0("Strength vs. Contrast for ", individual))
    #plot2 <- plot2 + expand_limits(x = 0, y = 0)
    return(plot2)
}

plotSnpIntensities <- function(path, individual, verbose) {
    # Population
    gnds.files <- list.files(path, full.names=T, pattern=".genotype.f.txt")
    
    dta <- mclapply(gnds.files, function(x) {
        message( x )
        gnds <- fread(x, header = TRUE)
        as.data.frame(gnds[gnds[, ProbeSetName==individual],])
    })
    names(dta) <- list.files("cychp/", pattern=".genotype.f.txt")
    
    dta <- data.table(do.call("rbind", dta))

    plot3 <- ggplot(data = dta, aes(x=log2(ASignal), y=log2(BSignal), color=Call)) + geom_point(shape = 20, size=1)
    plot3 <- plot3 + scale_color_manual(name = "GType", labels = c("AA", "AB", "BB", "NC"), values = c("#008000", "#FF6347", "#87CEEB", "black"))
    plot3 <- plot3 + ggtitle(paste0("ASignal vs. BSignal for ", individual))
    #plot3 <- plot3 + expand_limits(x = 0, y = 0)
    return(plot3)
}
