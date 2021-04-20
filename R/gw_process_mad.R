
# CreateBafLrrContainer
# -----------------------------------------------------------------------------
# 09/01/2014

CreateBafLrrContainer <- function(cnSet, oligoList) {
    ScreenAffy("Creating BAF/LRR/calls container.")
    res <- list()
    res$calls <- calls(cnSet)
    res$calls <- res$calls[grep("SNP", rownames(res$calls)), ]
    b <- baf(oligoList)
    l <- lrr(oligoList)
    p <- position(oligoList)
    n <- featureNames(oligoList)

    flatting <- function(ii, baf_lrr_list) {
        ds <- ffdf(baf_lrr_list[[ii]],
            ff(ii, length=length(rownames(baf_lrr_list[[ii]]))))
        rownames(ds) <- rownames(baf_lrr_list[[ii]])
        colnames(ds) <- c(colnames(baf_lrr_list[[ii]]), "chr")
        return(ds);
    }

    b <- sapply(1:length(b) , flatting, b)
    ScreenAffy("Filtering 'BAF'.");
    res$baf <- data.frame(b[[1]])
    for(ii in 2:length(b)) {
        res$baf <- rbind(res$baf, data.frame(b[[ii]]))
    }
    res$baf <- res$baf[rownames(res$baf) %in% rownames(res$calls), ]
    rm("b");

    ScreenAffy("Filtering 'LRR'.")
    l <- sapply(1:length(l), flatting, l);
    res$lrr <- data.frame(l[[1]])
    for(ii in 2:length(l)) {
        res$lrr <- rbind(res$lrr, data.frame(l[[ii]]))
    }
    res$lrr <- res$lrr[rownames(res$lrr) %in% rownames(res$calls), ]
    rm("l");

    ScreenAffy("Adding 'position' and 'name'.");
    res$pos <- vector()
    res$names <- vector()
    for(ii in 1:length(p)) {
        res$pos <- c(res$pos, p[[ii]])
        res$names <- c(res$names, n[[ii]])
    }
    names(res$pos) <- res$names
    res$names <- res$names[res$names %in% rownames(res$calls)]
    res$pos <- res$pos[names(res$pos) %in% rownames(res$calls)]
    rm("p", "n")

    ScreenAffy("Naming 'calls'.")
    c <- res$calls;
    l <- res$lrr;
    b <- res$baf;
    for(ii in 1:ncol(res$calls)) {
        res$calls[ , ii] <- c("AA", "AB", "BB")[c[ , ii]]
        res$lrr[ , ii] <- l[ , ii]/100
        res$baf[ , ii] <- NormalizeBaf(b[ , ii]/1000,
            res$calls[res$names , ii])
    }
    rm("c", "l", "b")

    res$calls <- res$calls[res$names, ]
    names(res) <- c("calls", "baf", "lrr", "pos", "names")

    return(res);
}


# NormalizeBaf
# -----------------------------------------------------------------------------
# 09/01/2014

NormalizeBaf <- function(baf, gtype) {
    # gtype <- ds[ , "GType" ]
    # baf <- as.numeric( as.vector( ds[ , "B.Allele.Freq" ] ) )

    frS <- rep( NA, length( gtype ) )
    result <- rep( NA, length( gtype ) )

    result[ gtype == "AA" ] <- baf[ gtype == "AA" ]
    result[ gtype == "AB" ] <- baf[ gtype == "AB" ]
    result[ gtype == "BB" ] <- baf[ gtype == "BB" ]

    # Apply normalization to the AA and BB
    thAA <- 0.05;
    thAB <- 0.5;
    thBB <- 0.95;
    result[ result < thAA ] <- 0;
    result[ result > thBB ] <- 1;
    result[ thAA <= result & result < thAB & gtype != "AB" & !is.na( result ) ] <-
        ( 0.5 * ( result[ thAA <= result & result < thAB & gtype != "AB" & !is.na( result ) ] - thAA ) ) / ( thAB - thAA );
    result[ thAB <= result & result < thBB & gtype != "AB" & !is.na( result ) ] <-
        0.5 + ( ( 0.5 * ( result[ thAB <= result & result < thBB & gtype != "AB"& !is.na( result ) ] - thAB ) ) / ( thBB - thAB ) );

    return( result );
}


# GwProcess.Mad
# -----------------------------------------------------------------------------
# 09/01/2014

GwProcess.Mad <- function(data, output.path, fill.24, verbose) {
    if (length(data) < 1) {
        ScreenAffy("No files to be saved?")
        return(NA)
    }

    ncel <- ncol(data$calls);
    for(ii in 1:ncel) {
        file.name <- paste0(gsub(".CEL", "", colnames(data$baf)[ii]))
        file.name <- file.path(output.path, file.name,
            fsep=.Platform$file.sep)

        madfile <- data.frame(data$names, data$baf[ , ncel + 1], data$pos,
            data$lrr[ , ii], data$baf[ , ii], data$calls[ , ii])
        colnames(madfile) <- c("Name", "Chr", "Position", "Log.R.Ratio",
            "B.Allele.Freq", "GType")

        madfile[madfile$Chr==23, "Chr"] <- "X";

        if (fill.24) {
            ScreenAffy("Filling chromosome 24 for '", file.name, "'.")
            lny <- round(sum(madfile$Chr == "X") / 3)
            cn <- madfile[madfile$Chr == "X", ][1:lny, ]
            cn[ , "Chr"] <- "Y"
            madfile <- rbind(madfile, cn)
        }

        ScreenAffy("Saving file '", file.name, "'.")
        write.table(madfile, file=file.name, quote=FALSE, row.names=FALSE,
            sep="\t", na="NA")
    }
}
