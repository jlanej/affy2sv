
# GwProcess.PennCNV
# -----------------------------------------------------------------------------
# 09/01/2014

GwProcess.PennCNV <- function(data, output.name, verbose) {
    if (length(data) < 1) {
        ScreenAffy("No data to be saved?")
        return(NA)
    }

    ncel <- ncol(data$calls);

    ds <- data.frame(data$name, data$baf[ , ncel + 1], data$pos);
    names <- c("Name", "Chr", "Position");

    for (ii in 1:ncel) {
        ds <- cbind(ds, data$calls[ , ii], data$lrr[ , ii], data$baf[ , ii])

        cel.name <- gsub(".CEL", "", colnames(data$baf)[ii])
        names <- c(names, paste0(cel.name, ".GType"),
            paste0(cel.name, ".Log R Ratio"),
            paste0(cel.name, ".B Allele Freq"))
    }
    colnames(ds) <- names;

    ScreenAffy("Saving file '", output.name, ".penncnv'.")
    write.table(ds, file=paste0(output.name, ".penncnv"), quote=FALSE,
        row.names=FALSE, sep="\t", na="NA")
}
