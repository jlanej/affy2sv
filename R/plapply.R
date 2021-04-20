
# plapply
# -----------------------------------------------------------------------------
# Created: 24/12/2013

plapply <- function(X, FUN, ...) {
    if (exists("mclapply") & (Sys.info()["sysname"] != "Windows")) {
        if (length(X) < affy.mclMaxElements) {
            result <- mclapply(X, FUN, ..., mc.cores = affy.cores)
        }
        else {
            warning("More than ", affy.mclMaxElements, 
                " samples to be processed. 'mclapply' limited the transfer",
                " to 2GB from the background to foreground. In order to ",
                " use the 'mclapply' the set will be processed in batchs of ",
                affy.mclMaxElements, " samples.")
            result <- list()
            for(ii in seq(1, length(X), affy.mclMaxElements)) {
                endp <- ii + affy.mclMaxElements -1
                if (endp > length(X)) {
                    endp <- length(X)
                }
                ScreenAffy("Processing samples from ", ii, " to ", endp, ".")
                result <- c(result, mclapply(X[ii:endp], FUN, ..., mc.cores = affy.cores))
            }
        }
    } else {
        result <- lapply(X, FUN, ...)
    }
    return(result)
}