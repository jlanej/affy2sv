
# GwProcess.SnpMatrix 
# -----------------------------------------------------------------------------
# 09/01/2014

GwProcess.SnpMatrix <- function(cel.files) {
    
    ScreenAffy("Creating 'SnpSet'.")
    tryCatch({
            crlmm.out <<- crlmm::crlmm(cel.files)
        }, warnning=function(w){ 
            stop("Warning performing 'CRLMM'.\n", w)
        }, error=function(e){ 
            stop("Error performing 'CRLMM'.\n", e)
        }
    )

    ScreenAffy("Coercing calls to 'SnpMatrix'.")
    calls.ind <<- tryCatch({
            new("SnpMatrix", t(calls(crlmm.out)))# - 1L)
        }, warning=function(w){
            stop("Warning coercing to 'SnpMatrix'.\n", w)
        }, error=function(e){ 
            stop("Error coercing to SnpMatrix.\n", e)
        }
    )

    return(calls.ind)
}
