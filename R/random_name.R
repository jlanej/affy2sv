

# RandomName
# -----------------------------------------------------------------------------
# Created: 18/11/2013

RandomName <- function(prefix="", len=10) {
    numbers <- paste(sample(0:9, len, replace=TRUE), sep="", collapse="")
    return(paste(prefix, numbers, sep="", collapse=""))
}
