#' Number formating to string
#' 
#' This function converts a number to srting using "ten annotation".
#' 
#' @export TenNotation
TenNotation <- function(X, ndec=2) {
    exp <- floor(log10(X))
    exp.abs <- abs(exp)
    man <- X * 10**exp.abs
    man <- round(man, ndec)
    man <- str_pad(man, width=ndec, side="left", pad="0")
    return(paste0(man, "*10^", exp))
}
