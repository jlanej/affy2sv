#' Unload \code{affy2sv}
#' 
#' Shortcut to \code{detach( "package:affy2sv", unload=TRUE )}.
#' 
#' @export affy2sv.unload
affy2sv.unload <- function() {
    detach( "package:affy2sv", unload=TRUE )
}

#' Close the current R session
#' 
#' Shortcut to \code{q('no')}.
#' 
#' @export exit
exit <- function() {
    q('no')
}
