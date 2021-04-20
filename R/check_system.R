#' Function used to check the current system
#' 
#' This function checks than the current system is runs a 64bits OS. Then
#' checks fot python and pandas using \code{CheckPython}.
#' 
#' @export CheckSystem
CheckSystem <- function(bits=TRUE) {

    valid <- NULL

    if (bits & (length(intersect(Sys.info()["machine"], c("x86_64", "x86-64"))) == 0)) {
        valid <- paste0("Functions that need the execution of APT ",
            "required a 64bits system.")
    }

    #if (length(intersect(Sys.info()["sysname"], c("Linux", "Windows"))) == 0) {
    #    valid <- paste0("Version 1.0 of the package was only developed to be ",
    #        "used under Linux or Windows.")
    #}
    
    if (!CheckPython()) {
        valid <- paste0("Version 1.0 of the package needs a global/environment variable called 'python' to run the requested python scripts. Also the library 'numpy' and the library 'pandas'.")
    }

    if (!is.null(valid)) {
        OptionsAffy.FunName("CheckSystem")
        stop("[CheckSystem]: ", valid)
    }
}
