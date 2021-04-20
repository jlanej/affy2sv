#' Function used to check the presence of python in the current system
#' 
#' This function runs a command-line sentence to check the rpresence of
#' python and pandas in the current system. The run command follows:
#' \code{python -c "import pandas"}.
#' 
#' @export CheckPython
CheckPython <- function() {
    t1 <- try(system('python -c "import pandas"', intern = TRUE))
    if (length(t1) != 0) {
        return(FALSE)
    }
    return(TRUE)
} 
