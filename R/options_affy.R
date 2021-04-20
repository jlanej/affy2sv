
OptionsAffy.Init <- function(fun.name, verbose) {
    OptionsAffy.FunName(fun.name)
    OptionsAffy.Verbose(verbose)
    OptionsAffy.MclMaxElements(150)
    OptionsAffy.Error.Init()
}

OptionsAffy.FunName <- function(fun.name) {
    affy.fun.name <<- fun.name
}

OptionsAffy.Verbose <- function(verbose) {
    affy.verbose <<- verbose
}

OptionsAffy.Cores <- function(cores) {
    affy.cores <<- cores
    if (exists("mclapply")) {
        options(cores=affy.cores)
    }
}

OptionsAffy.Error.Init <- function() {
    errors <<- list()
}

OptionsAffy.AddError <- function(error.name, ...) {
    errors[error.name] <<- paste0(...)
}

OptionsAffy.CheckErrors <- function() {
    if (length(errors) != 0) {
        stop(
            "Some errors succeed. Type 'OptionsAffy.ShowErrors' to check them")
    }
}

#' @export
OptionsAffy.ShowErrors <- function() {
    for (err in names(errors)) {
        message(err)
        message("\t", errors[err])
    }
}

OptionsAffy.MclMaxElements <- function(melms) {
    affy.mclMaxElements <<- melms
}