
# CytoParser
# -----------------------------------------------------------------------------
# Created: 24/12/2013

CytoParser <- function(filen.name, output.path, short=FALSE) {
    pypa <- system.file("exec/parser.py", package="affy2sv")
    if (short) {
        pypa <- paste0( "python ", pypa, " -i ",filen.name,
            " -o ", output.path, " -s")# " -c -s" )
    }
    else {
        pypa <- paste0("python ", pypa, " -i ", filen.name,
            " -o ", output.path)#, " -c")
    }

    ScreenAffy("Calling 'parser':\n", pypa)
    system(pypa)

    return(NA)
}
