
# ScreenAffy
# -----------------------------------------------------------------------------
# 14/12/2013

ScreenAffy <- function( ... ) {
    if (affy.verbose) {
        message(paste0("[", affy.fun.name, "]: ", ...))
    }
}
