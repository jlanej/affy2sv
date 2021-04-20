
# value
# ------------------------------------------------------------------------------
# 17/07/2013

value <- function( variable, options ) {
    item <- match( variable, options )
    if( !is.na( item ) ) {
        item <- options[ item ];
    }
    return( item );
}
