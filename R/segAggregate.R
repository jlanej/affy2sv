
# segAggregate
# -----------------------------------------------------------------------------
# 22/05/2014
# FROM BY BLOG

#' @export
segAggregate <- function( x, colSam = 1, colChr = 2, colSegS = 3, colSegE = 4,
    sensLev = 500, aggFun = NULL ) {

    agg <- function( x, colId, colChr, colStart, colEnd, sensLev, aggFun ) {
        sensLev <- as.integer(sensLev)
        
        x <- x[ order( x[ , colChr ], x[ , colStart ] ), ]
        y <- x[ 1, ]

        if ( nrow( x ) < 2 ) {
            return( y )
        }

        yi <- 1
        for ( ii in 2:nrow( x ) ) {
            val <- abs( x[ ii, colStart ] - x[ ii -1 , colEnd ] )
            if ( val < sensLev ) {

                nx   <- x[ ii, c( colId, colChr, colStart, colEnd ) ]
                ji   <- 1
                elms <- 1:ncol( x )
                for( jj in elms[ -c( colId, colChr, colStart, colEnd ) ] ) {
                    nx <- cbind( aggFun[[ ji ]]( y[ yi, jj ], x[ ii, jj ]  ) )
                    ji <- ji + 1
                }

                y[ yi, colStart ] <- x[ ii, colEnd ]
            } else {
                y <- rbind( y, x[ ii, ] )
                yi <- yi + 1
            }
        }
        return( y )
    }
    
    if( is.null( aggFun ) | ( length( aggFun ) != ( ncol( x ) - 4 ) ) ) {
        message("Argument 'aggFun' is missing. Default will be used (numeric: mean, other: substitution).")
        type <- unlist( lapply( x[, -c( colSam, colChr, colSegS, colSegE ) ],
            class ) )
        aggFun <- list()
        ni <- 0
        for( ii in type ) {
            ni <- ni + 1
            if ( ii == "numeric" ) {
                aggFun[[ ni ]] <- mean
            } else {
                aggFun[[ ni ]] <- function( x, y ) { y }
            }
        }
    }

    y <- do.call( rbind, 
        lapply( unique( x[ , colSam ] ), function( sample_name ) {
            agg( x[ x[ , colSam ]  == sample_name, ], colSam, colChr, colSegS, 
                colSegE, sensLev, aggFun )
        })
    )
    return( y )
}