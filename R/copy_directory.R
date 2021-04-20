#' Function to recursive copy of a directory
#' 
#' This function copies a full directory to another location
#' 
#' @export CopyDirectory
CopyDirectory <- function(from, to, overwrite = TRUE, verbose = FALSE) {

    show <- function(...) {
        if (verbose) {
            message("[ CopyDirectory ]: ", ... )
        }
    }

    iCopyDirectory <- function(from, to, overwrite = TRUE, verbose = FALSE) {
        show("Creating new '", to, "'.")
        dir.create(to, showWarnings=FALSE)

        # Get and prepare the content to be copies
        from.dirs <- basename(list.dirs(from))
        from.dirs <- from.dirs[from.dirs != basename(from)]

        from.files <- list.files(from)
        from.files <- from.files[!(from.files %in% from.dirs)]

        # Copy files in this level
        d <- sapply(from.files, function(x){
                complete.from <- file.path(from, x, fsep=.Platform$file.sep)
                complete.to <- file.path(to, x, fsep=.Platform$file.sep)
                show("Coping file '", complete.from, "' to '",
                    complete.to, "'.")
                file.copy(complete.from, complete.to)
            }
        )
        rm("d")

        # Look for deeper levels and apply the same function
        d <- sapply(from.dirs, function(x) {
                complete.from <- file.path(from, x, fsep=.Platform$file.sep)
                complete.to <- file.path(to, x, fsep=.Platform$file.sep)
                iCopyDirectory(complete.from, complete.to, overwrite=FALSE,
                    verbose=verbose)
            }
        )
        rm('d')
    }

    # Clean 'from' and 'to' content
    from <- basename(from)
    to <- basename(to)

    # Remove first level
    if (overwrite) {
        show("Removing '", to, "'.")
        unlink(to, recursive = TRUE, force=TRUE)
    }

    # Start copying
    iCopyDirectory(from, to)
}
