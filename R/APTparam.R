#' Function to create an object used to call Affymetrix Power Tools
#' 
#' This function creates an object of class \code{APTCopyNumberParam},
#' \code{APTGenoQcParam} or \code{APTProbeSetGenotypeParam}. Each one 
#' is used to run a specific command-tool of the Affymetrix Power Tools.
#' 
#' @export APTparam
APTparam <- function(type, ... ) {
	if(type=="cytoscan"){
		APTcytoscan(...)
	} else if(type=="axiom") {
		APTaxiom(...)
	} else {
		stop("Invalid content of 'type'. It must be 'cytoscan' or 'axiom'.")
	}
}

APTcytoscan <- function(cel.list, output.path, ...) {
	new("APTCopyNumberParam",
		cel.list = cel.list,
		output.path = output.path,
		...
	)
}


APTaxiom <- function(cel.list, output.path, ...) {
	params <- list()
	params$qc <- new("APTGenoQcParam",
		cel.list = cel.list,
		output.path = output.path,
		...
	)
	params$genotype <- new("APTProbeSetGenotypeParam",
		cel.list = cel.list,
		output.path = output.path,
		...
	)
	return(params)
}
