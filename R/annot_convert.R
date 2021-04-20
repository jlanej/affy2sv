#' Function to update the MAP annotation of a \code{SnpMatrix Container}.
#' 
#' Given a \code{SnpMatrix Container}, this function returns a new 
#' \code{SnpMatrix Container} with its annotation (MAP slot) updates to
#' GRCh37.p13.
#' 
#' @export update_map
update_map <- function(smc) {
    if (is.na(match("map", names(smc)))) {
        stop("No 'map' object in given 'smc'.")
    }

    if (is.na(match("position", colnames(smc$map)))) {
        stop("Invalid 'map' in given 'smc', no 'position' field.")
    }

    if (is.na(match("snp.name", colnames(smc$map)))) {
        stop("Invalid 'map' in given 'smc', no 'snp.name' field.")
    }

    sel <- !smc$map$snp.name == "---"
    annot <- smc$map[sel,c(2,1)]

    mart <- useMart(biomart="snp", dataset="hsapiens_snp") ## GRCh37.p13
    
    results <- getBM(attributes=c("refsnp_id", "chr_name", "chrom_start"), filters="snp_filter", values=annot$snp.name, mart=mart)
    rm("annot", "sel", "mart")

    results <- results[!duplicated(results$refsnp_id), ]
    nmap <- smc$map[smc$map$snp.name %in% results$refsnp_id, ]
    nmap$position <- results$chrom_start
    nmap <- nmap[order(nmap$chromosome, nmap$position), ]

    return(list(genotype=smc$genotype, map=nmap))
}
