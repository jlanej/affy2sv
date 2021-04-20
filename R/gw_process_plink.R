

# GwProcess.SnpMatrix 
# -----------------------------------------------------------------------------
# 09/01/2014

GwProcess.Plink <- function(cel.files, common, allele.A, allele.B, 
    strand, output.name) {

    ScreenAffy("Creating 'SnpSet'.")
    tryCatch({
            crlmm.out <- crlmm::crlmm(cel.files)
        }, warnning=function(w){ 
            stop("Warning performing 'CRLMM'.\n", w)
        }, error=function(e){ 
            stop("Error performing 'CRLMM'.\n", e)
        }
    )


    ScreenAffy("Getting calls from 'SnpSet'.")
    calls <- calls(crlmm.out)
    
    ScreenAffy("Creating full 'call data-frame' from 'SnpSet'.")
    alleleA <- alleleB <- matrix("", nrow=nrow(calls), ncol=ncol(calls),
        dimnames=dimnames(calls))
    
    allA <- allele.A
    allComplA <- c("A","C","G","T")[match(allA, c("T", "G", "C", "A"))]
    forwA <- ifelse(strand=="+", allA, allComplA)
    
    allB <- allele.B
    allComplB <- c("A","C","G","T")[match(allB, c("T", "G", "C", "A"))]
    forwB <- ifelse(strand == "+", allB, allComplB)

    for(ii in 1:ncol(calls)) {
      alleleA[ ,ii] <- ifelse(calls[ ,ii] < 3, forwA, forwB)
      alleleB[ ,ii] <- ifelse(calls[ ,ii] < 2, forwA, forwB)
    }
    
    p1 <- seq(1, ncol(calls) * 2, 2)
    p2 <- seq(2, ncol(calls) * 2, 2)

    tped <- matrix("", ncol=ncol(calls) * 2, nrow=nrow(calls))
    for(ii in 1:ncol(calls)) {
        tped[ , p1[ii]] <- alleleA[ ,ii]
        tped[ , p2[ii]] <- alleleB[ ,ii]
    }

    rownames(tped) <- rownames(calls)

    ScreenAffy("Binding common data.")
    tped <- cbind(common[rownames(tped), ], tped)

    ScreenAffy("Writing 'TPED' file to disk '" , output.name, "'.")
    write.table(tped, file=output.name, quote=FALSE, row.names=FALSE, 
        col.names=FALSE)

    return(NA)

}
