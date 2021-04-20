
# CY2APT
# -----------------------------------------------------------------------------
# 14/12/2013

#' @export
SmcBind <- function(smc.a, smc.b) {
    if ((class(smc.a) != "list")|is.null(smc.a$genotype)|is.null(smc.a$map)) {
        stop("Invalid argument for 'SmcBind' in 'smc.a'. Expected 'SnpMatrix Container'.")
    }
    if ((class(smc.b) != "list")|is.null(smc.b$genotype)|is.null(smc.b$map)) {
        stop("Invalid argument for 'SmcBind' in 'smc.b'. Expected 'SnpMatrix Container'.")
    }
    if ((class(smc.a$map) != "data.frame")|
        (class(smc.a$genotype) != "SnpMatrix")) {
        stop("Invalid argument for 'SmcBind'. Elements of 'smc.a' have not the correct type.")
    }
    if ((class(smc.b$map) != "data.frame")|
        (class(smc.b$genotype) != "SnpMatrix")) {
        stop("Invalid argument for 'SmcBind'. Elements of 'smc.b' have not the correct type.")
    }

    if (dim(smc.b$map)[1] < dim(smc.a$map)[1]) {
        smc.c <- smc.b
        smc.b <- smc.a
        smc.a <- smc.c
        rm("smc.c")
    }

    na <- rownames(smc.a$map)
    nb <- rownames(smc.b$map)
    smc.a$genotype <- smc.a$genotype[,na]
    smc.b$genotype <- smc.b$genotype[,nb]

    map <- smc.a$map
    to.add <- smc.b$map[!(nb %in% na), ]
    map <- rbind(map, to.add)

    gen <- smc.a$genotype
    if (dim(to.add)[1] > 0) {
        for (ii in 1:dim(to.add)[1]) {
            val <- rep(NA, dim(gen)[1])
            names(val) <- rownames(gen)
            gen <- cbind2(gen, new("SnpMatrix", val))
        }
    }
    colnames(gen) <- c(colnames(smc.a$genotype), rownames(to.add))
    gen <- rbind2(gen, smc.b$genotype)
    gen <- new("SnpMatrix",gen)


    return(list(map=map, genotype=gen))
}
