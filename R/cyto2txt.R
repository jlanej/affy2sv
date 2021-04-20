
cychp2txt <- function(cychp.path, output.path) {
    
    unlink(output.path, recursive=TRUE, force=TRUE)
    dir.create(output.path, showWarnings=FALSE)

    cychp.files <- list.files(cychp.path, full.names = TRUE, pattern = ".cychp")

    papply(cychp.files, callAPT, output.path = output.path)
}


callAPT <- function(filename.full, output.path) {
    
    call <- NULL
    os <- NULL
    if(Sys.info()["sysname"] == "Linux") {
        os <- "linux"
        apt.txt <- system.file(paste0("exec", .Platform$file.sep, 
            "linux64.apt-chp-to-txt"), package="affy2sv")
    } else {
        if(Sys.info()[ "sysname" ] == "Windows") {
            os <- "windows"
            apt.txt <- system.file(paste0("exec", .Platform$file.sep, 
                "win64.apt-chp-to-txt"), package="affy2sv")
        } else {
            os <- "mac"
            apt.txt <- system.file(paste0("exec", .Platform$file.sep, 
                "mac64.apt-chp-to-txt"), package="affy2sv")
        }
    }

    call <- paste0(apt.txt , " ", filename.full, " -o ", output.path)
    ScreenAffy("Creating call to APT (windows).\n", call)
    system(call)

    return(NA)
}