
# CytoPreprocess
# -----------------------------------------------------------------------------
# Created: 24/12/2013

CytoProcess.SnpMatrix <- function(cychp.file, output.temp, ann.p) {

  name <- gsub(".cyhd.cychp.txt", "", basename(cychp.file))
  ScreenAffy("Calling individual '", name ,"'.")

  gtype <- data.frame( ProbeSetName=ann.p, GType=NA, stringsAsFactors=FALSE )
  rownames(gtype) <- gtype$ProbeSetName

  # Call the python parser's to clean the data
  CytoParser(cychp.file, output.temp, short=TRUE)
  # /Call the python parser's to clean the data

  # Create files' names
  genotype.file <- file.path(output.temp,
  paste0(name, ".cyhd.cychp.genotype.f.txt"), fsep=.Platform$file.sep)
  # /Create files' names

  # Load data
  ScreenAffy("Loading calls ('", name ,"').")
  continue <- TRUE
  gnds <- tryCatch({
    fread( genotype.file, header=TRUE )#, as.is=TRUE )
  }, warnning=function(w) {
    OptionsAffy.AddError(error.name = genotype.file,
      "Invalid file provided in 'genotype.file' (",
      genotype.file ,") some error from 'python'.")
      continue <- FALSE
  }, error=function(e) {
    OptionsAffy.AddError(error.name = genotype.file,
      "Invalid file provided in 'genotype.file' (",
      genotype.file ,") some error from 'python'.")
      continue <- FALSE
  })

  if (!continue) {
    ScreenAffy("Failed calling for individual '", name ,"'.")
    return(NA)
  }
  # /Load data

  # Process
  ScreenAffy("Creating individual calls set ('", name ,"').")
  tryCatch({
    rownames(gnds) <- gnds$ProbeSetName
    gtype[rownames(gnds), ] <- gnds$Call
    gtype[gtype=="AA"] <- 1
    gtype[gtype=="AB"] <- 2
    gtype[gtype=="BB"] <- 3
    gtype[gtype=="NC"] <- NA
  }, warning=function(w) {
    warning(w)
    return( rep(NA, length(ann.p) ) )
  }, error=function(e) {
    warning(e)
    return( rep(NA, length(ann.p) ) )
  })


  # tryCatch({
  # rownames(gnds) <- gnds$ProbeSetName;
  # calls <- gnds$Call;
  # calls[calls=="AA"] <- 1
  # calls[calls=="AB"] <- 2
  # calls[calls=="BB"] <- 3
  # names(calls) <- gnds$ProbeSetName;
  # }, warnning=function(w) {
  # OptionsAffy.AddError(error.name = genotype.file,
  # "Error on naming calls (",
  # genotypeset ,").")
  # continue <- FALSE
  # }, error=function(e) {
  # OptionsAffy.AddError(error.name = genotype.file,
  # "Error on naming calls  (",
  # genotypeset ,").")
  # continue <- FALSE
  # })

  # if (!continue) {
  # ScreenAffy("Failed naming calls for individual '", name ,"'.")
  # return(NA)
  # }

  # c <- list();
  # c$genotypes <- calls;
  # c$names <- gnds$ProbeSetName;

  # return(c);
  return(gtype[ann.p, 2])
  # /Process
  }
