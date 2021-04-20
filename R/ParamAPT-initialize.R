setMethod( "initialize", 
  signature("APTCopyNumberParam"),
    function(.Object, cel.list, output.path, level = "standard", 
                analysis.path = "", cdf = "", chrX = "", chrY = "", 
                qca = "", qcc = "", snp = "", annot.db = "",
                refmodel = "", params = "") { 
      ## Fill internals
      if(Sys.info()["sysname"] == "Linux") {
        .Object@tool <- "linux64.apt-copynumber-cyto"
        .Object@os   <- "GNU/Linux"
      } else {
        if(Sys.info()[ "sysname" ] == "Windows") {
          .Object@tool <- "win64.apt-copynumber-cyto.exe"
          .Object@os   <- "Windows"
        } else {
          .Object@tool <- "mac64.apt-copynumber-cyto"
          .Object@os   <- "MAC OS"
        }
      }
      # /

      ## Check validity
      if(missing(level) | gsub("\\s", "", level) == "") {
        .Object@level <- "standard"
        warning("'level' was not defined, set to 'standard' by default")
      } else {
        if(!tolower(level) %in% c("standard", "custom")) {
          stop("invalid content for 'level', it must be 'standard' or 'custom'")
        }
        .Object@level <- tolower(level)
      }

	  if(!missing(cel.list)){
		  cel.list <- gsub("\\s","", cel.list)
		  if(cel.list == "") {
			return("'cel.list' must be filled with the path where the .CEL files are stored")
		  }
		  .Object@cel.list <- cel.list
	  }
	  if(!missing(output.path)){
		  output.path <- gsub("\\s","", output.path)
		  if(output.path == "") {
			return("'output.path' must be filled with a path where the APT output files will be stored while processed by affy2sv")
		  }
		  .Object@output.path <- output.path
	  }
	  
      params        <- gsub("\\s","", params)
      analysis.path <- gsub("\\s","", analysis.path)
      cdf           <- gsub("\\s","", cdf)
      chrX          <- gsub("\\s","", chrX)
      chrY          <- gsub("\\s","", chrY)
      qca           <- gsub("\\s","", qca)
      qcc           <- gsub("\\s","", qcc)
      snp           <- gsub("\\s","", snp)
      annot.db      <- gsub("\\s","", annot.db)
      refmodel      <- gsub("\\s","", refmodel)

      if(.Object@level == "standard") {
        if(params != "") {
          stop("while 'level' is 'standard', 'params' cannot be filled")
        }
        .Object@analysis.path <- analysis.path
        .Object@cdf <- cdf
        .Object@chrX <- chrX
        .Object@chrY <- chrY
        .Object@qca <- qca
        .Object@qcc <- qcc
        .Object@snp <- snp
        .Object@annot.db <- annot.db
        .Object@refmodel <- refmodel
        .Object@params <- ""
      } else {
        if(params == "") {
          stop("while 'level' is 'custom', 'params' must be filled")
        }
        .Object@analysis.path <- ""
        .Object@cdf <- ""
        .Object@chrX <- ""
        .Object@chrY <- ""
        .Object@qca <- ""
        .Object@qcc <- ""
        .Object@snp <- ""
        .Object@annot.db <- ""
        .Object@refmodel <- ""
        .Object@params <- params
      }

      if(.Object@level == "custom") {
        if(analysis.path != "" | cdf != "" | chrX != "" | chrY != "" | qca != "" | qcc != "" | snp != "" | annot.db != "" | refmodel != "") {
          stop("while 'level' is 'custom', 'analysis.path', 'cdf', 'chrX', 'chrY', 'qca', 'qcc', 'snp', 'annot.db' and 'refmodel' cannot be filled")
        }
      } else {
        if(analysis.path == "" | cdf == "" | chrX == "" | chrY == "" | qca == "" | qcc == "" | snp == "" | annot.db == "" | refmodel == "") {
          stop("while 'level' is 'standard', 'analysis.path', 'cdf', 'chrX', 'chrY', 'qca', 'qcc', 'snp', 'annot.db' and 'refmodel' must be filled")
        }
      }
      ## /
      return(.Object)
    }
)

setMethod( "initialize", 
  signature("APTGenoQcParam"),
    function(.Object, cel.list, output.path, level = "standard", 
                analysis.path = "", xml.config = "", params = "") { 
      ## Filling internals
      if(Sys.info()["sysname"] == "Linux") {
        .Object@tool <- "linux64.apt-geno-qc"
        .Object@os   <- "GNU/Linux"
      } else {
        if(Sys.info()[ "sysname" ] == "Windows") {
          .Object@tool <- "win64.apt-geno-qc.exe"
          .Object@os   <- "GNU/Linux"
        } else {
          .Object@tool <- "mac64.apt-geno-qc"
          .Object@os   <- "MAC OS"
        }
      }
      if(gsub("\\s", "", level) == "") {
        .Object@level <- "standard"
        warning("'level' was not defined, set to 'standard' by default")
      }
      ## /

      ## Check validity
      if(!level %in% c("standard", "custom")) {
        stop("valid content for 'level' must be 'standard' or 'custom'")
      }

      cel.list <- gsub("\\s","", object@cel.list)
      output.path <- gsub("\\s","", object@output.path)

      if(cel.list == "") {
        return("'cel.list' must be filled with the path where the .CEL files are stored")
      }
      if(output.path == "") {
        return("'output.path' must be filled with a path where the APT output files will be stored while processed by affy2sv")
      }

      .Object@cel.list <- cel.list
      .Object@output.path <- output.path

      params        <- gsub("\\s","", params)
      analysis.path <- gsub("\\s","", analysis.path)
      xml.config    <- gsub("\\s","", xml.config)

      if(.Object@level == "standard") {
        if(params != "") {
          return("while 'level' is 'standard', 'params' cannot be filled")
        }
        .Object@analysis.path <- analysis.path
        .Object@xml.config <- xml.config
        .Object@params <- ""
      } else {
        if(params == "") {
          return("while 'level' is 'custom', 'params' must be filled")
        }
        .Object@analysis.path <- ""
        .Object@xml.config <- ""
        .Object@params <- params
      }

      if(.Object@level == "custom") {
        if(analysis.path != "" | xml.config != "" ) {
          return("while 'level' is 'standard', 'analysis.path' and 'xml.config' cannot be filled")
        }
      } else {
        if(analysis.path == "" | xml.config == "" ) {
          return("while 'level' is 'standard', 'analysis.path' and 'xml.config' must be filled")
        }
      }
      ## /
      return(.Object)
    }
)

setMethod( "initialize", 
  signature("APTProbeSetGenotypeParam"),
    function(.Object, cel.list, output.path, level = "standard", 
                analysis.path = "", xml.config = "", params = "") { 
      ## Filling internals
      if(Sys.info()["sysname"] == "Linux") {
        .Object@tool <- "linux64.apt-geno-qc"
        .Object@os   <- "GNU/Linux"
      } else {
        if(Sys.info()[ "sysname" ] == "Windows") {
          .Object@tool <- "win64.apt-geno-qc.exe"
          .Object@os   <- "GNU/Linux"
        } else {
          .Object@tool <- "mac64.apt-geno-qc"
          .Object@os   <- "MAC OS"
        }
      }
      if(gsub("\\s", "", level) == "") {
        .Object@level <- "standard"
        warning("'level' was not defined, set to 'standard' by default")
      }
      ## /

      ## Check validity
      if(!level %in% c("standard", "custom")) {
        stop("valid content for 'level' must be 'standard' or 'custom'")
      }

      cel.list <- gsub("\\s","", object@cel.list)
      output.path <- gsub("\\s","", object@output.path)

      if(cel.list == "") {
        return("'cel.list' must be filled with the path where the .CEL files are stored")
      }
      if(output.path == "") {
        return("'output.path' must be filled with a path where the APT output files will be stored while processed by affy2sv")
      }

      .Object@cel.list <- cel.list
      .Object@output.path <- output.path

      params        <- gsub("\\s","", params)
      analysis.path <- gsub("\\s","", analysis.path)
      xml.config    <- gsub("\\s","", xml.config)

      if(.Object@level == "standard") {
        if(params != "") {
          return("while 'level' is 'standard', 'params' cannot be filled")
        }
        .Object@analysis.path <- analysis.path
        .Object@xml.config <- xml.config
        .Object@params <- ""
      } else {
        if(params == "") {
          return("while 'level' is 'custom', 'params' must be filled")
        }
        .Object@analysis.path <- ""
        .Object@xml.config <- ""
        .Object@params <- params
      }

      if(.Object@level == "custom") {
        if(analysis.path != "" | xml.config != "" ) {
          return("while 'level' is 'standard', 'analysis.path' and 'xml.config' cannot be filled")
        }
      } else {
        if(analysis.path == "" | xml.config == "" ) {
          return("while 'level' is 'standard', 'analysis.path' and 'xml.config' must be filled")
        }
      }
      ## /
      return(.Object)
    }
)