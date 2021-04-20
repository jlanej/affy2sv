setMethod( "show", 
  signature("APTCopyNumberParam"),
  function(object) {
    cat("Object of type 'APTCopyNumberParam'\n")
    cat("-----------------------------------\n")
    cat(" . Tool: ....", object@tool, "\n")
    cat(" . Level: ...", object@level, "\n")
    cat(" . OS: ......", object@os, "\n")
    cat(" . CEL(s): ..", object@cel.list, "\n")
    cat(" . Output: ..", object@output.path, "\n")
    cat(" . A. Path: .", object@analysis.path, "\n")
  }
)

setMethod( "show", 
  signature("APTGenoQcParam"),
  function(object) {
    cat("Object of type 'APTGenoQcParam'\n")
    cat("-------------------------------\n")
    cat(" . Tool: ....", object@tool, "\n")
    cat(" . Level: ...", object@level, "\n")
    cat(" . OS: ......", object@os, "\n")
    cat(" . CEL(s): ..", object@cel.list, "\n")
    cat(" . Output: ..", object@output.path, "\n")
	cat(" . A. Path: .", object@analysis.path, "\n")
  }
)

setMethod( "show", 
  signature("APTProbeSetGenotypeParam"),
  function(object) {
    cat("Object of type 'APTProbeSetGenotypeParam'\n")
    cat("-----------------------------------------\n")
    cat(" . Tool: ....", object@tool, "\n")
    cat(" . Level: ...", object@level, "\n")
    cat(" . OS: ......", object@os, "\n")
    cat(" . CEL(s): ..", object@cel.list, "\n")
    cat(" . Output: ..", object@output.path, "\n")
	cat(" . A. Path: .", object@analysis.path, "\n")
  }
)