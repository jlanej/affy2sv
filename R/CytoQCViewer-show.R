setMethod( "show", 
  signature("CytoQCViewer"),
  function(object) {
    cat("Object of class Cyto2QCView\n")
    cat("Visualization '", object@type, "'\n")
    if (object@type == "snp") {
        cat(" for all individuals in '", object@location, "'\n")
    }
    else {
        cat(" for all individual '", object@individual, "'\n")
    }
  }
)