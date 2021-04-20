

# List.Intersect
# -----------------------------------------------------------------------------
# Update: 06/01/2014

List.Intersect <- function(list.vec) {
    return(Reduce(intersect, list.vec))
}
