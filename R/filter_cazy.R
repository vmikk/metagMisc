

## Filter CAZY table
filter_cazy <- function(x, min_abund = 15, min_occurence = 10){
  abund <- colSums(x)
  occur <- apply(x, MARGIN = 2, FUN = function(z){ sum(z>0) })
  to_preserve <- which(abund >= min_abund & occur >= min_occurence)
  res <- x[to_preserve, ]
  return(res)
}