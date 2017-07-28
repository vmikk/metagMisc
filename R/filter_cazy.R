
#' @title Filter CAZy table
#' @description This functions filters CAZy annotations based on the total abundance and occurence of CAZy genes across organisms. In other words, it removes CAZy annotations with low abundance.
#' @param x Data frame with gene counts (every column represents one CAZy annotation and every row represents one organism)
#' @param min_abund The minimum CAZy count across all organisms to be retained (default, 15)
#' @param min_occurence The minimum number of organisms with a given CAZy annotation to be retained (default, 10)
#'
#' @return Subsetted data frame
#' @export
#'
#' @examples
#'
filter_cazy <- function(x, min_abund = 15, min_occurence = 10){
  abund <- colSums(x)
  occur <- apply(x, MARGIN = 2, FUN = function(z){ sum(z>0) })
  to_preserve <- which(abund >= min_abund & occur >= min_occurence)
  res <- x[, to_preserve]
  return(res)
}
