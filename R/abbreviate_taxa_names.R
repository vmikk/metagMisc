
#' @title Abbreviate taxa names
#' @param names Character vector of species names
#' @param nlet Number of letters to take from genus-species parts
#' @param totl Number of letters in the final abbreviation
#' @param sep Which separator to use for genus-species case (default is underscore)
#' @param seconditem Logical, take the second item of the original name for abbreviation
#'
#' @description This function is based on vegan::\code{\link{make.cepnames}} (git 1b816c1; Aug 2, 2011).
#' @return Character vector of abbreviated taxa names
#' @export
#' @author Author of the original function is Jari Oksanen (\code{\link{make.cepnames}})
#' @seealso \code{\link{make.cepnames}}
#' @examples
#' x <- c("Laccaria laccata", "Meliniomyces bicolor",
#'   "Inocybe cincinnata", "Inocybe", "Tylospora asterophora",
#'   "Cadophora finlandica", "Saccharomycetales")
#'
#' abbreviate_taxa_names(x, nlet = 3, totl = 7, sep = "_")
#' abbreviate_taxa_names(x, nlet = 4, totl = 8, sep = "")   # same as vegan::make.cepnames
#' vegan::make.cepnames(x)
#'
abbreviate_taxa_names <- function(names, nlet = 3, totl = 7, sep="_", seconditem = F){

  ## make valid names
  names <- make.names(names, unique = FALSE)

  ## remove trailing and duplicated dots
  names <- gsub("\\.[\\.]+", ".", names)
  names <- gsub("\\.$", "", names)

  ## split by dots and take `nlet` letters of each element (if several)
  names <- lapply(
    strsplit(names, "\\."), function(x){
    if(length(x) > 1){
      substring(x, 1, nlet)
    } else {
      x
    }})

  ## Take first and last element or `totl` characters if only one element
  names <- unlist(lapply(
    names, function(x){
      if(length(x) > 1){
        paste(x[c(1, if(seconditem) 2 else length(x))], collapse = sep)
      } else {
        x
      }}))

  names <- abbreviate(names, minlength = totl)

  ## Final clean-up
  names <- make.names(names, unique = TRUE)
  return(names)
}
