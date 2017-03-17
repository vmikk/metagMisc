

## shorten species name
# based on vegan::make.cepnames (1b816c1; Aug 2, 2011)
abbreviate_taxa_names <- function(names, nlet = 3, totl = 7, sep="_", seconditem = F){
  # names = character vector of species names
  # nlet = number of letters to take from genus-species parts
  # totl = number of letters in the final abbreviation
  # sep = which separator to use for genus-species case
  # seconditem =  take the second item of the original name for abbreviation

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


x <- c("Laccaria laccata", "Meliniomyces bicolor",
       "Inocybe cincinnata", "Inocybe", "Tylospora asterophora",
       "Cadophora finlandica", "Saccharomycetales")
abbreviate_taxa_names(x, nlet = 3, totl = 7, sep = "_")
abbreviate_taxa_names(x, nlet = 4, totl = 8, sep = "")   # same as vegan::make.cepnames
vegan::make.cepnames(x)
