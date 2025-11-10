
#' @title Abbreviate taxa names
#' @param names Character vector of species names
#' @param minlengths Numeric vector of length 2: c(genus_length, epithet_length). Default c(4, 4)
#' @param seconditem Logical, take the second item of the original name for epithet
#' @param uniqgenera Logical, enforce unique genera abbreviations
#' @param named Logical, preserve original names as names attribute
#' @param method Character, abbreviation method passed to `abbreviate` function ("left.kept" or "both.sides")
#'
#' @description 
#' Creates abbreviated names for taxonomic species similar to \code{\link[vegan]{make.cepnames}}.
#' Removes hyphens, validates names, and provides control over abbreviation methods.
#'
#' @return Character vector of abbreviated taxa names
#' @export
#' @author Author of the original function is Jari Oksanen (\code{\link[vegan]{make.cepnames}})
#' @seealso \code{\link[vegan]{make.cepnames}}
#' @examples
#' x <- c("Laccaria laccata", "Meliniomyces bicolor",
#'   "Inocybe cincinnata", "Inocybe", "Tylospora asterophora",
#'   "Cadophora finlandica", "Saccharomycetales",
#'   "Auricularia auricula-judae")
#'
#' abbreviate_taxa_names(x, minlengths = c(3, 4))
#' abbreviate_taxa_names(x, minlengths = c(4, 4))
#'
abbreviate_taxa_names <- function(names, minlengths = c(4, 4), seconditem = FALSE,
                                  uniqgenera = FALSE, named = FALSE, method) {
  
  ## Preserve original names
  orignames <- names
  
  ## Handle empty input
  if (length(names) == 0) {
    return(character(0))
  }
  
  ## Remove hyphens
  names <- gsub("-", "", names)
  
  ## Make valid names
  names <- base::make.names(names, unique = FALSE, allow_ = FALSE)

  ## Remove trailing and duplicated dots
  names <- gsub("\\.[\\.]+", ".", names)
  names <- gsub("\\.$", "", names)

  dt <- data.table(
    original_names  = orignames,
    validated_names = names)
  
  ## Split by dots and extract genus and epithet
  dt[, names_split := strsplit(validated_names, ".", fixed = TRUE)]
  dt[, gen := sapply(names_split, function(x) x[1])]
  
  ## Extract epithet
  dt[, epi := sapply(names_split, function(x) {
    if (seconditem) {
      if (length(x) >= 2) x[2] else ""
    } else {
      if (length(x) > 1) x[length(x)] else ""
    }
  })]

  ## Set abbreviation parameters
  glen <- minlengths[1]
  nmlen <- sum(minlengths)
  if(missing(method)){ method <- "left.kept" }

  ## Abbreviate genus
  gen_abbrev_all <- abbreviate(
    abbreviate(dt$gen, glen, use.classes = FALSE, strict = !uniqgenera),
    glen, use.classes = TRUE, method = method)
  dt[, gen_abbrev := ifelse(epi != "", gen_abbrev_all, gen)]
  
  ## Combine genus and epithet
  dt[, combined := paste0(gen_abbrev, epi)]
  
  ## Abbreviate combined names
  dt[, abbreviated := abbreviate(combined, nmlen, use.classes = FALSE)]
  
  ## Final abbreviation pass
  dt[, final_names := abbreviate(abbreviated, nmlen, use.classes = TRUE, 
                                 method = method, named = FALSE)]

  ## Extract final names
  names <- dt$final_names
  
  ## Final clean-up for uniqueness
  names <- base::make.names(names, unique = TRUE)
  
  ## Add original names as names attribute for the resulting vector
  if(named){ names(names) <- orignames }
    
  return(names)
}
