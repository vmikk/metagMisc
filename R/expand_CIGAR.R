expand_CIGAR <- function(cigar) {

  ## Insert 1 before any letter not preceded by a number
  cc <- gsub(pattern = "(?<![0-9])([A-Z])", replacement = "1\\1", x = cigar, perl = TRUE)

  ## Parse CIGAR
  cx <- GenomicAlignments::explodeCigarOps(cc)
  cn <- GenomicAlignments::explodeCigarOpLengths(cc)

  ## Expand each letter
  expn <- function(id){
    res <- mapply(rep, cx[[ id ]], cn[[ id ]])
    res <- paste(unlist(res), collapse = "")
    return(res)
  }

  ## Apply the function to all CIGAR strings
  res <- plyr::aaply(
    .data = seq_along(cigar),
    .margins = 1,
    .fun = expn)

  return(res)
}

