
#' @title Expand CIGAR Strings
#' @description
#' This function expands compact CIGAR (Compact Idiosyncratic Gapped Alignment Report) strings
#' into a full-length sequence of CIGAR codes ("operations").
#'
#' @param cigar A character vector containing CIGAR strings to be expanded.
#'
#' @return A character vector with the expanded CIGAR strings.
#' @export
#' @examples
#' expand_CIGAR("5MI2D")
#' expand_CIGAR(c("5MI2D", "3M2D3M"))
expand_CIGAR <- function(cigar) {

  ## There are no CIGAR strings for cluster representatives,
  ## instead, there is an "*" symbol
  ## Need to replace it
  if("*" %in% cigar | any(is.na(cigar))){
    asterisk <- TRUE
    asterisk_ids <- which(cigar %in% "*")
    cigar[ asterisk_ids ] <- "X"
  } else {
    asterisk <- FALSE
  }

  ## Handle missing data
  if(any(is.na(cigar))){
    nas <- TRUE
    nas_ids <- which(is.na(cigar))
    cigar[ nas_ids ] <- "X"
  } else {
    nas <- FALSE
  }

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

  ## Recover asterisk symbols and NAs
  if(asterisk == TRUE){ res[ asterisk_ids ] <- "*" }
  if(nas == TRUE)     { res[ nas_ids ]      <- NA  }

  names(res) <- NULL
  return(res)
}

