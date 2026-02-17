
#' Compute reverse complement of DNA/RNA sequences
#'
#' @description
#' Computes the reverse complement for each sequence in a character vector,
#' supporting standard nucleotides and common IUPAC ambiguity codes.
#' Useful for strand-aware sequence processing in amplicon and marker-gene workflows.
#'
#' @param x Character vector of nucleotide sequences to reverse-complement.
#'   Supports upper- and lower-case symbols, including ambiguous IUPAC bases.
#'
#' @details
#' Notes: Non-IUPAC characters are preserved as-is.
#'
#' @return Character vector of the same length as \code{x}, containing
#'   reverse-complemented sequences.
#'
#' @importFrom stringi stri_reverse
#' @export
rc <- function(x){
  # x = character vector

  s1 <- "ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv"
  s2 <- "TACGAtacgaNnRrYySsWwMmKkVvHhDdBb"

  res <- stringi::stri_reverse( chartr(old = s1, new = s2, x) )
  # res <- stringi::stri_reverse( stringx::chartr2(x = x, pattern = s1, replacement = s2) )
  return(res)
}

