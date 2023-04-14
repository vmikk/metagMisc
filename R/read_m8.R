
#' @title Read BLAST output
#' @description
#' This function loads BLAST results saved in a tabular output format 6
#' (`-outfmt 6`; 12 columns).
#'
#' @param x Input file name
#' @param blast_colz Vector with column names (optional)
#' @param package Package to use for data preparation (default, "data.table")
#'
#' @return data.table / data.frame / tibble (depending on the `package` parameter)
#' @export
#'
#' @examples
#'
read_m8 <- function(x, blast_colz = NULL, package = "data.table"){

  ## Load data with built-in R commands
  if(package == "base"){
    if(grepl(pattern = "\\.gz$", x)){ x <- gzfile(x) }
    res <- read.delim(file = x, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }

  ## Load data with `data.table` package
  if(package == "data.table"){
    res <- fread(file = x, header = FALSE, sep = "\t")
  }

  ## Load data with `readr` package
  if(package == "readr"){
    res <- readr::read_tsv(file = x, col_names = FALSE)
  }

  ## Add header
  ## Standard BLAST m8 format (-outfmt 6)
  if(is.null(blast_colz)){
    blast_colz <- c(
      "QueryName", "MatchID", "SeqIdentity", "AlignLen", "MismatchN", "GapOpenings",
      "QueryStart", "QueryEnd", "TargetStart", "TargetEnd", "Evalue", "BitScore")
  } # otherwise use user-provided column names

  colnames(res) <- blast_colz
  return(res)
}
