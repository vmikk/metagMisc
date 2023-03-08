
## Function to read BLAST output
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
