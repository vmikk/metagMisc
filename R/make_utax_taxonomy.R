
#' @title Prepare taxonomy annotations in SINTAX or UTAX style.
#' @description This function takes taxonomy table and returns character string in USEARCH-style that can be used with SINTAX or UTAX taxonomy classificators.
#' @param x Data frame with one row (columns are ordered taxonomy levels)
#' @details Missing values are allowed and should be coded as 'NA'.
#' @return Character string.
#' @export
#' @seealso http://drive5.com/usearch/manual/tax_annot.html
#' @examples
#'
make_utax_taxonomy <- function(x){

  ## Default single letter codes specifying taxonomic levels
  taxranks <- data.frame(
    Rank = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    Abbr = c("d", "k", "p", "c", "o", "f", "g", "s"),
    stringsAsFactors = F)

  ## Remove non-taxonomic columns
  x <- x[, which(colnames(x) %in% taxranks$Rank)]

  ## Find missing taxonomy ranks (filled with NAs)
  missing_ranks <- is.na(x)

  if(any(missing_ranks)){

    ## All ranks are NAs
    if( sum(missing_ranks) == length(x) ){
      cat("Warning: all taxonomy ranks are missing.\n")
      res <- paste("tax=;", sep="")
      return(res)
    }

    ## If there is only one non-missing tax rank
    if( sum(missing_ranks) == length(x) - 1 ){
      # one_rank <- TRUE                                   # indicator variable
      one_rank_name <- colnames(x)[which(!missing_ranks)]  # which ranks is it?
    } else {
      # one_rank <- FALSE
    }

    ## Remove NAs
    x <- x[, -which(missing_ranks)]
  }

  ## Several non-missing tax ranks are available (normal case)
  if(length(x) > 1){
    ## Prepare rank names
    rnk <- taxranks[match(x = colnames(x), table = taxranks$Rank), "Abbr"]
  }

  ## Only one tax rank is present
  if(length(x) == 1){
    ## Prepare rank names
    rnk <- taxranks[match(x = one_rank_name, table = taxranks$Rank), "Abbr"]
  }

  ## Merge rank and taxonomy
  res <- paste(rnk, x, sep=":")

  ## Merge into one string
  res <- paste("tax=", paste(res, collapse = ","), ";", sep="")

  return(res)
}
