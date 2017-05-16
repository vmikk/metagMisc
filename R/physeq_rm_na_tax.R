

# Remove unused tax levels
physeq_rm_na_tax <- function(physeq){
  # rm_all <- function(x) { Filter(function(x)!all(is.na(x)), df) }
  rm_all <- function(df) { df[, !apply(is.na(df), 2, all)] }
  tax_table(physeq) <- rm_all( tax_table(physeq) )
  return(physeq)
}


