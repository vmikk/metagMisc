

## Reverse-complement DNA sequence (supporting IUPAC degeneracies)
rc <- function(x){
  # x = character vector

  s1 <- "ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv"
  s2 <- "TACGAtacgaNnRrYySsWwMmKkVvHhDdBb"

  res <- stringi::stri_reverse( chartr(old = s1, new = s2, x) )
  # res <- stringi::stri_reverse( stringx::chartr2(x = x, pattern = s1, replacement = s2) )
  return(res)
}

