
# Split a vector into chunks equal size
# x = vector; n = number of chunks
chunk <- function(x, n){
  if(n > 1) { res <- split(x, cut(seq_along(x), n, labels = FALSE)) }
  if(n == 1){ res <- list(); res[[1]] <- x }
  return(res)
}
