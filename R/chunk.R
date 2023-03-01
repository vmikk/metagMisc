
# Split a vector into chunks equal size
# x = vector; n = number of chunks
chunk <- function(x, n){
  if(n > 1) { res <- split(x, cut(seq_along(x), n, labels = FALSE)) }
  if(n == 1){ res <- list(); res[[1]] <- x }
  return(res)
}

## for integer vector
# n = 3
# split(x, sort(x%%n))


## Split data.table into chunks
chunk_table <- function(x, n){
  if(n == 1){ res <- list(); res[[1]] <- x }
  if(n > 1) {
    cc <- chunk(1:nrow(x), n = n)
    res <- llply(.data = cc, .fun = function(r){ x[ r, ] })
  }
  return(res)
}
