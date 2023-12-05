
#' Split a vector into chunks equal size
#'
#' @param x A vector (e.g., character or numeric vector)
#' @param n The number of chunks required (integer)
#'
#' @return List (of the length defined by the number of chunks) with vectors
#' @export
#'
#' @examples
#' chunk(x = letters, n = 4)
#'
chunk <- function(x, n){
  if(n > 1) { res <- split(x, cut(seq_along(x), n, labels = FALSE)) }
  if(n == 1){ res <- list(); res[[1]] <- x }
  return(res)
}

## for integer vector
# n = 3
# split(x, sort(x%%n))


#' Split data.table into chunks
#'
#' @param x A data.table
#' @param n The number of chunks required (integer)
#'
#' @return List (of the length defined by the number of chunks) with parts of the data.table
#' @export
#'
#' @examples
#' chunk_table(x = data.table(Letter = letters[1:12]), n = 4)
#'
chunk_table <- function(x, n){
  if(n == 1){ res <- list(); res[[1]] <- x }
  if(n > 1) {
    cc <- chunk(1:nrow(x), n = n)
    res <- llply(.data = cc, .fun = function(r){ x[ r, ] })
  }
  return(res)
}
