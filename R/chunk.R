
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
#' @param group_col Column name to split by (character)
#' @param to_list Logical, return a list of data.tables (default: FALSE)
#'
#' @return data.table or list (of the length defined by the number of chunks) with parts of the data.table
#' @export
#'
#' @examples
#' chunk_table(x = data.table(Letter = letters[1:12]), n = 4)
#'
chunk_table <- function(x, group_col, n, to_list = FALSE){

  ## Convert to data.table if not already
  if(! inherits(x = x, what = "data.table")){
    setDT(x)
  }

  ## Get group sizes
  group_sizes <- x[, .N, by = get(group_col)]
  setnames(group_sizes, "get", group_col)
  group_sizes[, cumsum_size := cumsum(N)]
  total_size <- group_sizes[, sum(N)]
  chunk_size <- ceiling(total_size / n)
  
  ## Assign chunk IDs
  group_sizes[, chunk_id := pmin(ceiling(cumsum_size / chunk_size), n)]

  ## Join back to original data
  x[group_sizes, chunk_id := i.chunk_id, on = group_col]
  
  ## Return as list or data.table based on to_list argument
  if(to_list) {
    ## Split into list and remove chunk_id column from each chunk
    res <- split(x, by = "chunk_id", keep.by = FALSE)
    return(res)
  } else {
    ## Return data.table with chunk_id column
    return(x[])
  }
}
