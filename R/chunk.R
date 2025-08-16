
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
#' Splits a data.table into a specified number of chunks. When a grouping 
#' column is specified, all rows with the same group value are guaranteed to remain 
#' within the same chunk, preserving group integrity.
#'
#' @param x A data.table to be split into chunks
#' @param n Integer. The number of chunks to create (must be >= 1)
#' @param group_col Character. Optional column name for group-preserving splits. 
#'   When specified, all rows sharing the same value in this column will be kept 
#'   together in the same chunk. When NULL (default), performs simple row-based 
#'   chunking without group consideration.
#' @param to_list Logical. If TRUE, returns a list of data.tables (one per chunk) 
#'   without chunk identifiers. If FALSE (default), returns the original data.table 
#'   with an added 'chunk_id' column indicating chunk membership.
#'
#' @return When \code{to_list = FALSE}: Original data.table with added 'chunk_id' 
#'   column (integer values 1 to n). When \code{to_list = TRUE}: A list of 
#'   data.tables, each containing one chunk of the original data without the 
#'   'chunk_id' column.
#'   
#' @details 
#' Group-based chunking may result in uneven chunk sizes to preserve the constraint 
#' that all instances of each unique group value remain in the same chunk.
#'
#' @export
#'
#' @examples
#' library(data.table)
#' 
#' # Simple row-based chunking
#' dt <- data.table(Letter = letters[1:12], Value = 1:12)
#' chunk_table(dt, n = 4)
#' 
#' # Group-preserving chunking
#' dt2 <- data.table(Group = c(rep("A", 8), rep("B", 3), "C"), Value = 1:12)
#' chunk_table(dt2, n = 3, group_col = "Group")
#' 
#' # Return as list
#' chunks <- chunk_table(dt, n = 4, to_list = TRUE)
#' sapply(chunks, nrow)
#'
chunk_table <- function(x, n, group_col = NULL, to_list = FALSE){

  ## Check if n is a positive integer
  if(! is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)){
    stop("Error: n must be a positive integer.\n")
  }

  ## Convert to data.table if not already
  if(! inherits(x = x, what = "data.table")){
    setDT(x)
  }

  ## Edge case with n == 1 (return entire data.table as single chunk)
  if(n == 1) {
    if(to_list == TRUE) {
      ## Return as list with single element
      res <- list(x)
      return(res)
    } else {
      ## Add chunk_id column with value 1 for all rows
      x[, chunk_id := 1L]
      return(x[])
    }
  }

  ## TODO - in unbalanced groups, the number of chunks in output can be less than n

  ## If `group_col` is NULL, use simple row-based chunking
  if(is.null(group_col)) {

    total_rows <- nrow(x)
    chunk_size <- ceiling(total_rows / n)
    x[, chunk_id := ceiling(.I / chunk_size)]
    x[, chunk_id := pmin(chunk_id, n)]  # Ensure max chunk_id doesn't exceed n

  } else {
    ## Group-based chunking
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
  }
  
  ## Return as list or data.table based on to_list argument
  if(to_list == TRUE) {
    ## Split into list and remove chunk_id column from each chunk
    res <- split(x, by = "chunk_id", keep.by = FALSE)
    return(res)
  } else {
    ## Return data.table with chunk_id column
    return(x[])
  }
}
