
#' @title Transfer rownames of data frame to a new column (for data export)
#'
#' @param x Data frame
#' @param name Name of the new column to be created
#' @param stringsAsFactors Logical; should character vectors be converted to factors? (default, FALSE)
#'
#' @return Data frame
#' @export
#'
#' @examples
#' x <- as.data.frame(matrix(sample(1:100, 100), nrow = 10, dimnames = list(paste("row", 1:10, sep=""), paste("col", 1:10, sep=""))))
#' x
#' dfRowName(x, name = "ColumnWithRowNames")
dfRowName <- function(x, name = "Rows", stringsAsFactors = FALSE){
  res <- data.frame(rownames(x), x, stringsAsFactors = stringsAsFactors)
  colnames(res)[1] <- name
  rownames(res) <- NULL
  return(res)
}
