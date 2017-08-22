
#' @title Convert distance matrix to data frame
#' @description This function takes a distance matrix (of class 'dist') and transforms it to a data.frame, where each row represents a single pairwise comparison.
#' @param dist Distance matrix (object of class 'dist')
#' @param tri Logical, if TRUE - only lower triangular part of dist will be returned
#' @return Data frame
#' @export
#' @seealso \code{\link[spaa]{dist2list}} from spaa package.
#' @author The code is based on \code{\link[spaa]{dist2list}} by Jinlong Zhang.
#'
#' @examples
#'
dist2list <- function (dist, tri=TRUE) {
  if (!class(dist) == "dist") { stop("Error: The input data must be a dist object.\n") }

  dat <- as.data.frame(as.matrix(dist))
  if (is.null(names(dat))) {
    rownames(dat) <- paste(1:nrow(dat))
  }
  value <- stack(dat)$values
  rnames <- rownames(dat)
  namecol <- expand.grid(rnames, rnames)
  colnames(namecol) <- c("col", "row")
  res <- data.frame(namecol, value)

  if(tri == TRUE){    # return only lower triangular part of dist
    res <- res[-which(upper.tri(as.matrix(dist), diag = T)), ]
  }

  return(res)
}
