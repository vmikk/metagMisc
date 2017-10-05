
#' @title Average multiple distance matrices.
#' @description This function can be used in order to avergage beta diversity over multiple rarefaction iterations.
#' @param dlist List of distance matrices (class 'dist' or 'matrix')
#'
#' @return Distance matrix of the same class as input matrices ('dist' or 'matrix') with averaged values.
#' @export
#'
#' @examples
#' library(plyr)
#' # Generate dummy data (list with 3 5x5 matrices)
#' ddd <- rlply(.n = 3, .expr = function(){ as.dist( matrix(sample(1:1000, size = 25, replace = T), ncol = 5)) })
#' ddd
#'
#' # Average matrices
#' mult_dist_average(ddd)
#'
mult_dist_average <- function(dlist){
  # dlist = list of `dist` or `matrix` objects

  ## Input = list of matrices
  if(class(dlist[[1]]) %in% "matrix"){
    res <- Reduce("+", dlist) / length(dlist)
  }

  ## Input = list of dist
  if(class(dlist[[1]]) %in% "dist"){
    res <- Reduce("+", llply(.data = dlist, .fun = as.matrix)) / length(dlist)
    res <- as.dist(res)
  }

  return(res)
}


## Compute beta diversity for each rarefaction iteration
mult_dissim <- function(x, method = "bray", average = T){
  # x = result of phyloseq_mult_raref (list of phyloseq objects)

  require(phyloseq)
  require(plyr)

  physeq_dissim <- llply(
    .data = x,
    .fun = function(z, ...){ phyloseq::distance(physeq = z, type = "samples", ...) },
    method = method,
    .progress = "text")

  if(average == TRUE){
    physeq_dissim_avg <- mult_dist_average(dlist = physeq_dissim)
    return(physeq_dissim_avg)
  } else {
    return(physeq_dissim)
  }
}

