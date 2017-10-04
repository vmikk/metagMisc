

## Avergage beta diversity over all rarefaction iterations
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

