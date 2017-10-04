

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


## Compute beta diversity for each rarefaction iteration
mult_dissim <- function(x, method = "bray", avergage = T){
  # x = result of phyloseq_mult_raref (list of phyloseq objects)

  require(phyloseq)
  require(plyr)

  physeq_dissim <- llply(
    .data = x, 
    .fun = function(z, ...){ phyloseq::distance(physeq = z, type = "samples", ...) }, 
    method = method, 
    .progress = "text")

  if(avergage == TRUE){
    physeq_dissim_avg <- mult_dist_average(dlist = physeq_dissim)
    return(physeq_dissim_avg)
  } else {
    return(physeq_dissim)
  }
}

