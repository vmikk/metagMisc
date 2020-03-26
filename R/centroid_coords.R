
## Function to find centroid coordinates
centroid_coords <- function(x, type = "centroid"){
  # x = data.frame with ordination scores

  ## Mean
  if(type == "centroid"){
    res <- plyr::numcolwise(mean)(x)
  }

  ## Spatial median (minimizes sum of distances of points from it)
  ## in 1d it is equal to the ordinary median (but not in >= 2d)
  if(type == "medoid"){
    res <- vegan:::ordimedian(x, groups = 1)
  }

  return(res)
}
