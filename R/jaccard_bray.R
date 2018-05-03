
## Convert between Bray-Curtis and Jaccard dissimilarity
jaccard_bray <- function(x, to = "bray"){
  # x = vector or dist with dissimilarity values
  # to = convertion direction ("bray", from Jaccard to Bray-Curtis; "jaccard", from Bray-Curtis to Jaccard)

  ## Convert Bray-Curtis to Jaccard dissimilarity
  if(to == "jaccard"){
    res <- (2*x)/(1+x)
  }

  ## Convert Jaccard to Bray-Curtis dissimilarity 
  if(to == "bray"){
    res <- x/(2-x)
  }

  return(res)
}
