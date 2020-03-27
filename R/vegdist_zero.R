
## Function to handle samples with zero total abundances
vegdist_zero <- function(x, method = "bray", double_zero = 0, ...){

  ## Find samples with zero abundance
  rss <- rowSums(x, na.rm = TRUE) == 0
  
  ## All samples are non-zero
  if(!any(rss)){
    dd <- vegan::vegdist(x, method = method, ...)

  ## There are some samples with zero-abundance
  } else {

    zss <- which(rss)

    if(length(zss) == nrow(x)){
      stop("Sums for all samples are zero. Please check the data.\n")
    }

    ## Remove zero-samples
    tmp <- x[-zss, ]

    ## Estimate dissimilarity on non-zero samples
    dd <- vegan::vegdist(tmp, method = method, ...)

    ## How to substitute zero?
    if(!is.numeric(double_zero)){
      if(double_zero == "max"){ double_zero <- max(dd) 
      } else if(double_zero == "min"){ double_zero <- min(dd) 
      } else { stop("Error: Unknown 'double_zero' argument value.\n") }
    }

    ## Add zero samples to the distance matrix
    ## Dissimilarity between zero- and non-zero entries = 1
    ## Dissimilarity between zero- and zero entries = "double_zero"
    dd <- as.matrix(dd)
    dd <- cbind(dd, matrix(data = 1, nrow = nrow(dd), ncol = length(zss)))
    dd <- rbind(dd, matrix(data = 1, nrow = length(zss), ncol = ncol(dd)))

    rownames(dd) <- colnames(dd) <- c(rownames(tmp), names(zss))
    
    ## Find double zeros
    if(double_zero != 1){
      dd[names(zss), names(zss)] <- double_zero
    }

    ## Self-dissimilarity = 0
    diag(dd) <- 0

    ## Reorder samples as in original data
    dd <- dd[rownames(x), rownames(x)]
    dd <- as.dist(dd)
  } # end of rowSums == 0

  return(dd)
}
