

dissimilarity_to_distance <- function(datt, dist_type = "bray", dst = NULL, drop_species = F, importance_percentile = 0.02, show_plot = T, ...){
  # datt = initial data (species = columns, samples = rows)
  # dist_type = "bray"
  # dst = distance matrix if dist_type == "other"
  # drop_species = Logical, remove unimportant species
  # importance_percentile
  # ... passed to vegan::vegdist

  require(vegan)
  require(smacof)
  # require(ggplot2)

  ## Estimate dissimilarity
  if(dist_type != "other"){
    if(!is.null(dst)){ warning("Dissimilarity will be computed with ", dist_type, " method; user-provided matrix (dst) will be ignored.\n") }
    dst <- vegdist(x = datt, method = dist_type, ...)
  }
  if(dist_type == "other"){
    if(is.null(dst)){ stop("Error: dissimilarity matrix (dst) should be provided.\n") }
    if(!is.null(dst)){ if(!class(dst) %in% "dist"){ stop("Error: provided dissimilarity (dst) should be of class dist.\n") } }
    if(drop_species == TRUE){ stop("Error: unimportant species removal is not implemented for a custom dissimilarity matrix.\n") }
  }

  ## Determine weights to fit the original dissimilarities into Euclidean distances
  # smacofConstraint automatically standardizes dissimilarities to have sum of squares n(n-1)/2
  # The diagonal weights (given in cstr$C) are fitted to these standardized dissimilarities
  cstr <- smacofConstraint(
          delta = as.matrix(dst),
          constraint = "diagonal", 
          external = datt, 
          constraint.type = "ratio",
          eps = 1E-8, ndim = ncol(datt),
          verbose = FALSE)

  ## To get the weights C for the unstandardized dissimilarities a correction factor corfact is needed
  n <- nrow(datt)
  corfact <- (0.5 * sum(as.matrix(dst)^2)/(n*(n-1)/2))^0.5 
  C <- cstr$C * corfact
  
  ## Weights
  sp_weights <- abs(diag(C))        # the square roots of the w_j^2 in the article
  # sp_weights_2 <- sp_weights^2    # w_j^2

  ## Prepare species weights for export
  sp_weights_exp <- data.frame(
    Species = colnames(datt),
    Weight = sp_weights,
    stringsAsFactors = F)

  ## Add scaled weight [0-1]
  sp_weights_exp$WeightScaled <- with(sp_weights_exp, Weight / sum(Weight))


  ## Remove unimportant species that have little or no value in measuring difference between the samples
  if(drop_species == TRUE){

    ## Importance threshold
    imp_trsh <- quantile(x = sp_weights, probs = importance_percentile)

    ## Find species with importance below the threshold
    unimp_sp <- which(sp_weights < imp_trsh)
    cat("..Species selection step: ", length(unimp_sp), " species will be removed.\n")

    ## Add indicator column for species that were removed from data
    sp_weights_exp$Importance <- "important"
    sp_weights_exp$Importance[unimp_sp] <- "not_important_removed"

    ## Subset data
    datt <- datt[, -unimp_sp]

    ## Re-estimate dissimilarity
    dst <- vegdist(x = datt, method = dist_type, ...)

    ## Re-estimate transformation parameters
    cstr <- smacofConstraint(
          delta = as.matrix(dst),
          constraint = "diagonal", 
          external = datt, 
          constraint.type = "ratio",
          eps = 1E-8, ndim = ncol(datt),
          verbose = FALSE)

    ## Re-estimate the weights and correction factor
    corfact <- (0.5 * sum(as.matrix(dst)^2)/(n*(n-1)/2))^0.5 
    C <- cstr$C * corfact
    sp_weights <- abs(diag(C))
  }

  ## Estimate weighted Euclidean distances
  WEdist <- dist(as.matrix(datt) %*% diag(sp_weights))
  
  ## Stress 1 measure = loss of the variance due to distance transformation (stress * 100 = %)
  stress <- cstr$stress

  ## The same as 
  # D1 <- as.numeric(dst); D2 <- as.numeric(WEdist)
  # sqrt(sum((D1-D2)^2) / sum((D1)^2))

  ## Compare dissimilarities to weighted Euclidean distances
  if(show_plot == TRUE){

    require(ggplot2)

    plt <- data.frame(Dissim = as.numeric(dst), Euclid = as.numeric(WEdist))
    
    dissim_plot <- ggplot(data = plt, aes(x = Dissim, y = Euclid)) + 
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=Inf, color="grey", linetype = "longdash") +
      geom_point() + 
      labs(x = "Original dissimilarity", y="Weighted Euclidean distance") + 
      ggtitle(paste("Loss of the variance due to distance transformation = ", round(stress * 100, 3), "%", sep=""))

    print(dissim_plot)
  }

  ## Prepare output
  RES <- list()
  RES$WEdist <- WEdist                 # weighted Euclidean distance
  RES$sp_weights <- sp_weights_exp     # species weights
  RES$stress <- stress                 # stress value

  if(show_plot == TRUE){ RES$dissim_plot <- dissim_plot }

  return(RES)
}
