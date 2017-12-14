
#' @title Transform (non-metric) dissimilarity matrix to a weighted Euclidean distance (metric).
#' @description This function constructs a weighted Euclidean distance that optimally approximates the dissimilarities.
#' Advantages of the Euclidean approach include the neat decomposition of variance and the ordination’s optimal biplot display (Greenacre, 2017).
#' @param datt Data frame with species abundance (species = columns, samples = rows)
#' @param dist_type Dissimilarity index (for the supported methods see \code{\link[vegan]{vegdist}}) or "other"
#' @param dst Distance matrix (of class 'dist') if dist_type == "other"
#' @param drop_species Logical; TRUE indicates removal of "unimportant" species which doesn't contribute to the sample differentiation
#' @param importance_percentile Percentile value for importance below which species are considered unimportant
#' @param show_plot Logical; if TRUE, plot of original dissimilarities vs. the obtained weighted Euclidean distances will be shown
#' @param ndim Number of dimensions; NULL by default, number of dimesions will be automatically determined
#' @param ... Additional arguments may be passed to \code{\link[vegan]{vegdist}})
#' @details
#' The code of the function is based on the work of Prof. Michael Greenacre (2017).
#'
#' It is possible to eliminate species that have little or no value in measuring difference between the samples (with 'drop_species = TRUE').
#' For this one need to specify a threshold value for species importance ('importance_percentile').
#' By default, 'importance_percentile = 0.02', which indicates that all species with importance below the 2nd percentile of the species importance distribution will be removed.
#'
#' Pre-calculated distance matrix can be passed as input to this function with 'dst' argument ('dist_type' should be set to "other"). However, species removal ('drop_species') will be impossible in this case.
#' @return Function 'dissimilarity_to_distance' returns a list with the following components:
#' \itemize{
#'   \item WEdist. Weighted Euclidean distance (class 'dist');
#'   \item sp_weights. Data frame with species weights;
#'   \item stress. Stress 1 measure which corresponds to the loss of the variance due to distance transformation (see \code{\link[smacof]{stress0}});
#'   \item dissim_plot. (if 'show_plot = TRUE') gglot object with the corresponding plot.
#' }
#' @export
#' @seealso \code{\link[smacof]{smacofConstraint}}, \code{\link[vegan]{vegdist}}, \code{\link{dissimilarity_to_distance_importance_plot}}
#' @references Greenacre, M. (2017), Ordination with any dissimilarity measure: a weighted Euclidean solution. Ecology, 98: 2293–2300. doi:10.1002/ecy.1937
#' @examples
#' library(vegan); library(smacof)
#' data(varespec)
#'
#' ## Estimate weighted Euclidean distance that optimally approximates Bray-Curtis dissimilarity
#' bc_to_eucl <- dissimilarity_to_distance(varespec, dist_type = "bray")
#' # ade4::is.euclid(bc_to_eucl$WEdist)
#'
dissimilarity_to_distance <- function(datt, dist_type = "bray", dst = NULL, drop_species = F, importance_percentile = 0.02, show_plot = T, ndim = NULL, ...){

  # require(vegan)
  # require(smacof)
  # require(ggplot2)

  ## Estimate dissimilarity
  if(dist_type != "other"){
    if(!is.null(dst)){ warning("Dissimilarity will be computed with ", dist_type, " method; user-provided matrix (dst) will be ignored.\n") }
    dst <- vegan::vegdist(x = datt, method = dist_type, ...)
  }
  if(dist_type == "other"){
    if(is.null(dst)){ stop("Error: dissimilarity matrix (dst) should be provided.\n") }
    if(!is.null(dst)){ if(!class(dst) %in% "dist"){ stop("Error: provided dissimilarity (dst) should be of class dist.\n") } }
    if(drop_species == TRUE){ stop("Error: unimportant species removal is not implemented for a custom dissimilarity matrix.\n") }
  }

  ## Estimate number of dimesions
  if(is.null(ndim)){
    ndim <- min(nrow(datt), ncol(datt)) - 1
  }

  ## Determine weights to fit the original dissimilarities into Euclidean distances
  # smacofConstraint automatically standardizes dissimilarities to have sum of squares n(n-1)/2
  # The diagonal weights (given in cstr$C) are fitted to these standardized dissimilarities
  cstr <- smacof::smacofConstraint(
          delta = as.matrix(dst),
          constraint = "diagonal",
          external = datt,
          constraint.type = "ratio",
          eps = 1E-8, ndim = ndim,
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
    dst <- vegan::vegdist(x = datt, method = dist_type, ...)

    ## Re-estimate the number of dimesions
    ndim <- min(nrow(datt), ncol(datt)) - 1

    ## Re-estimate transformation parameters
    cstr <- smacof::smacofConstraint(
          delta = as.matrix(dst),
          constraint = "diagonal",
          external = datt,
          constraint.type = "ratio",
          eps = 1E-8, ndim = ndim,
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

    dissim_plot <- ggplot2::ggplot(data = plt, aes(x = Dissim, y = Euclid)) +
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



#' @title Display species importance plot based on a weighted Euclidean distance transformation of dissimilarities.
#' @details Species importance measures the contribution of individual species to the overall dissimilarity between the samples.
#' @param x Result of \code{\link{dissimilarity_to_distance}} function
#' @param scaled_weights Logical; if TRUE, scaled weights will be plotted
#' @param mark_removed Logical; if TRUE, species that were removed from distance will be highlighted with a separate symbol on a plot
#'
#' @return ggplot-object.
#' @export
#' @seealso \code{\link{dissimilarity_to_distance}}
#' @examples
#' library(vegan); library(smacof)
#' data(varespec)
#'
#' # Estimate weighted Euclidean distance
#' bc_to_eucl <- dissimilarity_to_distance(varespec, dist_type = "bray", drop_species = T, importance_percentile = 0.05, show_plot = F)
#'
#' # Plot species weigths
#' dissimilarity_to_distance_importance_plot(bc_to_eucl)
#'
dissimilarity_to_distance_importance_plot <- function(x, scaled_weights = TRUE, mark_removed = TRUE){

  require(ggplot2)

  ## Extract data
  dtt <- x$sp_weights

  ## Reorder according to a species importance
  dtt <- dtt[order(dtt$Weight, decreasing = F), ]
  dtt$Species <- factor(dtt$Species, levels=dtt$Species)

  ## Prepare a plot
  if(scaled_weights == FALSE){ xx <- "Weight" }
  if(scaled_weights == TRUE) { xx <- "WeightScaled" }

  if(mark_removed == FALSE){
    plt <- ggplot(data = dtt, aes_string(x = xx, y = "Species"))
  }
  if(mark_removed == TRUE){
    plt <- ggplot(data = dtt, aes_string(x = xx, y = "Species", shape = "Importance")) +
             scale_shape_manual(
               values = c(16, 8),
               breaks = c("important", "not_important_removed"),
               labels = c("important", "removed"))
  }

  plt <- plt +
    geom_point(size = 2) +
    labs(x = "Weight", y = "Species")

  return(plt)
}


## Internal function to characterize ordination solution from its eigenvalues
.eigen_stats <- function(eig, colname = "Value"){
  # eig = vector of eigenvalues
  # colname = name of the column with resulting values

  ## Number of positive eigenvalues
  positives <- which(eig > 0)

  ## Variance of the Euclidean part
  tot_eucl <- (100 * sum(eig[positives]) / sum(abs(eig)))

  ## Percentages of first two eigenvalues relative to Euclidean part
  rel_eucl <- (100 * eig[1:2] / sum(eig[positives]))

  ## Percentages of first two eigenvalues relative to sum of all eigenvalues (McArdle, Anderson, 2001)
  rel_all <-  (100 * eig[1:2] / sum(eig))

  ## Prepare output table
  res <- rbind(
    data.frame(
      Parameter = "Number of eigenvalues",
      Value = length(eig),
      stringsAsFactors = F),
    data.frame(
      Parameter = "Number of positive eigenvalues",
      Value = length(positives),
      stringsAsFactors = F),
    data.frame(
      Parameter = "Number of negative eigenvalues",
      Value = length(eig) - length(positives),
      stringsAsFactors = F),
    data.frame(
      Parameter = "Percentage of Euclidean variance",
      Value = tot_eucl,
      stringsAsFactors = F),
    data.frame(
      Parameter = "Percentage of Non-Euclidean variance",
      Value = 100 - tot_eucl,
      stringsAsFactors = F),
    data.frame(
      Parameter = "First eigenvalue relative to Euclidean part, %",
      Value = rel_eucl[1],
      stringsAsFactors = F),
    data.frame(
      Parameter = "Second eigenvalue relative to Euclidean part, %",
      Value = rel_eucl[2],
      stringsAsFactors = F),
    data.frame(
      Parameter = "First two eigenvalues relative to Euclidean part, %",
      Value = sum(rel_eucl),
      stringsAsFactors = F),
    data.frame(
      Parameter = "First eigenvalue relative to sum of all eigenvalues, %",
      Value = rel_all[1],
      stringsAsFactors = F),
    data.frame(
      Parameter = "Second eigenvalue relative to sum of all eigenvalues, %",
      Value = rel_all[2],
      stringsAsFactors = F),
    data.frame(
      Parameter = "First two eigenvalues relative to sum of all eigenvalues, %",
      Value = sum(rel_all),
      stringsAsFactors = F)
  )

  colnames(res)[2] <- colname

  return(res)
}
