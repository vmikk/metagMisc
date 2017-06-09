
#' @title Estimate multivariate standard error.
#' @description This function estimates a multivariate standard error as the residual mean square error from a PERMANOVA for a one-way model, including double resampling.
#' @param D A distance matrix among all samples
#' @param group The grouping vector
#' @param nresamp The number of re-samples (default, 10 000)
#' 
#' @details The routine calculates the means using the permutation approach, while the lower and upper quantiles are obtained using the bootstrapping approach including an adjustment for the bias in the bootstrap.
#' The main distinction from the original function (provied in the supplementary data in Anderson & Santana-Garcon, 2013) is that estimation is done not for increasing range of sample sizes, but only for one sample size common for all groups.
#' @return Data frame with the multivariate standard error averaged over resampling interations.
#' @references Anderson M.J., Santana-Garcon J. Measures of precision for dissimilarity-based multivariate analysis of ecological communities. Ecology Letters (2015) 18: 66-73. DOI 10.1111/ele.12385
#' @examples
#'
multSE <- function(D, group, nresamp = 10000) {

  require(plyr)

  # Ensure distance matrix is in the form of a matrix (rather than a "distance" object)
  D <- as.matrix(D)

  # Remove groups with only a single replicate (adapted from https://github.com/jslefche/multSE/)
  grp <- table(group)
  if(min(grp) == 1){
    groups_to_remove <- names(which(grp == 1))
    D <- D[!group %in% groups_to_remove, !group %in% groups_to_remove]
    group <- group[!group %in% groups_to_remove]
    warning("Groups with 1 replicate have been removed from the analysis!")
  }

  # Some necessary preliminary functions:
  MSE <- function(D, group) {
    D <- as.matrix(D)
    N <- dim(D)[1]
    g <- length(levels(factor(group)))
    X <- model.matrix(~factor(group))     # Model matrix
    H <- X %*% solve(t(X)%*%X) %*% t(X)   # Hat matrix
    I <- diag(N)                          # Identity matrix
    A <- -0.5*D^2
    G <- A - apply(A,1,mean) %o% rep(1,N) - rep(1,N) %o% apply(A,2,mean) + mean(A)
    MSE <- sum(diag((I-H) %*% G)) / (N-g)
    return(MSE)
  }

  quant.upper <- function(x) quantile(x, prob = 0.975, na.rm = TRUE)
  quant.lower <- function(x) quantile(x, prob = 0.025, na.rm = TRUE)

  # Getting parameters of the problem
  group <- factor(group)
  ng <- length(levels(group))   # number of groups
  n.i <- table(group)           # number of samples per group
  nmax <- min(n.i)              # smallest sample size across all groups
  N <- sum(n.i)                 # total number of samples
  index <- 1:N                  # sample index

  ## Double resampling function
  resamp <- function(nsub){
    ivec.p <- sample(index[group == levels(group)[1]], size = nsub, replace = FALSE)
    ivec.b <- sample(index[group == levels(group)[1]], size = nsub, replace = TRUE)
    for (i in 2:ng) {
       ivec.p <- c(ivec.p, sample(index[group == levels(group)[i]], size = nsub, replace = FALSE))
       ivec.b <- c(ivec.b, sample(index[group == levels(group)[i]], size = nsub, replace = TRUE))
    }
    group.resamp <- factor(rep(1:ng, each = nsub))
    D.perm <- D[ivec.p, ivec.p]
    D.boot <- D[ivec.b, ivec.b]
    res <- data.frame(
              multSE.store.p = sqrt(MSE(D.perm, group.resamp)/nsub),  # values under permutation resampling
              multSE.store.b = sqrt(MSE(D.boot, group.resamp)/nsub)   # values under bootstrap resampling
            )
    return(res)
  }

  # Repeat resampling for the smallest sample size across all groups
  permvals <- rdply(.n = nresamp, .expr = resamp(nmax), .progress = "text")

  # Estimate means and quantiles over resampling interations
  means <- mean(permvals$multSE.store.p)
  means.b <- mean(permvals$multSE.store.b)
  bias <- means - means.b
  lower <- quant.lower(permvals$multSE.store.b) + bias
  upper <- quant.upper(permvals$multSE.store.b) + bias

  # Calculation of bias and completion of output
  res <- data.frame(
    N = nmax,
    Mean = means,
    Lower.CI = lower,
    Upper.CI = upper
    )
  rownames(res) <- NULL

  return(res)
}
