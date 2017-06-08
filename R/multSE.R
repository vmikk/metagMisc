#######################################################################################
#
# Function for calculating multivariate standard error as the
# residual mean square error from a PERMANOVA for a one-way model, including double resampling.
# (Note also that sample sizes may vary among groups).

# The input will be a distance matrix among all samples, the grouping vector, and the number of
# re-samples (nresamp).
#
# The routine then calculates the means using the permutation approach, while the
# lower and upper quantiles are obtained using the bootstrapping approach
# including an adjustment for the bias in the bootstrap.
#
# The output consists of three vectors in a matrix called "output":
# 1. "means" = the means for each sample size
# 2. "upper" the upper 0.975 quantile for each sample size
# 3. "lower" the lower 0.0.025 quantile for each sample size
#
#######################################################################################


MSE.d = function (D, group, nresamp = 1000, ...) {

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
  ng <- length(levels(group))
  n.i <- table(group)
  nmax <- min(n.i)
  N <- sum(n.i)
  index <- 1:N

  # Ensure distance matrix is in the form of a matrix (rather than a "distance" object)
  D <- as.matrix(D)

  # Setting up the vectors for the results
  # Note that these do not have to be separate whole matrices now.
  means <- means.b <- lower <- upper <- rep(0,nmax)

     # One matrix is used to store the values under permutation resampling for each sample size
     multSE.store.p <- matrix(rep(0,nresamp*nmax), ncol = nmax, nrow = nresamp)
     # One matrix is used to store the values under bootstrap resampling for each sample size
     multSE.store.b <- matrix(rep(0,nresamp*nmax), ncol = nmax, nrow = nresamp)

        # Resampling loop for each sample size.
        for (nsub in 2:nmax) {
           for (iresamp in 1:nresamp) {
              ivec.p <- sample(index[group == levels(group)[1]], size = nsub, replace = FALSE)
              ivec.b <- sample(index[group == levels(group)[1]], size = nsub, replace = TRUE)
              for (i in 2:ng) {
                 ivec.p <- c(ivec.p, sample(index[group == levels(group)[i]], size = nsub, replace = FALSE))
                 ivec.b <- c(ivec.b, sample(index[group == levels(group)[i]], size = nsub, replace = TRUE))
              }
              group.resamp <- factor(rep(1:ng, each = nsub))
              D.perm <- D[ivec.p, ivec.p]
              D.boot <- D[ivec.b, ivec.b]
              multSE.store.p[iresamp, nsub] <- sqrt(MSE(D.perm, group.resamp)/nsub)
              multSE.store.b[iresamp, nsub] <- sqrt(MSE(D.boot, group.resamp)/nsub)
           }
        }

  # Means and quantiles
  means <- colMeans(multSE.store.p)
  means.b <- colMeans(multSE.store.b)
  upper <- apply(multSE.store.b, MARGIN = 2, quant.upper)
  lower <- apply(multSE.store.b, MARGIN = 2, quant.lower)

  # Calculation of bias and completion of output
  bias <-  means - means.b
  lower <- lower + bias
  upper <- upper + bias
  output <- cbind(1:nmax, means, lower, upper)
  colnames(output) <- c("n","mean", "lower.025", "upper.975")
  return(as.data.frame(output))
}
