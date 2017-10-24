
#' @title Estimate standardized effect sizes (SES) of phylogenetic diversity metrics (PD, MPD, MNTD, VPD) using randomization-based approach.
#' @description This function utilizes randomization-based approach to computet SES values
#' (in contrast to the exact standardisation approach used in \code{\link{phyloseq_phylo_div}} with 'standardize' option), 
#' which could be more computation intensive. Note that despite more freedom in defining the null-model with randomization-based approach 
#' it is debatable that this method can provide a reasonable approximation by selecting only a limited number of possible combinations 
#' (especially for large phylogenetic trees).
#'
#' Currently only non-abundance-weighted estimates are implemented.
#' @param physeq A phyloseq-class object (phylogenetic tree is required)
#' @param measures Character vector with diversity indices names ("PD", "MPD", "MNTD")
#' @param null_model Character string indicating which null model to use (for the supported list of models see \code{\link[picante]{ses.pd}})
#' @param nsim Number of randomizations for null-distribution generation
#' @param swapiter Number of iterations for independentswap or trialswap algorithms
#' @param verbose Logical; if TRUE, progress messages from the function will be printed
# @param ...
#' @details Values of SES above zero indicate that the species pool of a habitat is more diverse than the regional species pool.
#' @return Data frame.
#' @export
#' @seealso \code{\link{phyloseq_phylo_div}}, \code{\link{phyloseq_randomize}}, \code{\link[picante]{randomizeMatrix}}, \code{\link[picante]{ses.pd}}, \code{\link[picante]{ses.mpd}}, \code{\link[picante]{ses.mntd}}
#' @examples
#'
phyloseq_phylo_ses <- function(physeq, measures=c("PD", "MPD", "MNTD", "VPD"), 
  null_model, nsim = 1000, swapiter = 1000, verbose=TRUE, ...){

  require(plyr)
  require(picante)         # for null models
  require(PhyloMeasures)   # for phylogenetic diversity

  ## Data validation
  if( is.null(phy_tree(physeq, errorIfNULL=F)) ){ stop("Phylogenetic tree is missing in physeq.\n") }

  ## Print when analysis started
  if(verbose == TRUE){
    progr <- "text"
    cat("SES analysis started at ", format(Sys.time(), "%X"), "\n")
  } else {
    progr <- "none"
  }

  ## Check the orientation of the OTU table
  trows <- taxa_are_rows(physeq)

  ## Extact OTU table
  comm <- as(object = otu_table(physeq), Class = "matrix")

  ## Transpose OTU table (species should be columns for picante::randomizeMatrix)
  if(trows == TRUE){ comm <- t(comm) }

  ## Scale OTU abundance to presence/absence
  comm <- ifelse(comm > 0, 1, 0)

  ## Extact phylogenetic tree
  phy <- phy_tree(physeq)


  ## Function to estimate diversity with PhyloMeasures
  pdiv <- function(comm, phy, pdiv_measures){

    rez <- vector("list")  # initialize results

    if("PD" %in% pdiv_measures){
      rez <- c(rez, list(PD = pd.query(tree = phy, matrix = comm) ))
    }
    if("MPD" %in% pdiv_measures){
      rez <- c(rez, list(MPD = mpd.query(tree = phy, matrix = comm) ))
    }
    if("MNTD" %in% pdiv_measures){
      rez <- c(rez, list(MNTD = mntd.query(tree = phy, matrix = comm) ))
    }
    if("VPD" %in% pdiv_measures){
      rez <- c(rez, list(VPD = vpd(samp = comm, dis = phy) ))
    }

    rez <- do.call("cbind", rez)
    rez <- data.frame(SampleID = rownames(comm), rez)
    return(rez)
  }

  ## Observed diversity
  div.obs <- pdiv(comm, phy, pdiv_measures = measures)

  ## Function to randomize the observed  data
  null_mod_fun <- function(comm, phy, model, iter = swapiter){
    # model = models from picante::randomizeMatrix
    # iter = number of iterations for independentswap or trialswap algorithms

    if(model == "taxa.labels"){
      res_comm <- comm
      res_phy <- picante::tipShuffle(phy)
    }

    if(model %in% c("richness", "sample.pool")){  # same models?
      res_comm <- randomizeMatrix(comm, null.model="richness")
      res_phy <- phy
    }

    if(model == "frequency"){
      res_comm <- randomizeMatrix(comm, null.model="frequency")
      res_phy <- phy
    }

    if(model == "frequency"){
      res_comm <- randomizeMatrix(comm, null.model="frequency")
      res_phy <- phy
    }

    if(model == "phylogeny.pool"){
      res_comm <- randomizeMatrix(comm, null.model="richness")
      res_phy <- picante::tipShuffle(phy)
    }

    if(model == "independentswap"){
      res_comm <- randomizeMatrix(comm, null.model="independentswap", iterations = iter)
      res_phy <- phy
    }

    if(model == "trialswap"){
      res_comm <- randomizeMatrix(comm, null.model="trialswap", iterations = iter)
      res_phy <- phy
    }

    ## Return list with comm and phy
    rez <- list()
    rez$comm <- res_comm
    rez$phy <- res_phy
    return(rez)
  }

  ## Simulate multiple randomized communities
  ## TO DO ------- add parallel option here
  if(verbose == T){ cat("..Randomizing data with '", null_model, "' algorithm\n", sep = "") }
  nmods <- rlply(.n = nsim, .expr = null_mod_fun(comm, phy, model = null_model), .progress = progr)

  ## Esimate diversity metrics for each of the randomized communities
  ## TO DO ------- add parallel option here
  if(verbose == T){ cat("..Estimating phylogenetic diversity for the randomized data\n") }
  div.rnd <- ldply(
    .data = nmods,
    .fun = function(z, ...){ pdiv(comm = z$comm, phy = z$phy, ...) },
    pdiv_measures = measures,
    .progress = progr)

  ## Summarize null-distribution for each community (mean, SD and rank)
  if(verbose == T){ cat("..Estimating effect size\n") }
  rnd_mean <- ddply(.data=div.rnd, .variables="SampleID", .fun=numcolwise(mean, na.rm = T))
  rnd_sd <- ddply(.data=div.rnd, .variables="SampleID", .fun=numcolwise(sd, na.rm = T))
  rnd_rank <- ddply(.data = rbind(div.obs, div.rnd), .variables = "SampleID", numcolwise(.fun = function(z){ rank(z)[1] }))

  ## Rename columns
  colnames(rnd_mean)[-1] <- paste(colnames(rnd_mean)[-1], ".rand.mean", sep="")
  colnames(rnd_sd)[-1] <- paste(colnames(rnd_sd)[-1], ".rand.sd", sep="")
  colnames(rnd_rank)[-1] <- paste(colnames(rnd_rank)[-1], ".rank", sep="")

  ## Reorder observed diversity estimates to match the randomized order
  div.obs <- div.obs[match(x = div.obs$SampleID, table = rnd_mean$SampleID), ]

  ## Merge results
  res <- cbind(div.obs, rnd_mean[-1], rnd_sd[-1], rnd_rank[-1])

  ## Estimate Z-score and P-value
  if("PD" %in% measures){
    res$PD.z <- with(res, (PD - PD.rand.mean)/PD.rand.sd )
    res$PD.p <- with(res, PD.rank / (nsim + 1) )
  }
  if("MPD" %in% measures){
    res$MPD.z <- with(res, (MPD - MPD.rand.mean)/MPD.rand.sd )
    res$MPD.p <- with(res, MPD.rank / (nsim + 1) )
  }
  if("MNTD" %in% measures){
    res$MNTD.z <- with(res, (MNTD - MNTD.rand.mean)/MNTD.rand.sd )
    res$MNTD.p <- with(res, MNTD.rank / (nsim + 1) )
  }
  if("VPD" %in% measures){
    res$VPD.z <- with(res, (VPD - VPD.rand.mean)/VPD.rand.sd )
    res$VPD.p <- with(res, VPD.rank / (nsim + 1) )
  }

  ## Ending progress message
  if(verbose == TRUE){
    cat("..Done\n")
    cat("Analysis finished at ", format(Sys.time(), "%X"), "\n")
  }

  return(res)
}



## Internal function to estimate variance of pairwise distances separating taxa in a community
## Based on picante::mpd
vpd <- function(samp, dis, abundance.weighted=FALSE){
  # x = community data matrix (species as columns)
  # y = interspecific distance matrix (matrix or dist class) or phylogenetic tree of class phylo

  ## Convert phylogenetic tree to the pairwise distances
  if(class(dis) %in% "phylo"){ dis <- ape::cophenetic.phylo(dis) }
  if(class(dis) %in% "dist"){ dis <- as.matrix(dis) }

  N <- dim(samp)[1]
  res <- numeric(N)
  for(i in 1:N){

    ## Subset data
    sppInSample <- names(samp[i, samp[i, ] > 0])
    
    if(length(sppInSample) > 1){
      sample.dis <- dis[sppInSample, sppInSample]
      
      if(abundance.weighted == TRUE){
        sample.weights <- t(as.matrix(samp[i, sppInSample, drop=FALSE])) %*% as.matrix(samp[i, sppInSample, drop=FALSE])
        res[i] <- weighted.var(sample.dis, sample.weights)
      } else {
        res[i] <- var(sample.dis[lower.tri(sample.dis)])
      }
    }
    else{
      res[i] <- NA
    }
  }
  return(res)
}
# data(phylocom)
# vpd(phylocom$sample, phylocom$phylo, abundance.weighted=FALSE)
# vpd(phylocom$sample, phylocom$phylo, abundance.weighted=TRUE)


## Internal function to estimate weighted variance
weighted.var <- function(x, w = NULL, normwt = FALSE, na.rm = FALSE){
  # x = numeric vector
  # w = numeric vector of weights or NULL
  # normwt = Logical, to make weights sum to length(x) after deletion of NAs
  # na.rm = Logical, remove NAs

  ## If no weights are provied - just return ordinary variance
  if(is.null(w)){ return(var(x, na.rm = na.rm)) }

  ## Validate data
  if(length(x) != length(w)){ stop("Error: 'x' and 'w' must have the same length.\n")}

  ## Remove missing values
  if(na.rm == TRUE){
    nas <- is.na(x)
    if(any(nas)){
      w <- w[which(nas)]
      x <- x[which(nas)]
    }
  }
  
  ## Normalize weights
  if(normwt == TRUE){ w <- w * length(x) / sum(w) }

  sumw <- sum(w)
  wtm <- sum(w * x) / sumw    ## Weighted mean

  res <- sum(w *(x-wtm)^2)*(sumw / (sumw^2-sum(w^2)))
  # res <- (sum(w*x^2) * sumw - sum(w*x)^2) / (sumw^2 - sum(w^2))
  # res <- sum(w*((x - wtm)^2)) / (sumw - 1)

  return(res)
}
# x <- runif(500); wts <- sample(1:6, 500, TRUE)
# weighted.var(x, wts)

