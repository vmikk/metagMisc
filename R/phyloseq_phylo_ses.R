
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
#' # Load data
#' data(esophagus)
#'
#' # Estimate phylogenetic diversity index (PDI)
#' phyloseq_phylo_ses(esophagus, measures = "PD", null_model = "taxa.labels", nsim = 100, verbose = F)
#'
#' # Estimate phylogenetic diversity, mean pairwise distance,
#' # mean nearest taxon distance and variance of phylogenetic distances;
#' # use independent swap algorithm to generate the null distribution
#' phyloseq_phylo_ses(esophagus,
#'                    measures = c("PD", "MPD", "MNTD", "VPD"),
#'                    null_model = "independentswap",
#'                    nsim = 200, swapiter = 200)  # NB. increase the number of iterations!
#'
phyloseq_phylo_ses <- function(physeq, measures = c("PD", "MPD", "MNTD", "VPD"),
  null_model, package = "picante", abundance_weighted = FALSE, nsim = 1000, swapiter = 1000, verbose = TRUE, ...){

  # require(plyr)
  # require(picante)         # for null models
  # require(PhyloMeasures)   # for phylogenetic diversity

  ## Data validation
  if( is.null(phyloseq::phy_tree(physeq, errorIfNULL=F)) ){
    stop("Phylogenetic tree is missing in physeq.\n")
  }

  ## Print when analysis started
  if(verbose == TRUE){
    progr <- "text"
    cat("SES analysis started at ", format(Sys.time(), "%X"), "\n")
  } else {
    progr <- "none"
  }

  ## Check the orientation of the OTU table
  trows <- phyloseq::taxa_are_rows(physeq)

  ## Extact OTU table
  comm <- as(object = phyloseq::otu_table(physeq), Class = "matrix")

  ## Transpose OTU table (species should be columns for picante::randomizeMatrix)
  if(trows == TRUE){ comm <- t(comm) }

  ## Scale OTU abundance to presence/absence
  if(abundance_weighted == FALSE){
    comm <- ifelse(comm > 0, 1, 0)
  }

  ## Extact phylogenetic tree
  phy <- phyloseq::phy_tree(physeq)


  ## Function to estimate diversity with PhyloMeasures or picante
  pdiv <- function(comm, phy, pdiv_measures, method = "picante", abund_wei){

    ## Get pairwise distances from a phylogenetic tree
    if(method == "picante"){ dis <- ape::cophenetic.phylo(phy) }

    ## Initialize results
    rez <- vector("list")

    ## Estimate diversity with PhyloMeasures
    if(method == "PhyloMeasures"){
      if("PD" %in% pdiv_measures){
        rez <- c(rez, list(PD = PhyloMeasures::pd.query(tree = phy, matrix = comm) ))
      }
      if("MPD" %in% pdiv_measures){
        rez <- c(rez, list(MPD = PhyloMeasures::mpd.query(tree = phy, matrix = comm) ))
      }
      if("MNTD" %in% pdiv_measures){
        rez <- c(rez, list(MNTD = PhyloMeasures::mntd.query(tree = phy, matrix = comm) ))
      }
    }

    ## Estimate diversity with picante
    if(method == "picante"){
      if("PD" %in% pdiv_measures){
        if(abund_wei == FALSE){
          rez <- c(rez, list(PD = picante::pd(samp = comm, tree = phy, include.root = F)$PD ))
        }
        if(abund_wei == TRUE){
          rez <- c(rez, list(PD = pd_weighted(samp = comm, tree = phy) ))
        }
      }
      if("MPD" %in% pdiv_measures){
        rez <- c(rez, list(MPD = picante::mpd(samp = comm, dis = dis, abundance.weighted = abund_wei) ))
      }
      if("MNTD" %in% pdiv_measures){
        rez <- c(rez, list(MNTD = picante::mntd(samp = comm, dis = dis, abundance.weighted = abund_wei) ))
      }
    }

    ## Estimate VPD with internal function
    if("VPD" %in% pdiv_measures){
      rez <- c(rez, list(VPD = vpd(samp = comm, dis = phy, abundance.weighted = abund_wei) ))
    }

    rez <- do.call("cbind", rez)
    rez <- data.frame(SampleID = rownames(comm), rez)
    return(rez)
  }

  ## Observed diversity
  div.obs <- pdiv(comm, phy, pdiv_measures = measures, method = package, abund_wei = abundance_weighted)

  ## Function to randomize the observed  data
  null_mod_fun <- function(comm, phy, model, iter = swapiter){
    # model = models from picante::randomizeMatrix
    # iter = number of iterations for independentswap or trialswap algorithms

    if(model == "taxa.labels"){
      res_comm <- comm
      res_phy <- picante::tipShuffle(phy)
    }

    if(model %in% c("richness", "sample.pool")){  # same models?
      res_comm <- picante::randomizeMatrix(comm, null.model="richness")
      res_phy <- phy
    }

    if(model == "frequency"){
      res_comm <- picante::randomizeMatrix(comm, null.model="frequency")
      res_phy <- phy
    }

    if(model == "frequency"){
      res_comm <- picante::randomizeMatrix(comm, null.model="frequency")
      res_phy <- phy
    }

    if(model == "phylogeny.pool"){
      res_comm <- picante::randomizeMatrix(comm, null.model="richness")
      res_phy <- picante::tipShuffle(phy)
    }

    if(model == "independentswap"){
      res_comm <- picante::randomizeMatrix(comm, null.model="independentswap", iterations = iter)
      res_phy <- phy
    }

    if(model == "trialswap"){
      res_comm <- picante::randomizeMatrix(comm, null.model="trialswap", iterations = iter)
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
  nmods <- plyr::rlply(.n = nsim, .expr = null_mod_fun(comm, phy, model = null_model), .progress = progr)

  ## Esimate diversity metrics for each of the randomized communities
  ## TO DO ------- add parallel option here
  if(verbose == T){ cat("..Estimating phylogenetic diversity for the randomized data\n") }
  div.rnd <- plyr::ldply(
    .data = nmods,
    .fun = function(z, ...){ pdiv(comm = z$comm, phy = z$phy, ...) },
    pdiv_measures = measures,
    method = package, 
    abund_wei = abundance_weighted,
    .progress = progr)

  ## Summarize null-distribution for each community (mean, SD and rank)
  if(verbose == T){ cat("..Estimating effect size\n") }
  rnd_mean <- plyr::ddply(.data=div.rnd, .variables="SampleID", .fun=plyr::numcolwise(mean, na.rm = T))
  rnd_sd <- plyr::ddply(.data=div.rnd, .variables="SampleID", .fun=plyr::numcolwise(sd, na.rm = T))
  rnd_rank <- plyr::ddply(.data = rbind(div.obs, div.rnd), .variables = "SampleID", plyr::numcolwise(.fun = function(z){ rank(z)[1] }))

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




## Abundance weighted calculation of Faiths PD index
## based on lefse_0.5::weighted.faith by Nathan G. Swenson
## https://github.com/NGSwenson/lefse_0.5/commit/74814e2f4f8a1f0244da8e8a48bb17897718db67
pd_weighted <- function(samp, tree){

  ## Function for a single sample
  wpd <- function(smp){
    
    ## Extract the names of species in a community with an abundance greater than zero
    ## and make a pruned phylogeny for that community
    tmp.tree <- geiger::treedata(tree, smp[smp > 0], warnings=F)$phy

    ## Create empty branches matrix
    branches <- matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 4)

    ## Fill first two columns of the matrix with node numbers defining each edge
    branches[,1:2] <- tmp.tree$edge

    ## Fill the third column with the length of each branch
    branches[,3] <- tmp.tree$edge.length
    
    get.leaves <- function(x){ leaves.node <- geiger::tips(tmp.tree, x[2]) }
    
    ## Retrieve species names subtended by each branch (i.e. ## row) in the branches matrix
    leaves <- apply(X = branches, MARGIN = 1, FUN = get.leaves)
    
    ## Calculate the mean abundance (Ai) for species across each set of leaves
    for(i in 1:length(leaves)){
      branches[i, 4] <- mean(smp[leaves[[i]]], na.rm = T) 
    }

    ## Calculated the Weighted Faith’s Index
    rez <- nrow(tmp.tree$edge) * ((sum(branches[,3] * branches[,4])) / sum(branches[,4]))
    return(rez)
  }

  ## Estimate weithed PD for each sample
  res <- apply(X = samp, MARGIN = 1, FUN = wpd)
  res
}
