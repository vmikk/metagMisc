
#' @title Pairwise comparisons for permutational multivariate analysis of variance using distance matrices
#'
#' @param x Sample meta-data (data frame for the independent variables)
#' @param dd Dissimilarity matrix between samples
#' @param group.var Name of the independent variable to test (RHS in adonis formula)
#' @param permut Number of permutations required
#' @param p.adj Logical, adjust P-values for multiple comparisons
#' @param adj.meth Correction method from \code{\link{p.adjust}} ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
#' @param all_results Logical, return results of adonis and data subsets for each pairwise comparison
#' @param comparison_sep Character string to separate the levels of independent variable the in the pairwise comparison names (default, ".")
#' @param ... Additional arguments will be passed to \code{\link{adonis}}
#'
#' @return List with ......
#' @export
#' @seealso \code{\link{adonis}}
#' @examples
#' library(vegan)
#' data(dune)
#' data(dune.env)
#'
#' # Compare all Management levels
#' adonis(dune ~ Management, data = dune.env)
#'
#' # Pairwise comparisons between Management levels
#' ad <- adonis_pairwise(x = dune.env, dd = vegdist(dune), group.var = "Management")
#' ad$Adonis.tab
#'
adonis_pairwise <- function(x, dd, group.var, add_permdisp = TRUE, permut = 999,
  p.adj=T, adj.meth="fdr", all_results = T, comparison_sep = ".", permdisp_type = "median", ...){

  # require(vegan)

  ## Test if the selected predictor is in the meta-data
  if(!group.var %in% colnames(x)){
    stop("Check the 'group.var' value: independent variable '", group.var, "' is missing in the sample meta-data.\n")
  }

  VV <- which(colnames(x) %in% group.var)

  ## Combinations to test
  cc <- combn(levels(x[,VV]), 2)

  dd <- as.matrix(dd)

  ## Function to subset distance matrix
  prep.dist <- function(x, dd, sub1, sub2){
    subs <- which(x[,VV] == sub1 | x[,VV] == sub2)
    gr <- x[,VV][subs]
    res <- list()
      res$dd <- as.dist(dd[subs, subs])
      res$gr <- gr
    return(res)
  }

  ## Prepare data for pairwise comparisons
  dd.subs <- list()
  dd.groups <- list()
  for(i in 1:ncol(cc)){
    tmp <- prep.dist(x, dd, sub1=cc[1,i], sub2=cc[2,i])
    dd.subs[[i]] <- tmp$dd
    dd.groups[[i]] <- tmp$gr
    rm(tmp)
  }

  ## Pairwise tests
  adon <- list()              # placeholders
  if(add_permdisp == TRUE){ 
    permd <- list()   # for betadisper
    permdt <- list()  # for permutest
  }

  for(i in 1:length(dd.subs)){
    
    ## Multivariate analysis of variance (adonis)
    adon[[i]] <- vegan::adonis(dd.subs[[i]] ~ dd.groups[[i]], permutations = permut, ...)

    ## Multivariate homogeneity of groups dispersions (betadisper) and permutation test
    if(add_permdisp == TRUE){
      permd[[i]] <- vegan::betadisper(d = dd.subs[[i]], group = dd.groups[[i]], type = permdisp_type, ...)
      permdt[[i]] <- vegan::permutest(permd[[i]], permutations = permut, pairwise = FALSE)
    }
  }

  ## Prepare names of the pairwise comparisons
  foo <- vector()
  for(i in 1:ncol(cc)){
    foo[i] <- paste(cc[,i], collapse = comparison_sep)
  }
  names(dd.subs) <- names(dd.groups) <- names(adon) <- foo
  if(add_permdisp == TRUE){ names(permd) <- names(permdt) <- foo }

  ## Extract adonis results
  adon.extract.p <- function(adon){ adon$aov.tab$Pr[1] }
  adon.extract.F <- function(adon){ adon$aov.tab$F.Model[1] }
  adon.extract.df <- function(adon){ paste(adon$aov.tab$Df[1:2], collapse=";") }
  ad.t <- data.frame(Comparison = foo,
          F = unlist(lapply(adon, FUN = adon.extract.F)),
          df = unlist(lapply(adon, FUN = adon.extract.df)),
          p = unlist(lapply(adon, FUN = adon.extract.p)))
  rownames(ad.t) <- NULL

  ## Extract permutest results
  permutest_extract <- function(z){
    data.frame(F = z$tab$F[1], df = paste(z$tab$Df, collapse = ";"), p = z$tab$`Pr(>F)`[1])
  }
  pr.t <- plyr::ldply(.data = permdt, .fun = permutest_extract, .id = "Comparison")


  ## Adjust P-values
  if(p.adj == TRUE){
    ad.t$p.adj <- p.adjust(ad.t$p, method = adj.meth)

    if(add_permdisp == TRUE){
      pr.t$p.adj <- p.adjust(pr.t$p, method = adj.meth)
    }
  }

  ## Prepare the output
  res <- list()
  res$Adonis.tab <- ad.t
  res$Betadisper.tab <- pr.t

  ## Add additional data to the results
  if(all_results == TRUE){
    res$Adonis <- adon
    if(add_permdisp == TRUE){ res$Betadisper <- permd }
    res$Dist.subsets <- dd.subs
    res$Groups <- dd.groups
  }

  return(res)
}
