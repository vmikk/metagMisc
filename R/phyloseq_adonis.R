
#' @title Permutational multivariate analysis of variance using distance matrices (PERMANOVA) for phyloseq objects
#'
#' @param physeq A phyloseq-class object
#' @param group.var Variable name (contained in \code{\link{sample_data}}) to test (RHS in adonis formula)
#' @param method Distance/dissimilarity method (as character string; see \code{\link{distanceMethodList}})
#' @param pairwise Logical, perform pairwise adonis (dafault, FALSE)
#' @param add_permdisp Logical; if TRUE (default), results of tests for homogeneity of multivariate dispersions will be added to output (see \code{\link[vegan]{betadisper}})
#' @param permut Number of permutations required
#' @param permdisp_type Use the spatial median (default) or the group centroid for the analysis for homogeneity of multivariate dispersions (see \code{\link[vegan]{betadisper}})
#' @param p.adj Logical or character; if TRUE, adjust P-values for multiple comparisons with FDR; if character, specify correction method from \code{\link{p.adjust}} ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none")
#' @param all_results Logical, return results of adonis and data subsets for each pairwise comparison
#' @param comparison_sep Character string to separate the levels of independent variable the in the pairwise comparison names (default, ".")
#' @param ... Additional arguments will be passed to \code{\link[vegan]{adonis}} and \code{\link[vegan]{permutest}}
#'
#' @details Currently, only categorical variables (factors) are supported in `group.var`.
#' @return List with adonis and permdisp results
#' @export
#' @seealso \code{\link[vegan]{adonis}}, \code{\link[vegan]{betadisper}}, \code{\link[metagMisc]{adonis_pairwise}}
#' @examples
#' data("GlobalPatterns")
#' phyloseq_adonis(physeq = GlobalPatterns, group.var = "SampleType", all_results = FALSE, pairwise = TRUE, permut = 20)
#'
phyloseq_adonis <- function(physeq, group.var = NULL, method = "bray",
  pairwise = FALSE, add_permdisp = TRUE, permut = 999, permdisp_type = "median",
  p.adj = "fdr", all_results = TRUE, comparison_sep = ".", ...){

  ## Currently, independent variables can be only for factors !

  ## Extract sample metadata
  metad <- as(sample_data(physeq), "data.frame")
  metad[, group.var] <- factor(metad[, group.var])

  ## Estimate pairwise dissimilarity
  dd <- phyloseq::distance(physeq, method = method, type = "samples", ...)

  ## Perform pair-wise comparisons
  if(pairwise == TRUE){
    res <- adonis_pairwise(x = metad, dd = dd, group.var = group.var,
      add_permdisp = add_permdisp, permut = permut, p.adj = p.adj, adj.meth = adj.meth,
      all_results = all_results, comparison_sep = comparison_sep, permdisp_type = permdisp_type)

  ## Or perform a single test
  } else {

    ## Multivariate analysis of variance (adonis)
    adon <- vegan::adonis2(as.formula( paste("dd ~ ", group.var, sep = "") ),
      data = metad, permutations = permut)

    ## Multivariate homogeneity of groups dispersions (betadisper) and permutation test
    if(add_permdisp == TRUE){
      permd <- vegan::betadisper(d = dd, group = metad[, group.var], type = permdisp_type)
      permdt <- vegan::permutest(permd, permutations = permut, pairwise = FALSE)
    }

    ## Extract adonis results
    ad.t <- data.frame(
      Factor = group.var,
      R2 = adon$R2[1],
      F  = adon$F[1],
      df = paste(adon$Df[1:2], collapse=";"),
      p  = adon$Pr[1])

    ## Extract permutest results
    if(add_permdisp == TRUE){
      pr.t <- data.frame(
        Factor = group.var,
        F = permdt$tab$F[1],
        df = paste(permdt$tab$Df, collapse = ";"), p = permdt$tab$`Pr(>F)`[1])
    }

    ## Prepare the output
    res <- list()
    res$Adonis.tab <- ad.t
    if(add_permdisp == TRUE){ res$Betadisper.tab <- pr.t }

    ## Add additional data to the results
    if(all_results == TRUE){
      res$Adonis <- adon
      if(add_permdisp == TRUE){ res$Betadisper <- permd }
      res$Dist <- dd
    }
  } # end of pairwise == FALSE

  return(res)
}
