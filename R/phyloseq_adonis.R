## Adonis for phyloseq objects
phyloseq_adonis <- function(physeq, group.var = NULL, method = "bray",
  pairwise = FALSE, add_permdisp = TRUE, permut = 999, permdisp_type = "median",
  p.adj = TRUE, adj.meth = "fdr", all_results = TRUE, comparison_sep = ".", ...){

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
    adon <- vegan::adonis(as.formula( paste("dd ~ ", group.var, sep = "") ),
      data = metad, permutations = permut)

    ## Multivariate homogeneity of groups dispersions (betadisper) and permutation test
    if(add_permdisp == TRUE){
      permd <- vegan::betadisper(d = dd, group = metad[, group.var], type = permdisp_type)
      permdt <- vegan::permutest(permd, permutations = permut, pairwise = FALSE)
    }

    ## Extract adonis results
    ad.t <- data.frame(
      Factor = group.var,
      R2 = adon$aov.tab$R2[1], F = adon$aov.tab$F.Model[1],
      df = paste(adon$aov.tab$Df[1:2], collapse=";"), p = adon$aov.tab$Pr[1])

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
