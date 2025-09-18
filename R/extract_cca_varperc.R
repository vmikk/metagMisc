
## Function to extract percentage of explained variance from RDA, CCA, and capscale
extract_cca_varperc <- function(x){
  # x <- attr(x = ORDS_Abund$AB_synth, which = "Ordination")

  chi <- c(x$tot.chi, x$pCCA$tot.chi, x$CCA$tot.chi, x$CA$tot.chi)
  props <- chi/chi[1]
  rnk <- c(NA, x$pCCA$rank, x$CCA$rank, x$CA$rank)
  if(inherits(x, "dbrda") &&
      (!is.null(x$CCA) && x$CCA$poseig < x$CCA$qrank ||
       !is.null(x$CA) && x$CA$poseig < x$CA$rank)){
    poseig <- c(NA, if (!is.null(x$pCCA)) NA, x$CCA$poseig, x$CA$poseig)
  } else {
    poseig <- NULL
  }
  tbl <- cbind(chi, props, rnk, poseig)
  if(!is.null(x$CA$imaginary.chi)){
    tbl <- rbind(tbl, c(x$CA$imaginary.chi, x$CA$imaginary.chi/x$tot.chi, x$CA$imaginary.rank, NULL))
  }
  colnames(tbl) <- c("Inertia", "Proportion", "Rank", if (!is.null(poseig)) "RealDims")
  rn <- c("Total", "Conditional", "Constrained", "Unconstrained", "Imaginary")
  rownames(tbl) <- rn[ c(TRUE, !is.null(x$pCCA), !is.null(x$CCA), !is.null(x$CA), !is.null(x$CA$imaginary.chi)) ]
  if(is.null(x$CCA) && is.null(x$pCCA)){ tbl <- tbl[, -2] }

  tbl <- metagMisc::dfRowName(tbl, name = "VarianceType")
  tbl <- as.data.table(tbl)
  return(tbl)
}
