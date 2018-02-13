
## Average diversity estimates over all rarefaction iterations (for each sample)
phyloseq_mult_raref_div <- function(physeq, SampSize = min(sample_sums(physeq)), iter = 1000, 
  divindex = c("Observed", "Shannon"), parallel = FALSE, verbose = TRUE, ...){
  # divindex = alpha-diversity measures
  # ... additional arguments to rarefy_even_depth

  ## Progress indicator
  if(verbose == TRUE){
    progr <- "text"
  } else {
    progr <- "none"
  }

  ## Multiple rarefaction
  if(verbose == TRUE){ cat("..Multiple rarefaction\n") }
  phys_raref <- phyloseq_mult_raref(physeq, SampSize = SampSize, iter = iter, multithread = parallel, ...)

  ## Estimate richness and diversity
  phys_diversity <- function(x, indices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")){
    res <- estimate_richness(x, split = T, measures = indices)
    res <- dfRowName(res, name = "Sample")
    return(res)
  }

  if(verbose == TRUE){ cat("..Diversity esimation\n") }
  DIV <- ldply(
    .data = phys_raref,
    .fun = phys_diversity,
    .id = "Iteration",
    indices = divindex,
    .progress = progr)

  ## bug in ldply? column name is "X1" instead of "Iteration"
  colnames(DIV)[1] <- "Iteration"

  ## Reshape data to long-format
  DIV.resh <- melt(
    data = DIV,
    id.vars = c("Iteration", "Sample"),
    variable.name = "Index",
    value.name = "Value")

  ## Average diversity estimates over all iterations
  if(verbose == TRUE){ cat("..Averaging diversity estimates\n") }
  DIV.means <- ddply(
    .data=DIV.resh,
    .variables=c("Sample", "Index"),
    .fun=function(x) { mean_cl_boot(x$Value) })  # mean_cl_boot = ggplot2-wrapper for Hmisc

  colnames(DIV.means) <- c("Sample", "Index", "Estimate", "CI.lower", "CI.upper")

  return(DIV.means)
}
