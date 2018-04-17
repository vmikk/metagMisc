
# Estimate the observed abundance-based sample coverage for phyloseq object
phyloseq_coverage <- function(physeq, correct_singletons = FALSE){
  
  ## Prepare a list of OTU abundance vectors
  x <- prepare_inext(
        as.data.frame(otu_table(physeq, taxa_are_rows = T)),
        correct_singletons = correct_singletons)

  ## Estimate sample coverages
  res <- ldply(.data = x, .fun = function(z){ iNEXT:::Chat.Ind(z, sum(z)) })
  colnames(res) <- c("SampleID", "SampleCoverage")

  return(res)
}



## Estimate the required sample size for a particular coverage
## the code is based on iNEXT:::invChat.Ind by Johnson Hsieh (d76e3b8, Nov 12, 2016)
# https://github.com/JohnsonHsieh/iNEXT/blob/de46aeacb4433c539b2880df376e87b44bc1723c/R/invChat.R#L1
coverage_to_samplesize <- function(x, coverage = 0.95, add_attr = F){

  # iNEXT:::invChat.Ind(x, C = coverage)$m

  ## Total number of reads and the observed sample coverage
  n <- sum(x)
  refC <- iNEXT:::Chat.Ind(x, n)

  ## Interpolation
  f <- function(m, C) abs(iNEXT:::Chat.Ind(x, m) - C)
  if(refC > coverage){
    opt <- optimize(f, C = coverage, lower = 0, upper = sum(x))
    mm <- opt$minimum
    mm <- round(mm)
  }

  ## Extrapolation
  if(refC <= coverage){
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    if(f1>0 & f2>0)  {A <- (n-1)*f1/((n-1)*f1+2*f2)}
    if(f1>1 & f2==0) {A <- (n-1)*(f1-1)/((n-1)*(f1-1)+2)}
    if(f1==1 & f2==0){A <- 1}
    if(f1==0 & f2==0){A <- 1}
    mm <- (log(n/f1)+log(1-coverage))/log(A)-1
    mm <- n + mm
    mm <- round(mm)
  }

  ## Add attributes
  if(add_attr){
    if(refC > coverage) { attr(mm, "method") <- "interpolated" }
    if(refC <= coverage){ attr(mm, "method") <- "extrapolated" }
    attr(mm, "ObservedCoverage") <- refC
    attr(mm, "RequestedCoverage") <- coverage
  }
  
  return(mm)
}
# Example:  coverage_to_samplesize(x, coverage = 0.9, add_attr = T)

