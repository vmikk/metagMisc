
## RTK-based rarefaction function
raref_rtk <- function(physeq, SampSize = NULL, MinSizeTreshold = NULL, iter = 1000, trimOTUs = TRUE, ...){

  ## Filter samples by number of reads
  if(!is.null(MinSizeTreshold)){ x <- phyloseq::prune_samples(phyloseq::sample_sums(x) >= MinSizeTreshold, x) }

  ## Define rarefication depth
  if(is.null(SampSize)){ SampSize <- round( 0.9*min(phyloseq::sample_sums(x)) ) }

  ## Extract OTU abundance table (should be a matrix)
  abunds <- as( otu_table(physeq), "matrix" )

  ## Define margins
  TR <- taxa_are_rows(physeq)
  if(TR == TRUE) { MAR <- 2 }
  if(TR == FALSE){ MAR <- 1 }
  #   columns represent single samples (margin=2)
  #   rows are assumed to be samples (margin=1)

  ## Perform rarefaction
  rar <- rtk_mod(abunds, repeats = iter, depth = SampSize,
    margin = MAR, verbose = FALSE, threads = 1, ...)

  ## Function to replace OTU table in phyloseq object
  subst_otu_tab <- function(phys, newotu, taxrows, drop_zeros = FALSE){
    otu_table(phys) <- otu_table(newotu, taxa_are_rows = taxrows)

    ## Remove OTUs from the dataset that are no longer observed in any sample 
    if(drop_zeros == TRUE){
      phys <- prune_taxa(taxa_sums(phys) > 0, phys)
    }

    return(phys)
  }
  ## e.g. subst_otu_tab(physeq, newotu = rar[[1]])

  ## Re-create phyloseq object (single depth case)
  if(length(SampSize) == 1){
    res <- llply(
      .data = rar,
      .fun = function(z){ 
        subst_otu_tab(phys = physeq, newotu = z, taxrows = TR, drop_zeros = trimOTUs)
      })
  }

  ## Re-create phyloseq object (multiple depth case)
  # if(length(SampSize) > 1){
  #   ///////////////////       - TO DO
  # }


  ## Add rarefaction parameters as attributes to the phyloseq object
  attr(res, which = "RarefactionDepth") <- SampSize
  attr(res, which = "taxa_are_rows") <- TR

  return(res)
}
