

## Count number of OTUs by tax rank for each sample
phyloseq_ntaxa_by_tax <- function(x, TaxRank = "Phylum", relative = F, add_meta_data = T){

    require(phyloseq)
    require(plyr)

    ## Melt phyloseq data object into large data.frame
    mm <- psmelt(physeq)

    ## Count number of OTUs for each sample
    count_otus <- function(z, TaxRank = TaxRank, relative = relative){

      ## Remove zero-OTUs
      zero_otu <- z$Abundance > 0
      if(any(zero_otu)){ z <- z[zero_otu, ] }

      ## Count number of OTUs per taxonomic rank selected
      rez <- as.data.frame(table(z[, TaxRank]), stringsAsFactors = F)
      colnames(rez) <- c(TaxRank, "N.OTU")

      ## Transform to relative abundance
      if(relative == TRUE){
        rez$N.OTU <- with(rez, N.OTU / sum(N.OTU) )
      }

      return(rez)
    }

    ## Count number of OTUs for each sample
    res <- ddply(.data = mm, .variables = "Sample", .fun = count_otus, TaxRank = TaxRank, relative = relative)

    ## Add meta-data
    if(add_meta_data == TRUE){
      ## Extract meta-data
      metad <- data.frame(Sample = sample_names(x), sample_data(x))

      ## Extract column names
      # main_cols <- c("OTU", "Sample", "Abundance", rank_names(x))          # 'standard' columns
      # meta_cols <- colnames(mm)[ which(!colnames(mm) %in% main_cols) ]     # meta-data columns
      main_cols <- c("Sample")                                               # 'standard' columns
      meta_cols <- colnames(metad)[ which(!colnames(metad) %in% main_cols) ] # meta-data columns

      res <- cbind(res, metad[match(x = res$Sample, table = metad$Sample), meta_cols] )
    }

    rownames(res) <- NULL
    return(res)
}
