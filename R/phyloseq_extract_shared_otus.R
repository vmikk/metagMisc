
## Extract species in common between two samples
phyloseq_extract_shared_otus <- function(x, samp1, samp2){
	# x = phyloseq object
	# samp1 and samp2 = sample names (character string)

	require(phyloseq)

	# test if the sample names are valid
	if( any(!c(samp1, samp2) %in% sample_names(x)) ){
		stop("Check the sample names, not all of them are present in the phyloseq object.\n")
	}

	# extract samples
	xx <- prune_samples(samples = c(samp1, samp2), x = x)

	# subset to OTUs that are present in both samples
	xx <- filter_taxa(xx, function(z){ sum(z >= 1) == 2 }, TRUE)

	return(xx)
}


## example
phyloseq_extract_shared_otus(esophagus, samp1 = "B", samp2 = "C")
