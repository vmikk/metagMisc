
## Merge samples by name
## NB. sample metadata of the first ID from `samples_to_merge` will be used as a representative for a group
phyloseq_merge_samples <- function(phys, samples_to_merge, new_id = NULL){

  ## Check if more than one sample name was provided
  if(length(unique(samples_to_merge)) <= 1){
    warning("Nothing to merge. Please provide several sample IDs that should be merged.\n")
    return(phys)
  }

  ## Check if valid sample names were provided
  if(any(!samples_to_merge %in% sample_names(phys))){ 
    stop("Provided sample IDs are missing in the phyloseq object.\n")
  }

  ## Check new sample ID
  if(!is.null(new_id) & length(new_id) != 1){ 
    stop("Error: please provide a single value for a new sample name.\n")
  }

  ## Split phyloseq into samples that should be merged and the rest
  ps_m <- prune_samples(sample_names(phys) %in% samples_to_merge, phys)
  ps_r <- prune_samples(!sample_names(phys) %in% samples_to_merge, phys)

  ## Merge chosen samples
  ps_mrg <- phyloseq_combine_samples(ps_m)

  ## Assign new sample name
  if(is.null(new_id)){
    new_id <- paste(samples_to_merge, collapse = "__")
  }
  sample_names(ps_mrg) <- new_id

  ## Recover sample metadata
  if(!is.null(sample_data(phys, errorIfNULL=F))){
    ps_meta <- prune_samples(sample_names(phys) %in% samples_to_merge[1], phys)
    sample_names(ps_meta) <- new_id
    sample_data(ps_mrg) <- sample_data(ps_meta)
    rm(ps_meta)
  }

  ## Merge phyloseq objects back
  res <- merge_phyloseq(ps_mrg, ps_r)

  return(res)
}
