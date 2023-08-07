
# metagMisc v.0.5.0 (git xxxx; Mar xx, 2023)

- New functions
   * make rarefaction based on [RTK](https://github.com/hildebra/Rarefaction)
   * `phyloseq_num_shared_otus.R` - estimation of the number of shared and non-shared OTUs between all pairwise combinations of samples
   * `phyloseq_bootstrap` - bootstraps samples in phyloseq object
   * `phyloseq_combine_samples` - sums OTU abundances of replicates or sample groups into a single sample
   * `phyloseq_coverage` - estimates the observed abundance-based sample coverage
   * `coverage_to_samplesize` - estimates the required sample size for a particular coverage
   * `phyloseq_coverage_raref` - coverage-based rarefaction
   * `phyloseq_effort_div` - computes OTU diversity for a particular sample size or coverage
   * `phyloseq_effort_div_rangeplot` - vizualization of OTU diversity and richness for a particular sequencing depth or sample coverage ranges
   * `phyloseq_otu_appearance` - determines which OTUs appeared or disappeared in comparison with reference sample group
   * `phyloseq_richness_filter` - removes samples from phyloseq object that have less than `n` taxa
   * `phyloseq_inext` - estimates interpolated and extrapolated Hill numbers and sample coverage and constructs rarefaction curve
   * `phyloseq_merge_samples` - merges samples by name
   * `phyloseq_mult_raref_dist` - performs multiple rarefactions and averages sample dissimilarity across rarefactions
   * `phyloseq_rename_with_tax` - renames phyloseq OTUs or species with taxonomy ranks
   * `phyloseq_sep_pairwise` - splits phyloseq-class object for pairwise comparisons by sample-level variable
   * `phyloseq_taxonomic_resolution` - summarizes taxonomic resolution of the data
   * `phyloseq_taxonomy_imputation` - replaces missing taxonomy with higher ranks
   * `phyloseq_replace_zero` - replaces zeros in OTU abundance with psudocount or min observed abundnce
   * `phyloseq_transform_aldex_clr` - performs ALDEx2-based centred log-ratio transformation of OTU table  
   * `get_max_taxonomic_rank_DT` - determines the lowest level of taxonomic classification (with `data.table`)  
   * `prepare_breakaway` - prepares frequency counts for `breakaway` package  
   * `raref_rtk` - `RTK`-based rarefaction
   * `read_m8` - reads the standard BLAST m8 format (`-outfmt 6`)
   * `blast_to_wide` - converts BLAST results to a wide table
   * `SES` - standardized effect size estimation
   * `centroid_coords` - coordinates of centroids and medoids
   * `chunk` - split vector into chunks of equal size
   * `chunk_table` - split data.table into chunks
   * `count_primers` - counts number of reads in which the primer is found
   * `jaccard_bray` - converts between Bray-Curtis and Jaccard dissimilarity indicies
   * `leading_zero` - adds or removes leading zeros
   * `zero_pad` - pads a number with leading zeros
   * `package_availability` - checks the availability of dependencies
   * `vegdist_zero` - Bray-Curtis index with double-zero support

- Bug fixes
   * `phyloseq_extract_shared_otus` and `phyloseq_extract_non_shared_otus` - issue warning in the case if no OTUs were found
   * `get_max_taxonomic_rank` - data types bug
   * `phyloseq_average` - progress bar fixed
   * `parse_uc` - fixed a column name problem  
   * `prevalence` was added to NAMESPACE (thanks @LukDrey for reporting the issue)  

- The other updates  
   * `phyloseq_to_df` - handles taxa as rows case and checks the presence of taxonomy in phyloseq object  
   * `parse_uc` - is much faster now with `data.table`  
   * `parse_uc` - new option to remove duplicated rows  
   * `parse_taxonomy_amptk` rewritten  
   * `parse_taxonomy_amptk_batch` - supports multithreading now  
   * `phyloseq_standardize_otu_abundance` - added `wisconsin` and `rhea` methods  
   * `phyloseq_filter_prevalence` - informative error messages added  
   * `phyloseq_sep_variable` - informative error messages added  
   * `phyloseq_summary` - add OTU occurrence stats  
   * `eig_perc` - added plots of explained variance percentages  
   * `adonis_pairwise`  
   * `phyloseq_ntaxa_by_tax`  
   * `some` - works with phyloseq now  
   * `taxonomy_to_qiime` - supports Domain rank and performs validation of user-supplied ranks  
   * `phyloseq_otu_occurrence` - added counts  
   * `phyloseq_filter_sample_wise_abund_trim` - added relative abundance option  
   * `phyloseq_phylo_ses` - added VPD metric  
   * `phyloseq_group_dissimilarity` - added comparison types (between/within group)  
   * `phyloseq_extract_shared_otus` - estimates percentage of shared OTUs  
   * `prevalence` is now much faster with `data.table`  

- New logo, designed by Olesya Dulya  
- Documentation updates  

# metagMisc v.0.0.4 (git b0802ea; Feb 13, 2018)

- New functions
   * `phyloseq_mult_raref_div` - for averaging of diversity metrics across multiple rarefaction iterations
   * `phyloseq_to_MetaCommunity` - converter of phyloseq class to MetaCommunity object from entropart package
   * `phyloseq_otu_occurrence` - estimates taxa occurrence or occurrence frequency
   * `phyloseq_standardize_otu_abundance` - wrapper of decostand function from vegan package
   * `phyloseq_filter_top_taxa` - extracts the most abundant taxa from phyloseq-objects
   * `phyloseq_add_max_tax_rank` - adds the lowest level of taxonomic classification to the taxonomy table
   * `taxonomy_to_qiime` - prepares taxonomy in QIIME-style

- New experimental functions
   * `phyloseq_average`
   * `phyloseq_mult_raref_avg`
   * This functions implements OTU abundance averaging following compositional data analysis workflow.

- Bug fixes
   * `phyloseq_mult_raref` - multithreading fix
   * `phyloseq_to_df` - fix bug with taxonomy order
   * `parse_taxonomy_amptk_batch` - fix column reordering
   * `phyloseq_transform_vst_blind` - add workaround for excessive zeros

- Documentation and some minor updates


# metagMisc v.0.0.3 (git 924976f; Nov 8, 2017)

- New functions
   * `phyloseq_phylo_div` - for phylogenetic diversity estimation
   * `phyloseq_phylo_ses` - standardized effect size for phylogenetic diversity
   * `phyloseq_randomize` - randomization of abundance table and phylogeny in phyloseq objects (for null model testing and simulation)
   * `mult_dissim` - compute beta diversity for each rarefaction iteration
   * `mult_dist_average` - average multiple distance matrices
   * `dissimilarity_to_distance` - transform (non-metric) dissimilarity matrix to a weighted Euclidean distance (metric), as in Greenacre 2017
   * `add_metadata` - wrapper to add sample metadata to data frames
   * `uc_to_sumaclust` - convert UC files (from USEARCH or VSEARCH) to Sumaclust OTU observation maps

- Enhancements
   * `phyloseq_mult_raref` is able to run in parallel now; custom seeds option and rarefaction attributes were added, fixed auto-detection of rarefaction depth
   * `phyloseq_compare` was completely rewritten to allow for one or multiple input phyloseq objects; additional abundance statistics were added; fixed bug with percentage estimation
   * `prevalence` - mean and median OTU abundance were added

- Documentation updates and some other improvements


# metagMisc v.0.0.2 (git fd2697b; Oct 4, 2017)

- New functions
   * `prepare_inext`
   * `phyloseq_filter_sample_wise_abund_trim`
   * `make_utax_taxonomy`
   * `make_utax_taxonomy_batch`
   * `check_tax_uniqueness`

- Enhancements:
   * `phyloseq_group_dissimilarity` - generalized to more than 2 groups
   * `phyloseq_filter_prevalence` - added conditional option (AND/OR) for abundance and prevalence filtering

- Documentation updates


# metagMisc v.0.0.1 (git 591cb8a; Aug 23, 2017)

- New data-handling functions
   * `phyloseq_filter_prevalence`
   * `phyloseq_prevalence_plot`
   * `physeq_rm_na_tax`
   * `phyloseq_sep_variable`

- New functions for dissimilarity analysis
   * `phyloseq_group_dissimilarity`
   * `MultSE`

- Functional diversity related functions
   * `predict_metagenomes`
   * `metagenome_contributions`
   * `filter_cazy`

- New general purpose functions
   * `some`
   * `dist2list`

- Documentation updates and bug fixes


# metagMisc v.0.0.0.9000 (git 5735960; May 4, 2017)

- The first alpha release of metagMisc.
