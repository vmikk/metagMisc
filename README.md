# metagMisc <img src='man/figures/Logo_metagMisc.png' align="right" height="139" />

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597622.svg)](https://doi.org/10.5281/zenodo.597622)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/metagMisc)](https://cran.r-project.org/package=metagMisc)
[![Latest Github release](https://img.shields.io/github/release/vmikk/metagMisc.svg)](https://github.com/vmikk/metagMisc/releases/latest)

<!-- badges: end -->

## metagMisc tools for metabarcoding and microbiome data analysis

`metagMisc` is an R package that provides a suite of functions for metabarcoding and microbiome data analysis (also suitable for processed metagenomic data), with particular focus on handling data in `phyloseq` format ([McMurdie and Holmes 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217)). The package helps with common workflows including data filtering, normalization, rarefaction, diversity estimation, and comparative analysis. It offers both traditional and cutting-edge approaches for handling amplicon sequencing data, making it an essential toolkit for microbial ecologists and bioinformaticians.


# Package features
* Multiple rarefaction
* OTU abundance averaging following CoDa (Compositional Data Analysis) workflow
* Phylogenetic diversity estimation (including standardized effect sizes)
* Pairwise dissimilarity boxplots
<img src="man/figures/Pairwise_dissimilarity_boxplot.png" width="500" title="Pairwise dissimilarity boxplots (enterotype data)" />

### Visualization capabilities

**Pairwise dissimilarity analysis**  
<img src="man/figures/Pairwise_dissimilarity_boxplot.png" width="500" title="Pairwise dissimilarity boxplots" />

**Prevalence patterns**  
<img src="man/figures/Prevalence_plots.png" width="500" title="OTU prevalence vs abundance plots" />

**Diversity profiles with Hill numbers**  
<img src="man/figures/Diversity_profile.png" width="500" title="Diversity profiles based on Hill numbers" />

**Taxonomic filtering and abundance analysis**  
<img src="man/figures/Filter_top_taxa.png" width="500" title="Top taxa filtering visualization" />

**Taxonomic annotation quality assessment**  
<img src="man/figures/Taxonomic_resolution.png" width="500" title="Taxonomic resolution analysis" />

## Installation

Install the development version from GitHub:

```r
# Install from GitHub
if (!requireNamespace("remotes", quietly = TRUE)){
  install.packages("remotes")
}

remotes::install_github("vmikk/metagMisc")
```

## Dependencies

Most dependencies will be installed automatically. For Bioconductor packages, ensure you have BiocManager installed:

```r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

# Core Bioconductor dependencies
BiocManager::install(c("phyloseq", "DESeq2", "ALDEx2", "metagenomeSeq"))

# Additional dependencies (installed automatically with metagMisc)
install.packages(c("vegan", "data.table", "ggplot2", "plyr", "iNEXT", "SRS"))

# Optional packages for specific analyses
remotes::install_github("cran/PhyloMeasures")  # for phylogenetic diversity
remotes::install_github("mikemc/speedyseq")    # for faster phyloseq operations
BiocManager::install("GenomicAlignments")      # for CIGAR string expansion
```

## Acknowledgements

metagMisc builds upon the excellent work of the R community, particularly the [vegan](https://github.com/vegandevs/vegan/), [phyloseq](https://github.com/joey711/phyloseq/), and [data.table](https://github.com/Rdatatable/data.table) packages. We are grateful to all package developers who make reproducible microbiome research possible.

When using metagMisc, please also cite the underlying methods and packages used in your analysis. The package vignettes provide guidance on appropriate citations for specific functions.

The development of this software was supported by RFBR grants 16-04-01259 and 15-29-02765.
