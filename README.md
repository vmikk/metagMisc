# metagMisc
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.571403.svg)](https://doi.org/10.5281/zenodo.571403)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/metagMisc)](https://cran.r-project.org/package=metagMisc)
[![Latest Github release](https://img.shields.io/github/release/vmikk/metagMisc.svg)](https://github.com/vmikk/metagMisc/releases/latest)

Miscellaneous functions for metagenomic analysis.

The repository is currently in **ALPHA** state. Nothing is guaranteed and the material is subject to change without a notice (e.g., function names or arguments).

# Getting started

Vignette is under construction.

# Package features
- Multiple rarefaction
- OTU abundance averaging following CoDa (Compositional Data Analysis) workflow
- Phylogenetic diversity estimation
- Pairwise dissimilarity boxplots
![alt text](vignettes/Pairwise_dissimilarity_boxplot.png "Pairwise dissimilarity boxplots (enterotype data)")
- Prevalence plots (total OTU abundance vs OTU prevalence)
![alt text](vignettes/Prevalence_plots.png "Prevalence plots (GlobalPatterns data)")
- Diversity profiles based on Hill numbers (with `entropart` package)
![alt text](vignettes/Diversity_profile.png "Diversity profiles (esophagus data)")
- Extraction of the most abundant OTUs
![alt text](vignettes/Filter_top_taxa.png "Top taxa (GlobalPatterns data)")


# Installation
```
devtools::install_github("vmikk/metagMisc")
```

## Dependencies

`source("http://bioconductor.org/biocLite.R")`
* phyloseq: `biocLite("phyloseq")`
* dada2: `biocLite("dada2")`
* ALDEx2: `biocLite("ALDEx2")`
* metagenomeSeq: `biocLite("metagenomeSeq")`
* DESeq2: `biocLite("DESeq2")`
* vegan: `install.packages("vegan")`
* ggplot2
* plyr
* openssl

# Acknowledgements
The development of this software was supported by RFBR grant 16-04-01259.
