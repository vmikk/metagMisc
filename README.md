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
```

## Dependencies

`source("http://bioconductor.org/biocLite.R")`
* phyloseq: `biocLite("phyloseq")`
* dada2: `biocLite("dada2")`
* ALDEx2: `biocLite("ALDEx2")`
* metagenomeSeq: `biocLite("metagenomeSeq")`
* DESeq2: `biocLite("DESeq2")`
* vegan: `install.packages("vegan")`
* PhyloMeasures: `remotes::install_github("cran/PhyloMeasures")`
* speedyseq: `remotes::install_github("mikemc/speedyseq")`
* ggplot2
* plyr
* openssl

# Acknowledgements
`metagMisc` stands on the shoulders of numerous R-packages (see Dependencies). In particular, it would not have happened without [phyloseq](https://github.com/joey711/phyloseq/) and [vegan](https://github.com/vegandevs/vegan/) packages. Please cite R and R packages when you use them for data analysis. 

The development of this software was supported by RFBR grants 16-04-01259 and 15-29-02765.
