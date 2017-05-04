# metagMisc
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.571403.svg)](https://doi.org/10.5281/zenodo.571403)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/metagMisc)](https://cran.r-project.org/package=metagMisc)
Miscellaneous functions for metagenomic analysis.

The repository is currently in **ALPHA** state. Nothing is guaranteed.

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
