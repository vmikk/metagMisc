% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq_otu_occurrence.R
\name{phyloseq_otu_occurrence}
\alias{phyloseq_otu_occurrence}
\title{Estimate species occurrence (within groups of samples)}
\usage{
phyloseq_otu_occurrence(
  physeq,
  variable = NULL,
  taxa_frequency = "percentage",
  drop_zeroes = FALSE,
  justdf = FALSE,
  long = FALSE
)
}
\arguments{
\item{physeq}{A phyloseq-class object}

\item{variable}{Character string defining a variable name of the sample groups (this variable should be present in \code{\link[phyloseq]{sample_data}}) or NULL (no sample groups)}

\item{taxa_frequency}{Logical; if TRUE (default), relative frequency of species occurence within a sample group will be returned; if FALSE, number of samples within each sample group will be returned for each taxa}

\item{drop_zeroes}{Logical; if TRUE, taxa with total zero abundance will be removed}

\item{justdf}{Logical; if TRUE, return only a data frame with taxa occurrences; if FALSE (default), modified phyloseq object will be returned}

\item{long}{Logical; if TRUE, data frame with taxa occurrences will be returned in long format (with a single column defining the sample group); if FALSE (default), species occurrences will be returned in a wide format (sample groups as columns)}
}
\value{
A phyloseq-class object (if justdf = FALSE) or a data frame (justdf = TRUE).
}
\description{
Estimate species occurrence (within groups of samples)
}
\examples{
# Load data
data("GlobalPatterns")

## Return phyloseq-object with sample groups instead of samples
# With absolute counts (e.g., number of samples with the species)
gpa <- phyloseq_otu_occurrence(GlobalPatterns, variable = "SampleType", taxa_frequency = F)
gpa
head( otu_table(gpa) )

# With relative frequency of species occurence within a sample group
gpr <- phyloseq_otu_occurrence(GlobalPatterns, variable = "SampleType", taxa_frequency = T)
gpr
head( otu_table(gpr) )


## Return just a data frame with taxa occurrences
# In wide format (sample groups as columns)
gpw <- phyloseq_otu_occurrence(GlobalPatterns, variable = "SampleType", taxa_frequency = T, drop_zeroes = F, justdf = T, long = F)
head(gpw)

# In long format (with a single column with sample type)
gpl <- phyloseq_otu_occurrence(GlobalPatterns, variable = "SampleType", taxa_frequency = T, drop_zeroes = F, justdf = T, long = T)
head(gpl)

}
