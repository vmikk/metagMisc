% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metagenome_contributions.R
\name{metagenome_contributions}
\alias{metagenome_contributions}
\title{Partition metagenome functional contributions according to function, OTU, and sample.}
\usage{
metagenome_contributions(
  otu_tab,
  func_tab,
  tax_tab = NULL,
  features = NULL,
  NSTI_present = TRUE,
  rel_abund = TRUE,
  remove_zero_contributions = TRUE
)
}
\arguments{
\item{otu_tab}{Data frame with OTU abundances (rows = OTUs, columns = Samples, first column = OTU names)}

\item{func_tab}{Data frame with precalculated function predictions on per OTU basis (rows = OTUs, columns = feature counts, first column = OTU names)}

\item{tax_tab}{(Optional) Data frame with OTU taxonomy (rows = OTUs, columns = taxonomy ranks, first column = OTU names). If provided, taxonomy ranks will be added to the resulting table for each OTU}

\item{features}{Character vector with function names; if provided, results will be limited to only the specified functions}

\item{NSTI_present}{Logical; idnicating weather NSTI values are present in the last column of func_tab}

\item{rel_abund}{Logical; if TRUE, OTU counts will be transformed to relative abundances}

\item{remove_zero_contributions}{Logical; if TRUE, OTUs with zero contribution will be removed from results}
}
\value{
Data frame
}
\description{
This function partitions the predicted contribution to the metagenomes from each organism in the given OTU table for each function (e.g., KO) in each sample.
}
\details{
This function is analogous to the `metagenome_contributions.py` from PICRUSt.
Each line in the results relates how much a single OTU (third column) contributes to a single KO (first column) within a single sample (second column).
The fifth column contains the actual relative abundance contributed by this OTU, and the other columns contain other information about the abundance of the OTU and the percentage contribution of this OTU.
The last columns provide the taxonomy information for the OTU.
}
\examples{
## Create dummy data
set.seed(111)
NSAMP=10    # number of samples
NSPEC=30    # number of species/OTUs
NGENES=20   # number of features

dummy_name <- function(len = 5){ paste(sample(letters, size = len, replace = T), collapse = "") }

# Table with gene counts per OTU
func_tab <- data.frame(
              OTU_ID = replicate(n = 100, expr = dummy_name()),
              matrix(data = sample(0:100, size = 100*NGENES, replace = T), nrow = 100),
              stringsAsFactors = F)
colnames(func_tab)[-1] <- paste("F", 1:(ncol(func_tab)-1), sep="")

# Table with OTU abundances
otu_tab <- data.frame(
              OTU = sample(func_tab$OTU_ID, size = NSPEC),
              matrix(data = sample(1:10000, size = NSPEC*NSAMP, replace = T), nrow = NSPEC),
              stringsAsFactors = F)
colnames(otu_tab)[-1] <- paste("Samp", 1:(ncol(otu_tab)-1), sep="")

# Table with OTU taxonomy annotation
tax_tab <- data.frame(OTU_ID = func_tab$OTU_ID,
              Kingdom = sort(sample(LETTERS[1:4], size = nrow(func_tab), replace = T)),
              Phylum = sort(sample(LETTERS[5:20], size = nrow(func_tab), replace = T)),
              stringsAsFactors = F)

metag <- metagenome_contributions(otu_tab, func_tab, tax_tab,
                     features = c("F1", "F2", "F3"), NSTI_present = F, rel_abund = FALSE,
                     remove_zero_contributions = TRUE)
head(metag)

}
\references{
https://picrust.github.io/picrust/scripts/metagenome_contributions.html
}
