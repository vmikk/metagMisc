% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq_num_shared_otus.R
\name{phyloseq_num_shared_otus}
\alias{phyloseq_num_shared_otus}
\title{Estimate the number of shared and non-shared species (OTUs/ESvs) between all pairwise comparisons of samples.}
\usage{
phyloseq_num_shared_otus(physeq, use_Matrix = TRUE)
}
\arguments{
\item{physeq}{A phyloseq-class object}

\item{use_Matrix}{Logical, use sparse matrix class to perform the analysis (default, TRUE)}
}
\value{
The function's output is a list containing three matrices
(rows and columns correspond to samples):

\itemize{
 \item{\strong{"shared"}}{ - The number of shared OTUs between the samples}

 \item{\strong{"nonshared_total"}}{ - The number of non-shared OTUs between the samples.
  This matrix contains the total number of non-shared OTUs between the two samples.
  The matrix is symmetric, meaning that the values in the matrix are equal on both sides of the diagonal.}

 \item{\strong{"nonshared_asymmetric"}}{ - The number of non-shared OTUs between the samples.
   The matrix is asymmetric, with the upper and lower triangular parts of
   the matrix representing the number of unique OTUs present only in one of the two samples.}
}
}
\description{
This function counts the number of OTUs that are common (shared) and
distinct (non-shared) between a pair of samples. All pairwise sample
combinations are considered, and the output is stored as a matrix.
}
\details{
This function uses the sparse matrix class implemented in the package
\code{\link[Matrix]{Matrix}} to improve computational speed and
decrease memory footprint.
}
\examples{
data(esophagus)

# Count number of shared and non-shared OTUs
ps <- phyloseq_num_shared_otus(esophagus)
ps

# The same, but using base-R matrix class (not recommended for large data)
ps <- phyloseq_num_shared_otus(esophagus, use_Matrix = FALSE)
ps

# To obtain a list of shared or non-shared OTUs between a pair of samples, use
phyloseq_extract_shared_otus(x = esophagus, samp_names = c("B", "C"))
phyloseq_extract_non_shared_otus(x = esophagus, samp_names = c("B", "C"))

}
\seealso{
\code{\link[metagMisc]{phyloseq_extract_shared_otus}},
\code{\link[metagMisc]{phyloseq_extract_non_shared_otus}}
}
