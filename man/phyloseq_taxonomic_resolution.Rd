% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq_taxonomic_resolution.R
\name{phyloseq_taxonomic_resolution}
\alias{phyloseq_taxonomic_resolution}
\title{Estimate taxonomic resolution of data}
\usage{
phyloseq_taxonomic_resolution(physeq, add_counts = TRUE, justDF = FALSE)
}
\arguments{
\item{physeq}{A phyloseq-class object}

\item{add_counts}{Logical, add read and OTU counts on top of the bars (default, TRUE)}

\item{justDF}{Logical, just return the data.frame, do not plot (default, FALSE)}
}
\value{
ggplot-object (if justDF = FALSE) or a data.frame (if justDF = TRUE).
}
\description{
Estimate taxonomic resolution of data
}
\details{
This function counts the number of OTUs and reads annotated 
at the lowest level of taxonomic classification and shows it on a barplot.
}
\examples{
data(GlobalPatterns)
gpt <- phyloseq_taxonomic_resolution(GlobalPatterns)
gpt
gpt$data_wide

}
\seealso{
\code{\link{get_max_taxonomic_rank}}
}
