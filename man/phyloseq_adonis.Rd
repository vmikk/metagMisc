% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq_adonis.R
\name{phyloseq_adonis}
\alias{phyloseq_adonis}
\title{Permutational multivariate analysis of variance using distance matrices (PERMANOVA) for phyloseq objects}
\usage{
phyloseq_adonis(
  physeq,
  group.var = NULL,
  method = "bray",
  pairwise = FALSE,
  add_permdisp = TRUE,
  permut = 999,
  permdisp_type = "median",
  p.adj = "fdr",
  all_results = TRUE,
  comparison_sep = ".",
  ...
)
}
\arguments{
\item{physeq}{A phyloseq-class object}

\item{group.var}{Variable name (contained in \code{\link{sample_data}}) to test (RHS in adonis formula)}

\item{method}{Distance/dissimilarity method (as character string; see \code{\link{distanceMethodList}})}

\item{pairwise}{Logical, perform pairwise adonis (dafault, FALSE)}

\item{add_permdisp}{Logical; if TRUE (default), results of tests for homogeneity of multivariate dispersions will be added to output (see \code{\link[vegan]{betadisper}})}

\item{permut}{Number of permutations required}

\item{permdisp_type}{Use the spatial median (default) or the group centroid for the analysis for homogeneity of multivariate dispersions (see \code{\link[vegan]{betadisper}})}

\item{p.adj}{Logical or character; if TRUE, adjust P-values for multiple comparisons with FDR; if character, specify correction method from \code{\link{p.adjust}} ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none")}

\item{all_results}{Logical, return results of adonis and data subsets for each pairwise comparison}

\item{comparison_sep}{Character string to separate the levels of independent variable the in the pairwise comparison names (default, ".")}

\item{...}{Additional arguments will be passed to \code{\link[vegan]{adonis}} and \code{\link[vegan]{permutest}}}
}
\value{
List with adonis and permdisp results
}
\description{
Permutational multivariate analysis of variance using distance matrices (PERMANOVA) for phyloseq objects
}
\details{
Currently, only categorical variables (factors) are supported in `group.var`.
}
\examples{
data("GlobalPatterns")
phyloseq_adonis(physeq = GlobalPatterns, group.var = "SampleType", all_results = FALSE, pairwise = TRUE, permut = 20)

}
\seealso{
\code{\link[vegan]{adonis}}, \code{\link[vegan]{betadisper}}, \code{\link[metagMisc]{adonis_pairwise}}
}
