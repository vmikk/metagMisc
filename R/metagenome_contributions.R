
#' @title Partition metagenome functional contributions according to function, OTU, and sample.
#' @description This function partitions the predicted contribution to the metagenomes from each organism in the given OTU table for each function (e.g., KO) in each sample.
#' @param otu_tab Data frame with OTU abundances (rows = OTUs, columns = Samples, first column = OTU names)
#' @param func_tab Data frame with precalculated function predictions on per OTU basis (rows = OTUs, columns = feature counts, first column = OTU names)
#' @param tax_tab (Optional) Data frame with OTU taxonomy (rows = OTUs, columns = taxonomy ranks, first column = OTU names). If provided, taxonomy ranks will be added to the resulting table for each OTU
#' @param features Character vector with function names; if provided, results will be limited to only the specified functions
#' @param NSTI_present Logical; idnicating weather NSTI values are present in the last column of func_tab
#' @param rel_abund Logical; if TRUE, OTU counts will be transformed to relative abundances
#' @param remove_zero_contributions Logical; if TRUE, OTUs with zero contribution will be removed from results
#'
#' @details
#' This function is analogous to the `metagenome_contributions.py` from PICRUSt.
#' Each line in the results relates how much a single OTU (third column) contributes to a single KO (first column) within a single sample (second column).
#' The fifth column contains the actual relative abundance contributed by this OTU, and the other columns contain other information about the abundance of the OTU and the percentage contribution of this OTU.
#' The last columns provide the taxonomy information for the OTU.
#'
#' @return Data frame
#' @export
#' @references https://picrust.github.io/picrust/scripts/metagenome_contributions.html
#'
#' @examples
#' ## Create dummy data
#' set.seed(111)
#' NSAMP=10    # number of samples
#' NSPEC=30    # number of species/OTUs
#' NGENES=20   # number of features
#'
#' dummy_name <- function(len = 5){ paste(sample(letters, size = len, replace = T), collapse = "") }
#'
#' # Table with gene counts per OTU
#' func_tab <- data.frame(
#'               OTU_ID = replicate(n = 100, expr = dummy_name()),
#'               matrix(data = sample(0:100, size = 100*NGENES, replace = T), nrow = 100),
#'               stringsAsFactors = F)
#' colnames(func_tab)[-1] <- paste("F", 1:(ncol(func_tab)-1), sep="")
#'
#' # Table with OTU abundances
#' otu_tab <- data.frame(
#'               OTU = sample(func_tab$OTU_ID, size = NSPEC),
#'               matrix(data = sample(1:10000, size = NSPEC*NSAMP, replace = T), nrow = NSPEC),
#'               stringsAsFactors = F)
#' colnames(otu_tab)[-1] <- paste("Samp", 1:(ncol(otu_tab)-1), sep="")
#'
#' # Table with OTU taxonomy annotation
#' tax_tab <- data.frame(OTU_ID = func_tab$OTU_ID,
#'               Kingdom = sort(sample(LETTERS[1:4], size = nrow(func_tab), replace = T)),
#'               Phylum = sort(sample(LETTERS[5:20], size = nrow(func_tab), replace = T)),
#'               stringsAsFactors = F)
#'
#' metag <- metagenome_contributions(otu_tab, func_tab, tax_tab,
#'                      features = c("F1", "F2", "F3"), NSTI_present = F, rel_abund = FALSE,
#'                      remove_zero_contributions = TRUE)
#' head(metag)
#'
metagenome_contributions <- function(otu_tab, func_tab, tax_tab = NULL, features = NULL, NSTI_present = TRUE, rel_abund = TRUE, remove_zero_contributions = TRUE){

  require(vegan)
  require(reshape2)
  require(plyr)

  # Subset OTUs for the features present in the functional table
  otu_vs_func <- otu_tab[, 1] %in% func_tab[, 1]
  if(any(!otu_vs_func)){
    cat("Warning: OTUs that are not present in the table with functional features will be removed.\n")
    otu_tab <- otu_tab[otu_vs_func, ]
  }
  otus <- as.matrix( otu_tab[, -1] )
  rownames(otus) <- otu_tab[, 1]

  # Extract OTUs from functional table
  match_otus <- match(x = otu_tab[, 1], table = func_tab[,1])
  if(any(is.na(match_otus))){ match_otus <- na.omit(match_otus) }
  funcs <- as.matrix( func_tab[match_otus, -1] )
  rownames(funcs) <- otu_tab[, 1]

  # Remove NSTI scores if present
  if(NSTI_present == TRUE){ funcs <- funcs[, -ncol(funcs)] }

  # Subset functional table
  if(!is.null(features)){
    features_OK <- colnames(funcs) %in% features
    if(any(features_OK)){
        cat("Filtering the genome table to include only user-specified functions.\n")
        if(sum(features_OK) == 1){                    ## only 1 feature - we need to re-create dataframe
            fname <- colnames(funcs)[features_OK]
            funcs <- data.frame(funcs[, which(features_OK)])
            colnames(funcs) <- fname
        }
        if(sum(features_OK) > 1) { funcs <- funcs[, which(features_OK)] }
    } else {
        stop("Error: specified features are absent in the functional table.\n")
    }
  }

  # Standardize OTU counts to relative abundance
  if(rel_abund == TRUE){ otus <- decostand(otus, method = "total", MARGIN = 2) }

  # Reshape data functional table
  if(ncol(funcs) == 1){   # only one feature
    funcs_long <- data.frame(OTU = rownames(funcs), Gene = colnames(funcs), GeneCountPerGenome = funcs[,1], stringsAsFactors = F)
  } else {
    funcs_long <- melt(data = funcs, varnames = c("OTU", "Gene"), value.name = "GeneCountPerGenome")
  }

  # Reshape OTU table
  otus_long <- melt(data = otus, varnames = c("OTU", "Sample"), value.name = "OTUAbundanceInSample")

  # Combine data to a PICRUSt-like table
  res <- plyr::join(x = funcs_long, y = otus_long, by = "OTU", type = "full")

  # Reorder columns as in PICRUSt
  res <- data.frame(
            Gene = res$Gene,
            Sample = res$Sample,
            OTU = res$OTU,
            GeneCountPerGenome = res$GeneCountPerGenome,
            OTUAbundanceInSample = res$OTUAbundanceInSample,
            stringsAsFactors = F)

  # Contribution
  res$CountContributedByOTU <- with(res, OTUAbundanceInSample * GeneCountPerGenome)

  # Remove zero-contributors
  if(remove_zero_contributions == TRUE){
    zeroes <- res$CountContributedByOTU == 0
    if(any(zeroes)){
        res <- res[-which(zeroes), ]
    }
  }

  # Estimate realtive contribution
  rel_contrib <- function(x){ x$CountContributedByOTU / sum(x$CountContributedByOTU) }

  # Percentage of OTU contribution to the gene count within a sample
  res <- ddply(.data=res, .variables=c("Gene", "Sample"), .fun=function(z){
    data.frame(z, ContributionPercentOfSample = rel_contrib(z), stringsAsFactors = F)
    })

  # Percentage of OTU contribution over all samples
  res <- ddply(.data=res, .variables="Gene", .fun=function(z){
    data.frame(z, ContributionPercentOfAllSamples = rel_contrib(z), stringsAsFactors = F)
    })

  # Add taxonomy information for each OTU
  if(!is.null(tax_tab)){
    otu_in_tax <- unique(res$OTU) %in% tax_tab[,1]
    if( !any(otu_in_tax)){
        cat("Warning: not all OTUs from the OTU abundace table are present in the taxonomy table.\n")

        number_of_missing <- sum(!otu_in_tax)
        missing_names <- unique(res$OTU)[ !otu_in_tax ]
        number_of_tax_ranks <- ncol(tax_tab) - 1

        # Add NAs for the missing OTUs
        otu_no_tax <- matrix(rep(NA, number_of_missing*number_of_tax_ranks),
                        nrow = number_of_missing, ncol = number_of_tax_ranks)

        otu_no_tax <- data.frame(OTU = missing_names, otu_no_tax)
        colnames(otu_no_tax) <- colnames(tax_tab)
        tax_tab <- rbind(tax_tab, otu_no_tax)
    }

    res <- cbind(res,
            tax_tab[match(x = res$OTU, table = tax_tab[,1]), -1]
            )
  }  # end of tax_tab

  # Drop rownames
  rownames(res) <- NULL

  return(res)
}
