
#' @title Parse amplicon table
#'
#' @param amp File name with OTU table (generated with VSEARCH or USEARCH)
#' @param tax File name with OTU taxonomic annotations
#' @param tax_type Type of taxonomic annotation ("BLAST", "SINA", "STAMPA")
#' @param add.tax.to.data logical, add taxonomic annotations to the OTU table?
#' @param tax.db SILVA database with taxonomy strings and corresponding tax ranks (e.g. tax_slv_ssu_128.txt)
#'
#' @return list
#' @export
#'
#' @examples
#'
parse_amplicon_table <- function(amp=otu.tab, tax=blast, tax_type="BLAST", add.tax.to.data=TRUE, tax.db){

  # load OTU table
  datt <- read.delim(amp, header = T, stringsAsFactors = F)
  rownames(datt) <- datt$X.OTU.ID
  datt <- datt[, -which(colnames(datt) %in% c("X.OTU.ID"))]       # remove OTU names and total

  # sort by OTU abundance
  datt <- datt[order(rowSums(datt, na.rm = T), decreasing = T), ]

  # load taxonomic annotations (from BLAST)
  if( !is.null(tax) ){

    taxx <- read.delim(tax, header = F, stringsAsFactors = F)

    # extract sequence name from blast.txt table
    tt <- do.call(rbind, strsplit(x = taxx$V1, split = ";"))[,1]
    ttm <- match(rownames(datt), tt)

    ## Parse BLAST results
    if(tax_type == "BLAST"){

      # extract intact taxonomy string (will be splitted to tax ranks)
      tax.string <- taxx$V2[ ttm ]

      # prepare the same looking taxonomy as in VSEARCH
      ttt <- gsub(pattern = " ", replacement = "_", x = taxx$V2)    # replace whitespace
      # ttt <- gsub(pattern = ";", replacement = "|", x = ttt)      # replace semicolon

      # paste seq ID and taxonomy
      ttp <- paste(taxx$V4, ttt, sep=" ")

      # Prepare taxonomic annotation corresponiding to the OTU table
      taxx.OK <- data.frame(Taxon=ttp[ ttm ], E_value = taxx$V3[ ttm ], stringsAsFactors = F)
    }

    ## Parse SINA results
    if(tax_type == "SINA"){

      # extract intact taxonomy string (will be splitted to tax ranks)
      tax.string <- taxx$V3[ ttm ]

      # prepare the same looking taxonomy as in VSEARCH
      ttt <- gsub(pattern = " ", replacement = "_", x = taxx$V3)    # replace whitespace

      # Prepare taxonomic annotation corresponiding to the OTU table
      taxx.OK <- data.frame(Taxon=ttt[ ttm ], Identity = taxx$V2[ ttm ], stringsAsFactors = F)
    }

    ## Parse STAMPA results
    if(tax_type == "STAMPA"){

      # extract taxonomy string (will be splitted to tax ranks) and reorder it
      tax.string <- taxx$V4[ ttm ]

      # replace "N|o|_|h|i|t"
      tax.string[which(tax.string == "N|o|_|h|i|t")] <- "No_hit"

      # prepare the same looking taxonomy as in VSEARCH
      tax.string <- gsub(pattern = "\\+", replacement = " ", x = tax.string)   #
      tax.string <- gsub(pattern = "\\|", replacement = ";", x = tax.string)   # replace field separator | to ;

      # Prepare taxonomic annotation corresponiding to the OTU table
      taxx.OK <- data.frame(Taxon=tax.string, Identity = taxx$V3[ ttm ], stringsAsFactors = F)

      # remove ;* (stampa ambiguous assignment) prior to parsing
      tax.string <- gsub(pattern = ";\\*", replacement = "", x = tax.string)
    }

    if(add.tax.to.data == TRUE){
      datt <- data.frame(taxx.OK, datt, stringsAsFactors = F)
    }

    # Parse taxonomy
    cat("..Parsing taxonomy\n")
    tax.table <- silva_tax_parse_batch(tax.string, tax.db)
    rownames(tax.table) <- rownames(datt)
    cat("...Done\n")
  }

  res <- list()
  res$datt <- datt
  res$taxonomy <- tax.table
  res$tax.string <- tax.string
  return(res)
}
### Example:
## BLAST
# zz <- parse_amplicon_table(
#             amp="~/04.VSEARCH/_OTU_table.txt",
#             tax="~/05.TaxAnnot/BLAST_SILVA128/_OTU_centroids_noChimera_tax_assignments.txt")
# str(zz)

## SINA
# zz <- parse_amplicon_table(
#             amp="~/04.VSEARCH/_OTU_table.txt",
#             tax="~/05.TaxAnnot/SINA_SILVA/VSEARCH_SINA_LCA_res.txt",
#             tax_type="SINA")
# str(zz)

## STAMPA
# zz <- parse_amplicon_table(
#             amp="~/04.VSEARCH/_OTU_table.txt",
#             tax="~/05.TaxAnnot/STAMPA_SILVA/representatives.results",
#             tax_type="STAMPA")
# str(zz)
