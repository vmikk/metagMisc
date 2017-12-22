
#' @title Parse taxonomy string from AMPtk (for single OTU)
#'
#' @param x Character vector of length 1 with taxonomic annotation of single OTU/species.
#' @details This function will split
#' @return Named character vector with taxonomic ranks of OTU or species.
#' @export
#' @seealso \code{\link{parse_taxonomy_qiime}}
#' @references AMPtk: Amplicon ToolKit for NGS data (formally UFITS): https://github.com/nextgenusfs/amptk
#' @examples
#' parse_taxonomy_amptk("k:Fungi,p:Zygomycota,o:Mortierellales,f:Mortierellaceae,g:Mortierella,s:Mortierella parvispora")
#' parse_taxonomy_amptk("k:Fungi,p:Ascomycota,g:Chalara")
#'
parse_taxonomy_amptk <- function(x){
    
    # require(phyloseq)

    ## Convert ufits taxonomic ranks to QIIME-like style
    x <- gsub(pattern = ":", replacement = "__", x = x)
    x <- gsub(pattern = ",", replacement = ";", x = x)
    res <- phyloseq::parse_taxonomy_qiime(x)
    return(res)
}



#' @title Parse multiple taxonomy strings from AMPtk.
#'
#' @param x Vector of character strings from AMPtk taxonomy
#' @details
#' This function splits AMPtk (ex-UFITS) result of taxonomic annotation by taxonomic rank and combines it into one table.
#' The first part of taxonomy string is assumed to be an ID of the best hit (e.g., EF040844) or a taxonomic classificator (UTAX or SINTAX) followed by a semicolon.
#'
#' @return Data frame with OTUs or species as rows and their taxonomic ranks as columns.
#' @references AMPtk: Amplicon ToolKit for NGS data (formally UFITS): https://github.com/nextgenusfs/amptk
#' @seealso \code{\link{parse_taxonomy_amptk}} \code{\link{parse_taxonomy_qiime}}
#' @export
#' @examples
#' ## Taxonomy with sequence ID
#' tax <- c(
#' "UDB011192;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Russulaceae,g:Russula,s:Russula chloroides",
#' "SINTAX;k:Fungi,p:Chytridiomycota,c:Chytridiomycetes,o:Rhizophydiales",
#' "EF040844;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Atheliales,f:Atheliaceae,g:Piloderma",
#' "EU240039;k:Fungi,p:Zygomycota,o:Mortierellales,f:Mortierellaceae,g:Mortierella",
#' "UTAX;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Atheliales,f:Atheliaceae,g:Amphinema,s:Amphinema byssoides",
#' "SINTAX;k:Fungi,p:Chytridiomycota,c:Chytridiomycetes,o:Spizellomycetales"
#' )
#'
#' parse_taxonomy_amptk_batch(tax)
#'
#' ## Taxonomy without sequence ID
#' tax2 <- c(
#' "k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Russulaceae,g:Russula,s:Russula chloroides",
#' "k:Fungi,p:Chytridiomycota,c:Chytridiomycetes,o:Rhizophydiales",
#' "k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Atheliales,f:Atheliaceae,g:Piloderma",
#' "k:Fungi,p:Zygomycota,o:Mortierellales,f:Mortierellaceae,g:Mortierella",
#' "k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Atheliales,f:Atheliaceae,g:Amphinema,s:Amphinema byssoides",
#' "k:Fungi,p:Chytridiomycota,c:Chytridiomycetes,o:Spizellomycetales"
#' )
#'
#' parse_taxonomy_amptk_batch(tax2, withID = FALSE)
#'
parse_taxonomy_amptk_batch <- function(x, withID = TRUE){
    
    # require(plyr)

    if(withID == TRUE){
      ## Split OTUId to SequenceId (e.g., JQ976006) or MethodId (UTAX,SINTAX) + Taxonomy
      res <- strsplit(x, split = ";")
      res <- do.call(rbind, res)

      # Prepare list of taxonomic assignments
      res <- plyr::alply(.data = res[, 2], .margins = 1, .fun = parse_taxonomy_amptk)
    } else {
      # Prepare list of taxonomic assignments
      res <- plyr::alply(.data = x, .margins = 1, .fun = parse_taxonomy_amptk)
    }

    # Convert each vector to matrix
    res <- plyr::llply(.data = res, .fun = function(x){ t(as.matrix(x)) })

    # Prepare table
    res <- data.frame( do.call(rbind.fill.matrix, res), stringsAsFactors = F)

    # Sort columns
    tax.levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    tax.levels <- data.frame(Level = tax.levels, Rank = 1:length(tax.levels))

    mm <- match(x = colnames(res), table = tax.levels$Level)
    res <- res[, as.character(tax.levels$Level[sort(mm)]) ]

    return(res)
}
