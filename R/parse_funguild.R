
#' @title Parse FUNGuild database
#' @description FUNGuild = Fungi + fUNctional + Guild
#' @param url Funguild databse link
#' @param tax_name Logical, change numeral taxonomic coding with the corresponding taxonomic levels
#'
#' @return Data frame containing the entire FUNGuild database.
#' @export
#' @references
#' Nguyen NH, Song Z, Bates ST, Branco S, Tedersoo L, Menke J, Schilling JS, Kennedy PG. 2016. FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild. Fungal Ecology 20:241–248.
#' https://github.com/UMNFuN/FUNGuild
#' @examples
#' fg <- parse_funguild()
#'
parse_funguild <- function(url = 'http://www.stbates.org/funguild_db.php', tax_name = TRUE){

  require(XML)
  require(jsonlite)
  # require(RCurl)

  ## Parse data
  tmp <- htmlParse(url)                           # [XML]
  tmp <- xpathSApply(tmp, "//body", xmlValue)

  ## Read url and convert to data.frame
  db <- fromJSON(txt=tmp)

  ## Remove IDs
  db$`_id` <- NULL

  if(tax_name == TRUE){

    ## Code legend
    ## Taxon Level: A numeral corresponding the correct taxonomic level for the taxon
    taxons <- c(
      "keyword",                                                       # 0
      "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", # 3:8
      "Family", "Subfamily", "Tribe", "Subtribe", "Genus",             # 9:13
      "Subgenus", "Section", "Subsection", "Series", "Subseries",      # 15:19
      "Species", "Subspecies", "Variety", "Subvariety", "Form",        # 20:24
      "Subform", "Form Species")

    ## Table with coding
    taxmatch <- data.frame(
      TaxID = c(0, 3:13, 15:26),
      Taxon = factor(taxons, levels = taxons))

    ## Match taxon codes
    db$taxonomicLevel <- taxmatch[match(x = db$taxonomicLevel, table = taxmatch$TaxID), "Taxon"]
  }

  # remove rows with missing data
  # which(
  # 	with(db, trophicMode == "NULL" & guild == "NULL" & growthForm == "NULL" & trait == "NULL" & notes == "NULL")
  # 	)

  return(db)
}
