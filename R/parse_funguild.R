
#' @title Parse FUNGuild database
#' @description FUNGuild = Fungi + fUNctional + Guild
#' @param url Funguild databse link
#'
#' @return Data frame containing the entire FUNGuild database.
#' @export
#' @references
#' Nguyen NH, Song Z, Bates ST, Branco S, Tedersoo L, Menke J, Schilling JS, Kennedy PG. 2016. FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild. Fungal Ecology 20:241–248.
#' https://github.com/UMNFuN/FUNGuild
#' @examples
#' fg <- parse_funguild()
#'
parse_funguild <- function(url = 'http://www.stbates.org/funguild_db.php'){

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

  # remove rows with missing data
  # which(
  # 	with(db, trophicMode == "NULL" & guild == "NULL" & growthForm == "NULL" & trait == "NULL" & notes == "NULL")
  # 	)

  return(db)
}
