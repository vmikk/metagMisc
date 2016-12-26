## SILVA parsing functions


#' @title Split single SILVA taxonomic annotation to taxonomic ranks.
#'
#' @param x Character string with taxonomic annotation.
#' @param db Data frame with SILVA taxonomy.
#'
#' @return
#' @export
#'
#' @examples
#'
# Function to parse taxonomy from SILVA db
silva_tax_parse <- function(x, db){
  # x = character string from BLAST
  # db = data frame with SILVA taxonomy

  # Split the search string
  xs <- strsplit(x, split = ";")[[1]]

  # Generate path and add semicolon to the end of string
  pp <- vector()
  for(i in 1:length(xs)){
    pp[i] <- paste( paste(xs[1:i], collapse=";"), ";", sep="")
  }

  # Find path in db and assign taxonomic level
  tl <- aaply(.data = pp, .margins = 1, .fun = function(z){ tax.db$rank[ match(table = tax.db$path, x = z) ] })
  names(xs) <- tl

  # Assign the last column with 'species' if it is NA
  if( is.na(names(xs)[length(xs)]) ){ names(xs)[length(xs)] <- "species" }

  # Do we need to replace other NAs with something ??

  return(xs)
}
## Example
# x <- "Eukaryota;SAR;Alveolata;Ciliophora;Intramacronucleata;Conthreep;Oligohymenophorea;Peritrichia;Telotrochidium;uncultured eukaryote"
# silva.tax.parse(x, tax.db)



#' Title
#' @title Split vector of SILVA taxonomic annotations and prepare taxonomy table.
#' @param x Vector of character strings with taxonomic annotation.
#' @param db Data frame with SILVA taxonomy
#'
#' @return
#' @export
#'
#' @examples
#'
# Prepare taxonomy table
silva_tax_parse_batch <- function(x, db){
  # x = vector of character strings (taxonomy assigned by BLAST)
  # db = data frame with SILVA taxonomy

  # Prepare list of taxonomic assignments
  res <- alply(.data = x, .margins = 1, .fun = silva_tax_parse)

  # Convert each vector to matrix
  res <- llply(.data = res, .fun = function(x){ t(as.matrix(x)) })

  # Prepare table
  res <- data.frame( do.call(rbind.fill.matrix, res) )

  # Sort columns
  tax.levels <- c(
    "domain", "major_clade",
    "superkingdom", "kingdom", "subkingdom", "infrakingdom",
    "superphylum", "phylum", "subphylum", "infraphylum",
    "superclass", "class", "subclass", "infraclass",
    "superorder", "order", "suborder",
    "superfamily", "family", "subfamily",
    "genus", "species")

  tax.levels <- data.frame(Level = tax.levels, Rank = 1:length(tax.levels))

  mm <- match(x = colnames(res), table = tax.levels$Level)
  res <- res[, as.character(tax.levels$Level[sort(mm)]) ]

  return(res)
}
## Example
# x <- c("Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);Eumetazoa;Bilateria;Platyhelminthes;Catenulida;Catenulidae;Paracatenula galateia",
# "Eukaryota;SAR;Rhizaria;Cercozoa;Imbricatea;Silicofilosea;Euglyphida;Assulinidae;Assulina;Assulina muscorum",
# "Eukaryota;Opisthokonta;Nucletmycea;Fungi;Cryptomycota;LKM11;uncultured eukaryote",
# "Eukaryota;Opisthokonta;Nucletmycea;Fungi;Dikarya;Basidiomycota;Ustilaginomycotina;Exobasidiomycetes;Malasseziales;Incertae Sedis;Malassezia;Malassezia restricta CBS 7877",
# "Eukaryota;Opisthokonta;Nucletmycea;Fungi;Cryptomycota;LKM11;uncultured fungus",
# "Eukaryota;SAR;Alveolata;Ciliophora;Intramacronucleata;Conthreep;Oligohymenophorea;Peritrichia;Telotrochidium;uncultured eukaryote",
# "No blast hit", "Eukaryota;SAR;Alveolata;Ciliophora;Intramacronucleata;Armophorea;Armophorida;Metopus;Metopus striatus",
# "Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);Eumetazoa;Bilateria;Rotifera;Monogononta;Ploimida;Lepadella patella",
# "Eukaryota;SAR;Rhizaria;Cercozoa;Novel Clade 10;uncultured eukaryote"
# )
# silva_tax_parse_batch(x, tax.db)


