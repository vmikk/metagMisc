## SILVA parsing functions


#' @title Split single SILVA taxonomic annotation to taxonomic ranks.
#'
#' @param x Character string with taxonomic annotation from BLAST.
#' @param db Data frame with rank designations for all taxonomic paths in the SILVA taxonomy.
#'
#' @details This function will split the taxonomic annotation based on the SILVA database and assign a corresponding taxonomic rank to the obtained parts.
#' @return Named charcter vector with taxonomic annotation (names = rank designations).
#' @export
#' @seealso \code{\link{parse_silva_tax_batch}}
#' @examples
#' # Download taxonomic rank designations for all taxonomic paths used in the SILVA taxonomies
#' tax.db <- read.delim("https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/taxonomy/tax_slv_ssu_128.txt", header = F, stringsAsFactors = F)
#' colnames(tax.db) <- c("path", "taxid", "rank", "remark", "release")
#'
#' x <- "Eukaryota;SAR;Alveolata;Ciliophora;Intramacronucleata;Conthreep;Oligohymenophorea;Peritrichia;Telotrochidium;uncultured eukaryote"
#' parse_silva_tax(x, tax.db)
#'
# Function to parse taxonomy from SILVA db
parse_silva_tax <- function(x, db){

  # Split the search string
  xs <- strsplit(x, split = ";")[[1]]

  # Generate path and add semicolon to the end of string
  pp <- vector()
  for(i in 1:length(xs)){
    pp[i] <- paste( paste(xs[1:i], collapse=";"), ";", sep="")
  }

  # Find path in db and assign taxonomic level
  tl <- aaply(.data = pp, .margins = 1, .fun = function(z){ db$rank[ match(table = db$path, x = z) ] })
  names(xs) <- tl

  # Assign the last column with 'species' if it is NA
  if( is.na(names(xs)[length(xs)]) ){ names(xs)[length(xs)] <- "species" }

  # Do we need to replace other NAs with something ??

  return(xs)
}


#' @title Split vector of SILVA taxonomic annotations and prepare taxonomy table.
#' @param x Vector of character strings with taxonomic annotation (assigned by BLAST)
#' @param db Data frame with rank designations for all taxonomic paths in the SILVA taxonomy
#'
#' @return Data frame with OTUs or species as rows and their taxonomic ranks as columns.
#' @export
#' @seealso \code{\link{parse_silva_tax}}
#' @examples
#' # Download taxonomic rank designations for all taxonomic paths used in the SILVA taxonomies
#' tax.db <- read.delim("https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/taxonomy/tax_slv_ssu_128.txt", header = F, stringsAsFactors = F)
#' colnames(tax.db) <- c("path", "taxid", "rank", "remark", "release")
#'
#' x <- c("Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);Eumetazoa;Bilateria;Platyhelminthes;Catenulida;Catenulidae;Paracatenula galateia",
#' "Eukaryota;SAR;Rhizaria;Cercozoa;Imbricatea;Silicofilosea;Euglyphida;Assulinidae;Assulina;Assulina muscorum",
#' "Eukaryota;Opisthokonta;Nucletmycea;Fungi;Cryptomycota;LKM11;uncultured eukaryote",
#' "Eukaryota;Opisthokonta;Nucletmycea;Fungi;Dikarya;Basidiomycota;Ustilaginomycotina;Exobasidiomycetes;Malasseziales;Incertae Sedis;Malassezia;Malassezia restricta CBS 7877",
#' "Eukaryota;Opisthokonta;Nucletmycea;Fungi;Cryptomycota;LKM11;uncultured fungus",
#' "Eukaryota;SAR;Alveolata;Ciliophora;Intramacronucleata;Conthreep;Oligohymenophorea;Peritrichia;Telotrochidium;uncultured eukaryote",
#' "No blast hit", "Eukaryota;SAR;Alveolata;Ciliophora;Intramacronucleata;Armophorea;Armophorida;Metopus;Metopus striatus",
#' "Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);Eumetazoa;Bilateria;Rotifera;Monogononta;Ploimida;Lepadella patella",
#' "Eukaryota;SAR;Rhizaria;Cercozoa;Novel Clade 10;uncultured eukaryote"
#' )
#' parse_silva_tax_batch(x, tax.db)
#'
parse_silva_tax_batch <- function(x, db){

  require(plyr)

  # Prepare list of taxonomic assignments
  res <- alply(.data = x, .margins = 1, .fun = parse_silva_tax, db = db)

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
