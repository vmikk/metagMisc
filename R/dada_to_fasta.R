
#' @title Write DADA sequences to fasta.
#' @description Rename sequences to their hash values and write them to fasta file.
#' @param seqtab dada-class or derep-class objects
#' @param out Output file name (fasta)
#' @param hash Hash function to use: "sha1" (default), "sha256", "md5"
#' @param ... Additional parameters passed on to \code{\link{uniquesToFasta}}
#'
#' @details
#' This function relabels sequences using diffetent message digest algorithms applied to each sequence. This approach guarantees (with a very high probability) that FASTA entries from different projects with identical names will also have identical sequences.
#' MD5 algorithm generates a 128-bit digest that is represented by 32 hexadecimal characters.
#' SHA1 generates a 160-bit digest that is represented by 40 hexadecimal characters.
#' SHA256 generates a 256-bit digest that is represented by 64 hexadecimal characters.
#' The probability of a collision (two non-identical sequences resulting in the same digest) is smaller for the SHA-algorithms than it is for the MD5 algorithm.
#' Default hash function is SHA1 which should produce identical results with "--relabel_sha1" function of VSEARCH.
#'
#' @return Invisible returns sequence names in the VSEARCH/USEARCH style.
#' @export
#' @seealso \code{\link{uniquesToFasta}}, \code{\link{getUniques}}
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, err=tperr1)
#' seqtab <- getUniques(dada1)
#' dada_to_fasta(seqtab, out = "DADA2.fasta", hash = "sha1")
#'
dada_to_fasta <- function(seqtab, out = "DADA2.fasta", hash = "sha1", ...){

  require(dada2)
  require(openssl)

  # prepare sequence names in USEARCH and VSEARCH-style
  seq_uniq <- getUniques(seqtab)   # integer vector named by unique sequence and valued by abundance.

  if(hash == "sha1"){ hh <- sha1(names(seq_uniq)) }
  if(hash == "sha256"){ hh <- sha256(names(seq_uniq)) }
  if(hash == "md5"){ hh <- md5(names(seq_uniq)) }

  seq_names <- paste(as.character(hh),
                  ";size=",
                  seq_uniq,
                  ";",
                  sep="")

  # Export sequence as fasta
  uniquesToFasta(seq_uniq, fout = out, ids = seq_names, ...)

  invisible(seq_names)
}
