
## Count primers (counts number of reads in which the primer is found)
# https://benjjneb.github.io/dada2/ITS_workflow.html
count_primers <- function(fq,
  FWD = "GTGARTCATCGAATCTTTG", REV = "TCCTCCGCTTATTGATATGC",
  mismatch = 0){

  ## Function to load data
  read_fq <- function(fn){
    # e.g., fn = "R1.fastq.gz"

    ## Read FASTQ file (-> ShortReadQ) and convert to DNAStringSet
    fn <- ShortRead::sread( ShortRead::readFastq(fn) )
    return(fn)
  }

  ## Function to count number of matches
  primerHits <- function(primer, fn, mism = 0) {
    # e.g., primer = "GTGARTCATCGAATCTTTG", fn = read_fq("R1.fastq.gz")

    ## Search primer
    nhits <- Biostrings::vcountPattern(primer, fn, max.mismatch = mism, fixed = FALSE)
    
    return(sum(nhits > 0))
  }

  ## Create all orientations of the input sequence
  allOrients <- function(primer) {
    dna <- Biostrings::DNAString(primer)
    orients <- c(
      Forward = dna,
      Complement = complement(dna),
      Reverse = reverse(dna),
      RevComp = reverseComplement(dna)
      )

    ## Convert back to character vector
    return(sapply(orients, toString))  
  }

  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)

  ## Load data
  fq1 <- read_fq(fq[1])
  fq2 <- read_fq(fq[2])

  ## Count number of primer occurrences
  rez <- rbind(
    FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fq1), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fq2), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fq1), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fq2)
    )

  return(rez)
}
# e.g., count_primers(fq = c("s1_R1.fastq.gz", "s1_R2.fastq.gz"))