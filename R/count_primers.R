
## Count primers (counts number of reads in which the primer is found)
# https://benjjneb.github.io/dada2/ITS_workflow.html
count_primers <- function(fq,
  FWD = "GTGARTCATCGAATCTTTG", REV = "TCCTCCGCTTATTGATATGC",
  mismatch = 0, nreads = NULL){

  ## Function to load data
  read_fq <- function(fn, subs = NULL){
    # e.g., fn = "R1.fastq.gz"

    ## Load data
    if(grepl(pattern = "fastq$|fastq.gz$|fq$|fq.gz$", x = fn)){
      ## Read FASTQ file (-> ShortReadQ)
      ff <- ShortRead::readFastq(fn)
    }
    if(grepl(pattern = "fasta$|fasta.gz$|fa$|fa.gz$", x = fn)){
      ## Read FASTA file (-> ShortRead)
      ff <- ShortRead::readFasta(fn)
    }

    ## Subset reads
    if(!is.null(subs)){
      if(subs > length(ff)){ subs <- length(ff) } # limit number of reads to all reads
      ff <- ff[1:subs]
    }

    ## Convert to DNAStringSet
    ff <- ShortRead::sread(ff)

    return(ff)
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
  fq1 <- read_fq(fq[1], subs = nreads)
  if( length(fq) == 2 ){
    fq2 <- read_fq(fq[2], subs = nreads)
  }

  ## Count number of primer occurrences
  if( length(fq) == 2 ){
    rez <- rbind(
      FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fq1), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fq2), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fq1), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fq2)
      )
  }
  if( length(fq) == 1 ){
    rez <- rbind(
      FWD = sapply(FWD.orients, primerHits, fn = fq1), 
      REV = sapply(REV.orients, primerHits, fn = fq1)
      )
  }

  return(rez)
}
# e.g., count_primers(fq = c("s1_R1.fastq.gz", "s1_R2.fastq.gz"))
# e.g., count_primers(fq = c("s1_R1.fasta", "s1_R2.fasta"))
