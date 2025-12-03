###############################################################################
### 1. Fun to create all orientations of the input sequence
###############################################################################

all_orients <- function(PRIMER) {
  require(Biostrings)
  
  DNA <- DNAString(PRIMER)  # The Biostrings works w/ DNAString objects rather than character vectors
  
  orients <- c(Forward = DNA, 
               Complement = complement(DNA), 
               Reverse = reverse(DNA), 
               RevComp = reverseComplement(DNA))
  
  return(sapply(orients, toString))  # Convert back to character vector
}

###############################################################################
### 2. Fun to count the number of reads in which the primer is found
###############################################################################

primer_hits <- function(PRIMER, INPUT) {
  
  nhits <- vcountPattern(pattern = PRIMER, 
                         subject = sread(readFastq(INPUT)), 
                         fixed = FALSE)
  return(sum(nhits > 0))
}


primer_hits_fasta <- function(PRIMER, INPUT) {
  
  nhits <- vcountPattern(pattern = PRIMER, 
                         subject = sread(readFasta(INPUT)), 
                         fixed = FALSE)
  return(sum(nhits > 0))
}


###############################################################################
### 3. Fun to run cutadapt from R
###############################################################################

cutadapt <- "/home/epereira/.local/bin/cutadapt"

cutadapt_runner <- function(primer_fwd = PRIMER_FWD, 
                            primer_rev = PRIMER_REV, 
                            input_r1 = INPUT_R1, 
                            input_r2 = INPUT_R2, 
                            output_r1 = OUTPUT_R1, 
                            output_r2 = OUTPUT_R2,
                            nslots = NSLOTS) {

  primer_fwd_rc <- dada2::rc(primer_fwd) # rc is used instead of reverseComplement
  primer_rev_rc <- dada2::rc(primer_rev) # to avoid converting data to DNAstring
  
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1_flags <- paste("-g", primer_fwd, "-a", primer_rev_rc) 
  
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  R2_flags <- paste("-G", primer_rev, "-A", primer_fwd_rc)
  
  # Run Cutadapt
  for(i in seq_along(input_r1)) {
    
    system2(cutadapt, 
            args = c(R1_flags, R2_flags, 
                     "-n", 2,
                     "--minimum-length 1",
                     "--cores", nslots,
                     "-o", output_r1[i], 
                     "-p", output_r2[i], 
                     input_r1[i], 
                     input_r2[i]))
  }
}

###############################################################################
### 4. Fun to run blastn from R
###############################################################################

# "2>/dev/null"))

blastn <- "/home/bioinf/bin/blast/blast_v2.10.0+/ncbi-blast-2.10.0+/bin/blastn"

blastn_runner <- function(db = BLAST_DB, input_seqs = INPUT_SEQS, evalue = EVALUE, 
                          min_identity = MIN_IDENTITY, blout = BLOUT, nslots = NSLOTS) {
  
  system2(blastn,
          args = c("-db", db,
                   "-query", input_seqs,
                   "-evalue", evalue,
                   "-perc_identity", min_identity,
                   "-outfmt", 6,
                   "-out", blout,
                   "-num_alignments", 1,
                   "-num_threads", nslots))
                  
}

