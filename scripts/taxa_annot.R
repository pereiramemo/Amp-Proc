#!/usr/bin/env Rscript

# This code is a modified version of the DADA2 tutorial: 
# https://benjjneb.github.io/dada2/tutorial.html

###############################################################################
### 1. Def. env
###############################################################################

suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(tidyverse)))

# Get the directory where this script is located
script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  script_dir <- getwd()
}
toolbox <- file.path(script_dir, "toolbox.R")
source(toolbox)

###############################################################################
### 2. Parse command line arguments
###############################################################################

# Function to display help
show_usage <- function() {
  cat("Usage: ./taxa_annot.R <options>\n")
  cat("--help                          print this help\n")
  cat("--input_asv_table CHAR          asv table generated with DADA2 (required)\n")
  cat("--input_fasta CHAR              fasta file with sequences to be annotated (optional, will be generated from ASV table if not provided)\n")
  cat("--output_asv_table CHAR         output asv table with taxonomic annotation (required)\n")
  cat("--method CHAR                   method used to annotate: NBC, NBCandEM, BLAST (default: NBC)\n")
  cat("                                NBC: Naive Bayes Classifier; EM: Exact Matching\n")
  cat("--evalue NUM                    evalue used in BLAST search (default: 1e-10)\n")
  cat("--min_identity NUM              minimum identity used in BLAST search (default: 97)\n")
  cat("--train_db CHAR                 training database to run NBC (default: silva_nr_v138_train_seq.fa.gz)\n")
  cat("--ref_db CHAR                   reference database to run EM (default: silva_species_assignment_v138.fa.gz)\n")
  cat("--blast_db CHAR                 blast formatted database to run BLAST (default: SILVA_138_SSURef_NR99_tax_silva.fasta)\n")
  cat("--taxa_map CHAR                 tsv file mapping silva acc with taxonomy (used when running BLAST)\n")
  cat("--nslots NUM                    number of threads used (default: 12)\n")
  cat("--save_workspace                save R workspace image (default: TRUE)\n")
  cat("--no_save_workspace             disable saving workspace\n")
  cat("--overwrite                     overwrite previous output (default: FALSE)\n")
  quit(status = 0)
}

# Initialize parameters with defaults
INPUT_ASV_TABLE <- NULL
INPUT_FASTA <- NULL
OUTPUT_ASV_TABLE <- NULL
METHOD <- "NBC"
EVALUE <- "1e-10"
MIN_IDENTITY <- "97"
TRAIN_DB <- "silva_nr_v138_train_seq.fa.gz"
REF_DB <- "silva_species_assignment_v138.fa.gz"
BLAST_DB <- "SILVA_138_SSURef_NR99_tax_silva.fasta"
TAXA_MAP <- "taxmap_slv_ssu_ref_nr_138.txt"
NSLOTS <- 12
SAVE_WORKSPACE <- TRUE
OVERWRITE <- FALSE
BLOUT <- NULL
SEQ_MAP <- NULL

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  arg <- args[i]
  
  if (arg == "--help" || arg == "-h") {
    show_usage()
  } else if (arg == "--input_asv_table") {
    INPUT_ASV_TABLE <- args[i + 1]
    i <- i + 1
  } else if (arg == "--input_fasta") {
    INPUT_FASTA <- args[i + 1]
    i <- i + 1
  } else if (arg == "--output_asv_table") {
    OUTPUT_ASV_TABLE <- args[i + 1]
    i <- i + 1
  } else if (arg == "--method") {
    METHOD <- args[i + 1]
    i <- i + 1
  } else if (arg == "--evalue") {
    EVALUE <- args[i + 1]
    i <- i + 1
  } else if (arg == "--min_identity") {
    MIN_IDENTITY <- args[i + 1]
    i <- i + 1
  } else if (arg == "--train_db") {
    TRAIN_DB <- args[i + 1]
    i <- i + 1
  } else if (arg == "--ref_db") {
    REF_DB <- args[i + 1]
    i <- i + 1
  } else if (arg == "--blast_db") {
    BLAST_DB <- args[i + 1]
    i <- i + 1
  } else if (arg == "--taxa_map") {
    TAXA_MAP <- args[i + 1]
    i <- i + 1
  } else if (arg == "--nslots") {
    NSLOTS <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--save_workspace") {
    SAVE_WORKSPACE <- TRUE
  } else if (arg == "--no_save_workspace") {
    SAVE_WORKSPACE <- FALSE
  } else if (arg == "--overwrite") {
    OVERWRITE <- TRUE
  } else {
    cat(sprintf("Warning: Unknown option '%s'\n", arg))
  }
  
  i <- i + 1
}

###############################################################################
### 3. Validate required parameters and inputs
###############################################################################

# Validate required parameters
if (is.null(INPUT_ASV_TABLE)) {
  cat("Error: --input_asv_table is required\n")
  show_usage()
}

if (is.null(OUTPUT_ASV_TABLE)) {
  cat("Error: --output_asv_table is required\n")
  show_usage()
}

# Check input ASV table exists
if (!file.exists(INPUT_ASV_TABLE)) {
  cat(sprintf("Error: Input ASV table '%s' does not exist\n", INPUT_ASV_TABLE))
  quit(status = 1)
}

# Handle output file
if (file.exists(OUTPUT_ASV_TABLE)) {
  if (OVERWRITE) {
    cat(sprintf("Removing existing output file: %s\n", OUTPUT_ASV_TABLE))
    file.remove(OUTPUT_ASV_TABLE)
  } else {
    cat(sprintf("Error: Output file '%s' already exists. Use --overwrite to overwrite\n", OUTPUT_ASV_TABLE))
    quit(status = 1)
  }
}

# Validate method
if (!METHOD %in% c("NBC", "NBCandEM", "BLAST")) {
  cat(sprintf("Error: Invalid method '%s'. Must be one of: NBC, NBCandEM, BLAST\n", METHOD))
  quit(status = 1)
}

###############################################################################
### 4. Create fasta and seq map files if needed
###############################################################################

TMP_FASTA <- FALSE
TMP_BLOUT <- FALSE

if (is.null(INPUT_FASTA) || !file.exists(INPUT_FASTA)) {
  cat("Creating FASTA file from ASV table...\n")
  
  # Read ASV table to extract sequences
  asv_data <- read_csv(file = INPUT_ASV_TABLE, col_names = TRUE, show_col_types = FALSE)
  sequences <- asv_data[[1]]  # First column contains sequences
  
  # Create temporary fasta file
  INPUT_FASTA <- tempfile(pattern = "asvs_", fileext = ".fasta")
  TMP_FASTA <- TRUE
  
  # Write sequences to fasta
  fasta_lines <- character(length(sequences) * 2 - 2)  # -2 to skip header
  for (i in 2:length(sequences)) {  # Start from 2 to skip header
    fasta_lines[(i-1)*2 - 1] <- sprintf(">asv_%d", i-1)
    fasta_lines[(i-1)*2] <- sequences[i]
  }
  writeLines(fasta_lines, INPUT_FASTA)
  
  cat(sprintf("FASTA file created: %s\n", INPUT_FASTA))
}

# Create seq map file (maps sequence IDs to actual sequences)
SEQ_MAP <- tempfile(pattern = "seq_map_", fileext = ".tsv")
cat("Creating sequence map file...\n")

# Read fasta and create mapping
fasta_lines <- readLines(INPUT_FASTA)
headers <- fasta_lines[seq(1, length(fasta_lines), 2)]
sequences <- fasta_lines[seq(2, length(fasta_lines), 2)]
headers <- gsub("^>", "", headers)

seq_map_df <- data.frame(
  qseqid = headers,
  asv = sequences,
  stringsAsFactors = FALSE
)
write_tsv(seq_map_df, SEQ_MAP)

# Create temporary BLOUT file if not specified
if (is.null(BLOUT)) {
  BLOUT <- tempfile(pattern = "blast_", fileext = ".out")
  TMP_BLOUT <- TRUE
}

###############################################################################
### 5. Load and format data
###############################################################################

INPUT_ASV_TABLE <- read_csv(file = INPUT_ASV_TABLE, col_names = T, show_col_types = FALSE) %>%
                   column_to_rownames("X1") %>%
                   as.matrix %>%
                   t

INPUT_ASV_TABLE_tdf <- INPUT_ASV_TABLE %>% 
                       t %>% 
                       as.data.frame %>%
                       rownames_to_column("asv")

###############################################################################
### 6. Run taxa annot: Naive Bayes Classifier (NBC)
###############################################################################

if (METHOD == "NBC") {
  
  print("Running NBC ...")
  TAXA <- assignTaxonomy(seqs = INPUT_ASV_TABLE, 
                         refFasta = TRAIN_DB,
                         outputBootstraps = T,
                         multithread = NSLOTS)  
  
  TAXA <- TAXA %>%
          as.data.frame %>%
          rownames_to_column("asv")
  
  INPUT_ASV_TABLE_ANNOT <- right_join(x = TAXA, y = INPUT_ASV_TABLE_tdf, by = "asv")
  
}

###############################################################################
### 7. Run taxa annot: NBC + Exact Matching (EM)
###############################################################################

if (METHOD == "NBCandEM") {
  
  print("Running NBC ...")
  TAXA <- assignTaxonomy(seqs =  INPUT_ASV_TABLE, 
                         refFasta = TRAIN_DB, 
                         outputBootstraps = T,
                         multithread = NSLOTS)

  print("Running EM ...")
  SPECS <- addSpecies(taxtab = TAXA[[1]], 
                      refFasta = REF_DB)
  
  colnames(SPECS) <- paste("tax", colnames(SPECS), sep = ".")
  
  TAXA <- list(SPECS, TAXA[2]) %>%
          as.data.frame %>%
          rownames_to_column("asv")
  
  INPUT_ASV_TABLE_ANNOT <- right_join(x = TAXA, y = INPUT_ASV_TABLE_tdf, by = "asv")
  
}

###############################################################################
### 8. Run taxa annot: blastn
###############################################################################

if (METHOD == "BLAST") {
 
  print("Running BLAST ...")
  blastn_runner(db = BLAST_DB, 
                input_seqs = INPUT_FASTA,
                blout = BLOUT,
                evalue = EVALUE, 
                min_identity = MIN_IDENTITY, 
                nslots = NSLOTS)
  
  print("Mapping taxonomy to acc ...")
  
  # load blout as df
  BLOUT_DF <- read_tsv(BLOUT, col_names = F, show_col_types = FALSE)
  colnames(BLOUT_DF) <- c("qseqid", "sseqid", "pident", "length", 
                       "mismatch", "gapopen", "qstart", "qend", 
                       "sstart", "send", "evalue", "bitscore")
  
  # load seq map as df
  SEQ_MAP_DF <- read_tsv(SEQ_MAP, col_names = T, show_col_types = FALSE)
  
  # load tax map as df
  TAXA_MAP_DF <- read_tsv(TAXA_MAP, col_names = T, show_col_types = FALSE)
  TAXA_MAP_DF$sseqid <- paste(TAXA_MAP_DF$primaryAccession, TAXA_MAP_DF$start, TAXA_MAP_DF$stop, sep = ".")
  
  # Check if all BLOUT sseqids are in TAXA_MAP
  if (sum(BLOUT_DF$sseqid %in% TAXA_MAP_DF$sseqid) != length(BLOUT_DF$sseqid)) {
    stop("Not all sseqid in TAXA_MAP file")
  }
  
  # Check if all BLOUT qseqids are in SEQ_MAP
  if (sum(BLOUT_DF$qseqid %in% SEQ_MAP_DF$qseqid) != length(BLOUT_DF$sseqid)) {
    stop("Not all sseqid in SEQ_MAP file")
  }
  
  # cross tables
  BLOUT_TAXA_MAPPED <- left_join(x = BLOUT_DF, y = TAXA_MAP_DF, by = "sseqid") 
  SEQ_TAXA_MAPPED <- left_join(x = SEQ_MAP_DF, y = BLOUT_TAXA_MAPPED, by = "qseqid") %>%
                     select(asv, qseqid, sseqid, path, organism_name, pident) 
  
  # add taxonomy to asv table
  INPUT_ASV_TABLE_ANNOT <- right_join(x = SEQ_TAXA_MAPPED, y = INPUT_ASV_TABLE_tdf, by = "asv")
}

###############################################################################
### 9. Save asv annot table
###############################################################################

write.csv(x = INPUT_ASV_TABLE_ANNOT, file = OUTPUT_ASV_TABLE)
print("Output ASV table saved")

###############################################################################
### 10. Save R workspace
###############################################################################

if (SAVE_WORKSPACE) {
  output_dir <- dirname(OUTPUT_ASV_TABLE)
  filename_RData <- file.path(output_dir, ".RData")
  save.image(file = filename_RData)
  cat(sprintf("R workspace saved to: %s\n", filename_RData))
}

###############################################################################
### 11. Cleanup temporary files
###############################################################################

if (TMP_FASTA && file.exists(INPUT_FASTA)) {
  file.remove(INPUT_FASTA)
  cat("Temporary FASTA file removed\n")
}

if (TMP_BLOUT && file.exists(BLOUT)) {
  file.remove(BLOUT)
  cat("Temporary BLAST output file removed\n")
}

if (file.exists(SEQ_MAP)) {
  file.remove(SEQ_MAP)
  cat("Temporary sequence map file removed\n")
}

cat("Taxa annotation completed successfully!\n")
