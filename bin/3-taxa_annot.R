#!/usr/bin/env Rscript

# This code is a modified version of the DADA2 tutorial:
# https://benjjneb.github.io/dada2/tutorial.html

###############################################################################
### 1. Def. env
###############################################################################

suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(tidyverse)))

script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  script_dir <- getwd()
}
source(file.path(script_dir, "toolbox.R"))

###############################################################################
### 2. Parse command line arguments
###############################################################################

show_usage <- function() {
  cat("Usage: ./taxa_annot.R <options>\n")
  cat("--help                          print this help\n")
  cat("--input_asv_table CHAR          asv table generated with DADA2 (required)\n") # nolint
  cat("--input_fasta CHAR              fasta file with sequences to annotate (optional)\n") # nolint
  cat("--output_asv_table CHAR         output asv table with taxonomic annotation (required)\n") # nolint
  cat("--method CHAR                   annotation method: NBC, NBCandEM, BLAST (default: NBC)\n") # nolint
  cat("                                NBC: Naive Bayes Classifier; EM: Exact Matching\n") # nolint
  cat("--evalue NUM                    evalue used in BLAST search (default: 1e-10)\n") # nolint
  cat("--min_identity NUM              minimum identity used in BLAST search (default: 97)\n") # nolint
  cat("--train_db CHAR                 training database to run NBC (default: silva_nr_v138_train_seq.fa.gz)\n") # nolint
  cat("--ref_db CHAR                   reference database to run EM (default: silva_species_assignment_v138.fa.gz)\n") # nolint
  cat("--blast_db CHAR                 blast formatted database to run BLAST (default: SILVA_138_SSURef_NR99_tax_silva.fasta)\n") # nolint
  cat("--taxa_map CHAR                 tsv file mapping silva acc with taxonomy (used with BLAST)\n") # nolint
  cat("--nslots NUM                    number of threads used (default: 12)\n")
  cat("--save_workspace                save R workspace image (default: TRUE)\n") # nolint
  cat("--no_save_workspace             disable saving workspace\n") # nolint
  cat("--overwrite                     overwrite previous output (default: FALSE)\n") # nolint
  quit(status = 0)
}

input_asv_table <- NULL
input_fasta <- NULL
output_asv_table <- NULL
method <- "NBC"
evalue <- "1e-10"
min_identity <- "97"
train_db <- "silva_nr_v138_train_seq.fa.gz"
ref_db <- "silva_species_assignment_v138.fa.gz"
blast_db <- "SILVA_138_SSURef_NR99_tax_silva.fasta"
taxa_map <- "taxmap_slv_ssu_ref_nr_138.txt"
nslots <- 12
save_workspace <- TRUE
overwrite <- FALSE
blout <- NULL
seq_map <- NULL

args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  arg <- args[i]

  if (arg == "--help" || arg == "-h") {
    show_usage()
  } else if (arg == "--input_asv_table") {
    input_asv_table <- args[i + 1]
    i <- i + 1
  } else if (arg == "--input_fasta") {
    input_fasta <- args[i + 1]
    i <- i + 1
  } else if (arg == "--output_asv_table") {
    output_asv_table <- args[i + 1]
    i <- i + 1
  } else if (arg == "--method") {
    method <- args[i + 1]
    i <- i + 1
  } else if (arg == "--evalue") {
    evalue <- args[i + 1]
    i <- i + 1
  } else if (arg == "--min_identity") {
    min_identity <- args[i + 1]
    i <- i + 1
  } else if (arg == "--train_db") {
    train_db <- args[i + 1]
    i <- i + 1
  } else if (arg == "--ref_db") {
    ref_db <- args[i + 1]
    i <- i + 1
  } else if (arg == "--blast_db") {
    blast_db <- args[i + 1]
    i <- i + 1
  } else if (arg == "--taxa_map") {
    taxa_map <- args[i + 1]
    i <- i + 1
  } else if (arg == "--nslots") {
    nslots <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--save_workspace") {
    save_workspace <- TRUE
  } else if (arg == "--no_save_workspace") {
    save_workspace <- FALSE
  } else if (arg == "--overwrite") {
    overwrite <- TRUE
  } else {
    cat(sprintf("Warning: Unknown option '%s'\n", arg))
  }

  i <- i + 1
}

###############################################################################
### 3. Validate required parameters and inputs
###############################################################################

if (is.null(input_asv_table)) {
  cat("Error: --input_asv_table is required\n")
  show_usage()
}

if (is.null(output_asv_table)) {
  cat("Error: --output_asv_table is required\n")
  show_usage()
}

if (!file.exists(input_asv_table)) {
  cat(sprintf("Error: Input ASV table '%s' does not exist\n", input_asv_table))
  quit(status = 1)
}

if (file.exists(output_asv_table)) {
  if (overwrite) {
    cat(sprintf("Removing existing output file: %s\n", output_asv_table))
    file.remove(output_asv_table)
  } else {
    cat(sprintf(
      "Error: Output file '%s' already exists. Use --overwrite to overwrite\n",
      output_asv_table
    ))
    quit(status = 1)
  }
}

if (!method %in% c("NBC", "NBCandEM", "BLAST")) {
  cat(sprintf(
    "Error: Invalid method '%s'. Must be one of: NBC, NBCandEM, BLAST\n",
    method
  ))
  quit(status = 1)
}

###############################################################################
### 4. Create fasta and seq map files if needed
###############################################################################

tmp_fasta <- FALSE
tmp_blout <- FALSE

if (is.null(input_fasta) || !file.exists(input_fasta)) {
  cat("Creating FASTA file from ASV table...\n")

  asv_data <- read_csv(
    file = input_asv_table, col_names = TRUE, show_col_types = FALSE
  )
  sequences <- asv_data[[1]]

  input_fasta <- tempfile(pattern = "asvs_", fileext = ".fasta")
  tmp_fasta <- TRUE

  fasta_content <- character((length(sequences) - 1) * 2)
  for (idx in 2:length(sequences)) {
    fasta_content[(idx - 1) * 2 - 1] <- sprintf(">asv_%d", idx - 1)
    fasta_content[(idx - 1) * 2] <- sequences[idx]
  }
  writeLines(fasta_content, input_fasta)

  cat(sprintf("FASTA file created: %s\n", input_fasta))
}

seq_map <- tempfile(pattern = "seq_map_", fileext = ".tsv")
cat("Creating sequence map file...\n")

fasta_lines <- readLines(input_fasta)
headers <- fasta_lines[seq(1, length(fasta_lines), 2)]
sequences <- fasta_lines[seq(2, length(fasta_lines), 2)]
headers <- gsub("^>", "", headers)

seq_map_df <- data.frame(qseqid = headers, asv = sequences)
write_tsv(seq_map_df, seq_map)

if (is.null(blout)) {
  blout <- tempfile(pattern = "blast_", fileext = ".out")
  tmp_blout <- TRUE
}

###############################################################################
### 5. Load and format data
###############################################################################

asv_matrix <- read_csv(
  file = input_asv_table, col_names = TRUE, show_col_types = FALSE
) |>
  column_to_rownames("X1") |>
  as.matrix() |>
  t()

asv_tdf <- asv_matrix |>
  t() |>
  as.data.frame() |>
  rownames_to_column("asv")

###############################################################################
### 6. Run taxa annot: Naive Bayes Classifier (NBC)
###############################################################################

if (method == "NBC") {

  print("Running NBC ...")
  taxa <- assignTaxonomy(
    seqs = asv_matrix,
    refFasta = train_db,
    outputBootstraps = TRUE,
    multithread = nslots
  )

  taxa <- taxa |>
    as.data.frame() |>
    rownames_to_column("asv")

  asv_table_annot <- right_join(x = taxa, y = asv_tdf, by = "asv")

}

###############################################################################
### 7. Run taxa annot: NBC + Exact Matching (EM)
###############################################################################

if (method == "NBCandEM") {

  print("Running NBC ...")
  taxa <- assignTaxonomy(
    seqs = asv_matrix,
    refFasta = train_db,
    outputBootstraps = TRUE,
    multithread = nslots
  )

  print("Running EM ...")
  specs <- addSpecies(taxtab = taxa[[1]], refFasta = ref_db)

  colnames(specs) <- paste("tax", colnames(specs), sep = ".")

  taxa <- list(specs, taxa[2]) |>
    as.data.frame() |>
    rownames_to_column("asv")

  asv_table_annot <- right_join(x = taxa, y = asv_tdf, by = "asv")

}

###############################################################################
### 8. Run taxa annot: blastn
###############################################################################

if (method == "BLAST") {

  print("Running BLAST ...")
  blastn_runner(
    db = blast_db,
    input_seqs = input_fasta,
    blout = blout,
    evalue = evalue,
    min_identity = min_identity,
    nslots = nslots
  )

  print("Mapping taxonomy to acc ...")

  blout_df <- read_tsv(blout, col_names = FALSE, show_col_types = FALSE)
  colnames(blout_df) <- c(
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
  )

  seq_map_df <- read_tsv(seq_map, col_names = TRUE, show_col_types = FALSE)

  taxa_map_df <- read_tsv(taxa_map, col_names = TRUE, show_col_types = FALSE)
  taxa_map_df$sseqid <- paste(
    taxa_map_df$primaryAccession,
    taxa_map_df$start,
    taxa_map_df$stop,
    sep = "."
  )

  if (!all(blout_df$sseqid %in% taxa_map_df$sseqid)) {
    stop("Not all sseqid in TAXA_MAP file")
  }

  if (!all(blout_df$qseqid %in% seq_map_df$qseqid)) {
    stop("Not all qseqid in SEQ_MAP file")
  }

  blout_taxa_mapped <- left_join(x = blout_df, y = taxa_map_df, by = "sseqid")
  seq_taxa_mapped <- left_join(
    x = seq_map_df, y = blout_taxa_mapped, by = "qseqid"
  ) |>
    select(asv, qseqid, sseqid, path, organism_name, pident)

  asv_table_annot <- right_join(x = seq_taxa_mapped, y = asv_tdf, by = "asv")

}

###############################################################################
### 9. Save asv annot table
###############################################################################

write.csv(x = asv_table_annot, file = output_asv_table)
print("Output ASV table saved")

###############################################################################
### 10. Save session info
###############################################################################

output_dir <- dirname(output_asv_table)
filename_session_info <- file.path(output_dir, "session_info_taxa_annot.txt")

sink(filename_session_info)
on.exit(sink(), add = TRUE)

sep_major <- paste0(strrep("=", 80), "\n")
sep_minor <- paste0(strrep("-", 80), "\n")

cat(sep_major)
cat("Taxa Annotation Pipeline - Session Information\n")
cat(sep_major)
cat(sprintf("Date: %s\n\n", Sys.time()))

cat(sep_minor)
cat("INPUT FILES\n")
cat(sep_minor)
cat(sprintf("Input ASV table: %s\n", input_asv_table))
if (!is.null(input_fasta) && !tmp_fasta) {
  cat(sprintf("Input FASTA file: %s\n", input_fasta))
} else {
  cat("Input FASTA file: Generated from ASV table\n")
}
cat("\n")

cat(sep_minor)
cat("PARAMETERS USED\n")
cat(sep_minor)
cat(sprintf("Output ASV table: %s\n", output_asv_table))
cat(sprintf("Annotation method: %s\n", method))
cat(sprintf("Number of threads (nslots): %d\n", nslots))

if (method %in% c("NBC", "NBCandEM")) {
  cat(sprintf("Training database: %s\n", train_db))
}

if (method == "NBCandEM") {
  cat(sprintf("Reference database: %s\n", ref_db))
}

if (method == "BLAST") {
  cat(sprintf("BLAST database: %s\n", blast_db))
  cat(sprintf("Taxonomy map: %s\n", taxa_map))
  cat(sprintf("E-value threshold: %s\n", evalue))
  cat(sprintf("Minimum identity: %s%%\n", min_identity))
}

cat(sprintf("Workspace saved: %s\n\n", save_workspace))

cat(sep_minor)
cat("PACKAGE VERSIONS\n")
cat(sep_minor)
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("dada2: %s\n", packageVersion("dada2")))
cat(sprintf("tidyverse: %s\n\n", packageVersion("tidyverse")))

cat(sep_minor)
cat("FULL SESSION INFO\n")
cat(sep_minor)
print(sessionInfo())

cat("\n")
cat(sep_major)
cat("END OF SESSION INFO\n")
cat(sep_major)

sink()

cat(sprintf("Session info saved to: %s\n", filename_session_info))

###############################################################################
### 11. Save R workspace
###############################################################################

if (save_workspace) {
  filename_r_data <- file.path(output_dir, ".RData")
  save.image(file = filename_r_data)
  cat(sprintf("R workspace saved to: %s\n", filename_r_data))
}

###############################################################################
### 12. Cleanup temporary files
###############################################################################

if (tmp_fasta && file.exists(input_fasta)) {
  file.remove(input_fasta)
  cat("Temporary FASTA file removed\n")
}

if (tmp_blout && file.exists(blout)) {
  file.remove(blout)
  cat("Temporary BLAST output file removed\n")
}

if (file.exists(seq_map)) {
  file.remove(seq_map)
  cat("Temporary sequence map file removed\n")
}

cat("Taxa annotation completed successfully!\n")
