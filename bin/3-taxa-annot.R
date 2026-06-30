#!/usr/bin/env Rscript

# This code is a modified version of the DADA2 tutorial:
# https://benjjneb.github.io/dada2/tutorial.html

###############################################################################
### 1. Def. env
###############################################################################

suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(tidyverse)))

# Load toolbox.R (logging + utility functions) from the script's own directory.
# this.path::this.dir() resolves the running script's location across Rscript,
# source(), and RStudio, so toolbox.R is found regardless of the working dir.
suppressMessages(suppressWarnings(library(this.path)))
source(file.path(this.dir(), "toolbox.R"))

script_name <- "3-taxa-annot.R"
script_desc <- "Taxonomic annotation of a sequence-keyed count table (DADA2 NBC / NBC+EM)." # nolint

###############################################################################
### 2. Parse command line arguments
###############################################################################

show_usage <- function() {
  cat("Usage: ./3-taxa-annot.R <options>\n")
  cat("--help                          print this help\n")
  cat("--input_asv_table CHAR          sequence-keyed count table from DADA2/VSEARCH (required)\n") # nolint
  cat("--table_delim CHAR              delimiter for input table (default: csv)\n") # nolint
  cat("--output_dir CHAR               directory to output generated data (required)\n") # nolint
  cat("--method CHAR                   annotation method: NBC, NBCandEM (default: NBC)\n") # nolint
  cat("                                NBC: Naive Bayes Classifier; EM: Exact Matching\n") # nolint
  cat("--train_db CHAR                 training database to run NBC (default: silva_nr99_v138.1_train_set.fa.gz)\n") # nolint
  cat("--ref_db CHAR                   reference database to run EM (default: silva_species_assignment_v138.1.fa.gz)\n") # nolint
  cat("--nslots NUM                    number of threads used (default: 12)\n")
  cat("--save_workspace                save R workspace image (default: TRUE)\n") # nolint
  cat("--no_save_workspace             disable saving workspace\n")
  cat("--overwrite                     overwrite previous directory (default: FALSE)\n") # nolint
  quit(status = 0)
}

input_asv_table <- NULL
table_delim <- "csv"
output_dir <- NULL
method <- "NBC"
train_db <- "silva_nr99_v138.1_train_set.fa.gz"
ref_db <- "silva_species_assignment_v138.1.fa.gz"
nslots <- 12
save_workspace <- TRUE
overwrite <- FALSE

# Dev only — comment out before production use
# input_asv_table <- "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output_nf/2.1-dada2-piepeline-out/output/tables/asv_table.csv" # nolint
# input_asv_table <- "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output_nf/2.2.2-vsearch-pipeline-out/output/otu_table.tsv" # nolint
# output_dir <- "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/3-taxa_annot_output/" # nolint
# method <- "NBCandEM" # nolint

args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  arg <- args[i]

  if (arg == "--help" || arg == "-h") {
    show_usage()
  } else if (arg == "--input_asv_table") {
    input_asv_table <- args[i + 1]
    i <- i + 1
  } else if (arg == "--table_delim") {
    table_delim <- args[i + 1]
    i <- i + 1
  } else if (arg == "--output_dir") {
    output_dir <- args[i + 1]
    i <- i + 1
  } else if (arg == "--method") {
    method <- args[i + 1]
    i <- i + 1
  } else if (arg == "--train_db") {
    train_db <- args[i + 1]
    i <- i + 1
  } else if (arg == "--ref_db") {
    ref_db <- args[i + 1]
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

if (is.null(output_dir)) {
  cat("Error: --output_dir is required\n")
  show_usage()
}

if (!file.exists(input_asv_table)) {
  cat(sprintf("Error: Input ASV table '%s' does not exist\n", input_asv_table))
  quit(status = 1)
}

if (!method %in% c("NBC", "NBCandEM")) {
  cat(sprintf(
    "Error: Invalid method '%s'. Must be one of: NBC, NBCandEM\n",
    method
  ))
  quit(status = 1)
}

if (dir.exists(output_dir)) {
  if (overwrite) {
    cat(sprintf("Removing existing output directory: %s\n", output_dir))
    unlink(output_dir, recursive = TRUE)
  } else {
    cat(sprintf(
      "Error: Output directory '%s' already exists. Use --overwrite to overwrite\n", # nolint
      output_dir
    ))
    quit(status = 1)
  }
}

###############################################################################
### 4. Create output dirs
###############################################################################

results_dir <- file.path(output_dir, "output")
logs_dir    <- file.path(output_dir, "logs")
stats_dir   <- file.path(output_dir, "stats")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, showWarnings = FALSE)
dir.create(stats_dir, showWarnings = FALSE)
dir.create(file.path(results_dir, "tables"), recursive = TRUE, showWarnings = FALSE) # nolint

###############################################################################
### 5. Load and format data
###############################################################################

log_msg("Loading and formatting ASV table ...")

if (table_delim == "csv") {
  asv_data <- read_csv(
    file = input_asv_table, col_names = TRUE, show_col_types = FALSE
  )
} else if (table_delim == "tsv") {
  asv_data <- read_tsv(
    file = input_asv_table, col_names = TRUE, show_col_types = FALSE
  )
} else {
  log_error(sprintf("Unsupported table delimiter: '%s'", table_delim))
  quit(status = 1)
}
# The first column holds the sequences (DADA2-style sequence-keyed table). Its
# header is blank/auto-named, so normalize it by position before use.
colnames(asv_data)[1] <- "asv"

asv_matrix <- asv_data |>
  column_to_rownames("asv") |>
  as.matrix() |>
  t()

asv_tdf <- asv_matrix |>
  t() |>
  as.data.frame() |>
  rownames_to_column("asv")

sample_names <- colnames(asv_tdf)[-1]
log_msg(sprintf(
  "Loaded %d sequences across %d samples", nrow(asv_tdf), length(sample_names)
))

###############################################################################
### 6. Run taxa annot: Naive Bayes Classifier (NBC)
###############################################################################

if (method == "NBC") {

  train_db <- ensure_database(train_db)

  log_msg("Running NBC ...")
  taxa <- assignTaxonomy(
    seqs = asv_matrix,
    refFasta = train_db,
    outputBootstraps = TRUE,
    multithread = nslots
  )

log_msg("Running NBC worked")

  tax_mat <- taxa$tax

  taxa_df <- taxa |>
    as.data.frame() |>
    rownames_to_column("asv")

  asv_table_annot <- right_join(x = taxa_df, y = asv_tdf, by = "asv")

}

###############################################################################
### 7. Run taxa annot: NBC + Exact Matching (EM)
###############################################################################

if (method == "NBCandEM") {

  log_msg("Running NBC ...")
  train_db <- ensure_database(train_db)
  taxa <- assignTaxonomy(
    seqs = asv_matrix,
    refFasta = train_db,
    minBoot = 50,
    outputBootstraps = TRUE,
    multithread = nslots
  )

  log_msg("Running EM (addSpecies) ...")
  ref_db <- ensure_database(ref_db)
  specs <- addSpecies(taxtab = taxa$tax, refFasta = ref_db)

  tax_mat <- specs

  specs_df <- as.data.frame(specs) 
  colnames(specs_df) <- paste("tax", colnames(specs_df), sep = ".")
  boot_df <- as.data.frame(taxa$boot) 
  colnames(boot_df) <- paste("boot", colnames(boot_df), sep = ".")

  if (!all(rownames(specs_df) == rownames(boot_df))) {
    log_error("Row names of taxonomy and bootstrap data frames do not match")
    quit(status = 1)
  }
  taxa_df <- cbind(specs_df, boot_df) |>
    rownames_to_column("asv")

  asv_table_annot <- right_join(x = taxa_df, y = asv_tdf, by = "asv")

}

###############################################################################
### 8. Save asv annot table
###############################################################################

filename_annot <- file.path(results_dir, "tables", "asv_table_annot.csv")
write.csv(x = asv_table_annot, file = filename_annot, row.names = FALSE)
log_msg(sprintf("Annotated ASV table saved to: %s", filename_annot))

###############################################################################
### 9. Compute annotation statistics
###############################################################################

log_msg("Computing annotation statistics ...")

# Mean/SD of bootstrap support over the ASVs classified at a given rank.
# Return NA when there are too few classified ASVs for the statistic to be
# defined (mean needs >=1 value, sd needs >=2), so the stats file never carries
# NaN (mean of an empty vector) or an accidental NA (sd of a single value).
boot_mean <- function(boot, classified) {
  vals <- boot[classified]
  if (length(vals) == 0) NA else mean(vals, na.rm = TRUE)
}
boot_sd <- function(boot, classified) {
  vals <- boot[classified]
  if (length(vals) < 2) NA else sd(vals, na.rm = TRUE)
}

ranks <- c("Phylum", "Class", "Order", "Family", "Genus")

# Build one stats row over the subset of ASVs selected by `mask` (a logical
# vector over rows of asv_table_annot): number of ASVs, mean/sd bootstrap
# support per rank, and the percent annotated to species. Used once per sample
# (ASVs present in that sample) and once for the pooled all_samples total.
annot_stats_row <- function(sample, mask) {
  row <- data.frame(
    sample = sample, n_asvs = sum(mask), stringsAsFactors = FALSE
  )
  for (rank in ranks) {
    boot <- asv_table_annot[[paste0("boot.", rank)]]
    classified <- mask & !is.na(asv_table_annot[[paste0("tax.", rank)]])
    row[[paste0("mean_", tolower(rank), "_boot")]] <- as.numeric(boot_mean(boot, classified)) # nolint
    row[[paste0("sd_",   tolower(rank), "_boot")]] <- as.numeric(boot_sd(boot, classified))   # nolint
  }
  # Species annotation only comes from NBCandEM (addSpecies); NA under NBC.
  if (!is.null(asv_table_annot$tax.Species) && sum(mask) > 0) {
    row$perc_spec_annot <-
      sum(mask & !is.na(asv_table_annot$tax.Species)) / sum(mask) * 100
  } else {
    row$perc_spec_annot <- NA_real_
  }
  row
}

# One row per sample (ASVs with count > 0 in that sample), plus a pooled
# all_samples row over every ASV.
tax_stats <- do.call(rbind, c(
  lapply(sample_names, function(s) annot_stats_row(s, asv_table_annot[[s]] > 0)), # nolint
  list(annot_stats_row("all_samples", rep(TRUE, nrow(asv_table_annot))))
))

###############################################################################
### 10. Export stats
###############################################################################

filename_stats <- file.path(stats_dir,
                            paste0(sub(".R", "", script_name), "-stats.tsv")) # nolint
write.table(file = filename_stats, tax_stats,
            sep = "\t", row.names = FALSE, quote = FALSE)

log_msg(sprintf("Stats table saved to: %s", filename_stats))

###############################################################################
### 11. Save R workspace
###############################################################################

if (save_workspace) {

  filename_r_data <- file.path(results_dir, ".RData")
  save.image(file = filename_r_data)

}

###############################################################################
### 12. Write standardized log (general info + run log)
###############################################################################

log_msg("\033[0;32m3-taxa-annot.R completed successfully\033[0m")

cmd_executed <- paste(script_name, paste(args, collapse = " "))
filename_log <- file.path(logs_dir, paste0(sub(".R", "", script_name), ".log"))

log_text <- build_log(
  script_name = script_name,
  script_desc = script_desc,
  sample_name = paste(sample_names, collapse = ", "),
  inputs = c(
    sprintf("Input ASV table: %s", input_asv_table),
    sprintf("Sequences: %d", nrow(asv_tdf)),
    sprintf("Samples: %d", length(sample_names))
  ),
  params = c(
    sprintf("Method: %s", method),
    sprintf("Training database: %s", train_db),
    if (method == "NBCandEM") sprintf("Reference database: %s", ref_db) else NULL,
    sprintf("Threads: %d", nslots)
  ),
  outputs = c(
    sprintf("Annotated ASV table: %s", filename_annot),
    sprintf("Results directory: %s", results_dir),
    sprintf("Statistics: %s", filename_stats)
  ),
  command = cmd_executed,
  exit_status = 0,
  # R session info recorded at the end of the log (third-party tools section)
  tool_log = capture.output(sessionInfo())
)

writeLines(log_text, filename_log)
