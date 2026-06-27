#!/usr/bin/env Rscript

# This code is a modified version of the DADA2 tutorial:
# https://benjjneb.github.io/dada2/tutorial.html

###############################################################################
### 1. Def. env
###############################################################################

suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ShortRead)))

# Run-log buffer: log_msg() prints to console and accumulates lines for the log file # nolintr
log_lines <- character(0)
log_msg <- function(msg) {
  line <- paste0("[INFO] ", msg)
  message(line)
  log_lines[[length(log_lines) + 1]] <<- line
}

###############################################################################
### 2. Parse command line arguments
###############################################################################

show_usage <- function() {
  cat("Usage: ./dada2_pipeline.R <options>\n")
  cat("--help                          print this help\n")
  cat("--input_dir CHAR                directory with the input raw fastq files (required)\n") # nolintr
  cat("--output_dir CHAR               directory to output generated data (required)\n") # nolintr
  cat("--nslots NUM                    number of threads used (default: 12)\n") # nolintr
  cat("--trunc_r1 NUM                  number of nuc to remove in R1 from the 3' end (default: 250)\n") # nolintr
  cat("--trunc_r2 NUM                  number of nuc to remove in R2 from the 3' end (default: 200)\n") # nolintr
  cat("--pattern_r1 CHAR               pattern of R1 reads to load fastq files (default: _L001_R1_001.fastq.gz)\n") # nolintr
  cat("--pattern_r2 CHAR               pattern of R2 reads to load fastq files (default: _L001_R2_001.fastq.gz)\n") # nolintr
  cat("--bimeras_method CHAR           method to check bimeras: pooled, consensus, per-sample (default: consensus)\n") # nolintr
  cat("--min_overlap NUM               minimum number of nucleotides to overlap in merging (default: 12)\n") # nolintr
  cat("--pooled                        use pooled option when running dada2 (default: TRUE)\n") # nolintr
  cat("--no_pooled                     disable pooled option\n")
  cat("--qual_plot                     create quality plots (default: TRUE)\n")
  cat("--no_qual_plot                  disable quality plots\n")
  cat("--err_plot                      create error plots (default: TRUE)\n")
  cat("--no_err_plot                   disable error plots\n")
  cat("--save_workspace                save R workspace image (default: TRUE)\n") # nolintr
  cat("--no_save_workspace             disable saving workspace\n")
  cat("--overwrite                     overwrite previous directory (default: FALSE)\n") # nolintr
  quit(status = 0)
}

input_dir <- NULL
output_dir <- NULL
nslots <- 12
trunc_r1 <- 250
trunc_r2 <- 200
pattern_r1 <- "_L001_R1_001.fastq.gz"
pattern_r2 <- "_L001_R2_001.fastq.gz"
min_overlap <- 12
bimeras_method <- "consensus"
pool_option <- TRUE
qual_plot <- TRUE
err_plot <- TRUE
save_workspace <- TRUE
overwrite <- FALSE

# Dev only — comment out before production use
input_dir <- "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/03_primer_removal/trimmed" # nolintr
output_dir <- "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/2.1-dada2_pipeline_output/" # nolintr
pattern_r1 <- "_L001_R1_trimmed.fastq.gz" # nolintr
pattern_r2 <- "_L001_R2_trimmed.fastq.gz" # nolintr

args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  arg <- args[i]

  if (arg == "--help" || arg == "-h") {
    show_usage()
  } else if (arg == "--input_dir") {
    input_dir <- args[i + 1]
    i <- i + 1
  } else if (arg == "--output_dir") {
    output_dir <- args[i + 1]
    i <- i + 1
  } else if (arg == "--nslots") {
    nslots <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--trunc_r1") {
    trunc_r1 <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--trunc_r2") {
    trunc_r2 <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--pattern_r1") {
    pattern_r1 <- args[i + 1]
    i <- i + 1
  } else if (arg == "--pattern_r2") {
    pattern_r2 <- args[i + 1]
    i <- i + 1
  } else if (arg == "--min_overlap") {
    min_overlap <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--bimeras_method") {
    bimeras_method <- args[i + 1]
    i <- i + 1
  } else if (arg == "--pooled") {
    pool_option <- TRUE
  } else if (arg == "--no_pooled") {
    pool_option <- FALSE
  } else if (arg == "--qual_plot") {
    qual_plot <- TRUE
  } else if (arg == "--no_qual_plot") {
    qual_plot <- FALSE
  } else if (arg == "--err_plot") {
    err_plot <- TRUE
  } else if (arg == "--no_err_plot") {
    err_plot <- FALSE
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

if (is.null(input_dir)) {
  cat("Error: --input_dir is required\n")
  show_usage()
}

if (is.null(output_dir)) {
  cat("Error: --output_dir is required\n")
  show_usage()
}

if (!dir.exists(input_dir)) {
  cat(sprintf("Error: Input directory '%s' does not exist\n", input_dir))
  quit(status = 1)
}

if (dir.exists(output_dir)) {
  if (overwrite) {
    cat(sprintf("Removing existing output directory: %s\n", output_dir))
    unlink(output_dir, recursive = TRUE)
  } else {
    cat(sprintf(
      "Error: Output directory '%s' already exists. Use --overwrite to overwrite\n", # nolintr
      output_dir
    ))
    quit(status = 1)
  }
}

###############################################################################
### 3. Create output dirs
###############################################################################

results_dir <- file.path(output_dir, "output")
logs_dir    <- file.path(output_dir, "logs")
stats_dir   <- file.path(output_dir, "stats")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, showWarnings = FALSE)
dir.create(stats_dir, showWarnings = FALSE)
dir.create(file.path(results_dir, "plots"), recursive = TRUE, showWarnings = FALSE) # nolintr
dir.create(file.path(results_dir, "filtered"), showWarnings = FALSE)
dir.create(file.path(results_dir, "tables"), showWarnings = FALSE)

###############################################################################
### 4. Load data
###############################################################################

log_msg("Loading data ...")
raw_r1 <- sort(list.files(input_dir, pattern = pattern_r1, full.names = TRUE))
raw_r2 <- sort(list.files(input_dir, pattern = pattern_r2, full.names = TRUE))

sample_names <- basename(raw_r1) |>
  sub(pattern = pattern_r1, replacement = "")

###############################################################################
### 5. Create quality data plots
###############################################################################

if (qual_plot) {

  log_msg("Creating quality plots ...")
  pq_r1 <- plotQualityProfile(raw_r1)
  pq_r2 <- plotQualityProfile(raw_r2)

  filename_qual_r1 <- file.path(results_dir, "plots/quality_plot_R1.pdf")
  filename_qual_r2 <- file.path(results_dir, "plots/quality_plot_R2.pdf")

  ggsave(pq_r1,
         device = "pdf",
         filename = filename_qual_r1,
         width = 20,
         height = 20)
  ggsave(pq_r2,
         device = "pdf",
         filename = filename_qual_r2,
         width = 20,
         height = 20)
}

###############################################################################
### 6. Quality check
###############################################################################

log_msg("Quality check ...")
filt_r1 <- file.path(
  results_dir, "filtered", paste(sample_names, "R1_filt.fastq.gz", sep = "_")
)
filt_r2 <- file.path(
  results_dir, "filtered", paste(sample_names, "R2_filt.fastq.gz", sep = "_")
)

names(filt_r1) <- sample_names
names(filt_r2) <- sample_names

filter_and_trim_log <- filterAndTrim(
  fwd = raw_r1,
  filt = filt_r1,
  rev = raw_r2,
  filt.rev = filt_r2,
  truncLen = c(trunc_r1, trunc_r2),
  maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
  compress = TRUE,
  multithread = nslots
)

###############################################################################
### 7. Learn error rates
###############################################################################

log_msg("Learning errors ...")
err_r1 <- learnErrors(filt_r1, multithread = nslots)
err_r2 <- learnErrors(filt_r2, multithread = nslots)

###############################################################################
### 8. Create error plots
###############################################################################

if (err_plot) {

  log_msg("Creating error plots ...")
  pe_r1 <- plotErrors(err_r1, nominalQ = TRUE)
  pe_r2 <- plotErrors(err_r2, nominalQ = TRUE)

  filename_error_r1 <- file.path(results_dir, "plots/error_plot_R1.pdf")
  filename_error_r2 <- file.path(results_dir, "plots/error_plot_R2.pdf")

  ggsave(pe_r1,
         device = "pdf",
         filename = filename_error_r1,
         width = 10, height = 10)
  ggsave(pe_r2,
         device = "pdf",
         filename = filename_error_r2,
         width = 10,
         height = 10)

}

###############################################################################
### 9. Dereplicate
###############################################################################

# This step is deliberaly skept, since dereplication is done internally by the dada() function. # nolintr
# derep_r1 <- derepFastq(filt_r1, verbose = TRUE) # nolintr
# derep_r2 <- derepFastq(filt_r2, verbose = TRUE) # nolintr

###############################################################################
### 10. Apply sample inference algorithms
###############################################################################

log_msg("Finding ASVs ...")
dada_r1 <- dada(filt_r1, err = err_r1, multithread = nslots,
                pool = pool_option)
dada_r2 <- dada(filt_r2, err = err_r2, multithread = nslots,
                pool = pool_option)

###############################################################################
### 11. Merge paired reads
###############################################################################

log_msg("Merging ...")
mergers <- mergePairs(
  dada_r1, filt_r1,
  dada_r2, filt_r2,
  minOverlap = min_overlap,
  verbose = TRUE
)
# The output is a list of data.frames from each sample.
# Each data.frame contains the merged $sequence, its $abundance,
# and the indices of the $forward and $reverse sequence variants that were
# merged.

###############################################################################
### 12. Construct sequence table
###############################################################################

seqtab <- makeSequenceTable(mergers)
x <- dim(seqtab)
log_msg(paste("asv table dim:", x[1], "x", x[2], sep = " "))

###############################################################################
### 13. Remove chimeras
###############################################################################

log_msg("Bimeras check ...")
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = bimeras_method,
  multithread = nslots,
  verbose = TRUE
)

x <- dim(seqtab_nochim)
log_msg(paste("asv table (no bimeras) dim:", x[1], "x", x[2], sep = " "))

perc_bim <- (1 - sum(seqtab_nochim) / sum(seqtab)) * 100
log_msg(paste("bimeras: ", round(perc_bim, 4), "%", sep = ""))

###############################################################################
### 14. Save ASV table
###############################################################################

filename_asv <- file.path(results_dir, "asv_table.csv")

write.csv(x = t(seqtab_nochim), file = filename_asv)

###############################################################################
### 15. Track number of seqs
###############################################################################

log_msg("Creating n seq and length plots ...")
count_seqs <- function(x) {
  sum(getUniques(x))
}

track_n_seqs <- data.frame(
  samples = sample_names,
  raw = filter_and_trim_log[, 1],
  filtered = filter_and_trim_log[, 2],
  denoisedR1 = sapply(dada_r1, count_seqs),
  denoisedR2 = sapply(dada_r2, count_seqs),
  merged = sapply(mergers, count_seqs),
  nobim = rowSums(seqtab_nochim),
  stringsAsFactors = FALSE
)

col_totals <- colSums(select(track_n_seqs, -samples))
track_n_seqs["all", ] <- data.frame(
  samples = "all_samples",
  t(col_totals),
  stringsAsFactors = FALSE
)

track_n_seqs_long <- track_n_seqs |>
  gather(key = "var", value = "value", raw:nobim)

track_n_seqs_long$var <- factor(
  track_n_seqs_long$var,
  levels = c("raw", "filtered", "denoisedR1", "denoisedR2", "merged", "nobim")
)

###############################################################################
### 16. Track reads length
###############################################################################

track_mergers_length_long <- lapply(mergers, "[[", 1) |>
  lapply(nchar) |>
  plyr::ldply(cbind)

colnames(track_mergers_length_long) <- c("samples", "length")

###############################################################################
### 17. Create plots
###############################################################################

nseq_barplots <- ggplot(track_n_seqs_long, aes(x = var, y = value)) +
  facet_wrap(~samples, ncol = 7, scales = "free") +
  geom_bar(stat = "identity", fill = "gray50", alpha = 0.7) +
  ylab("Number of seqs") +
  theme_light() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )

seqlength_barplots <- ggplot(track_mergers_length_long, aes(y = length)) +
  facet_wrap(~samples, ncol = 5, scales = "free") +
  geom_boxplot(fill = "gray80", alpha = 0.7) +
  ylab("Read length") +
  theme_light() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", face = "bold"),
    axis.title.x = element_blank()
  )

###############################################################################
### 18. Save plots
###############################################################################

filename_nseq <- file.path(results_dir, "plots/nseq_barplot.pdf")
filename_seqlength <- file.path(results_dir, "plots/seq_length_hist.pdf")

ggsave(nseq_barplots, filename = filename_nseq,
  device = "pdf", width = 18, height = 30)

ggsave(seqlength_barplots, filename = filename_seqlength,
  device = "pdf", width = 10, height = 30)

###############################################################################
### 19. Save track stats tables
###############################################################################

filename_nseq <- file.path(results_dir, "tables/nseq_counts.csv")
write.csv(file = filename_nseq, track_n_seqs)

track_mergers_length_stats <- track_mergers_length_long |>
  group_by(samples) |>
  summarize(
    mean = mean(length),
    sd = sd(length),
    max = max(length),
    min = min(length)
  )

filename_seqlength_stats <- file.path(results_dir, "tables/seq_length_stats.csv")
write.csv(file = filename_seqlength_stats, track_mergers_length_stats)

###############################################################################
### 19b. Save standardized stats table (read tracking through the pipeline)
###############################################################################

filename_stats <- file.path(stats_dir, "dada2_stats.tsv")
write.table(track_n_seqs, file = filename_stats, sep = "\t",
  row.names = FALSE, quote = FALSE)
log_msg(sprintf("Stats table saved to: %s", filename_stats))

###############################################################################
### 20. Save session info
###############################################################################

filename_session_info <- file.path(results_dir, "session_info.txt")

sink(filename_session_info)

sep_major <- paste0(strrep("=", 80), "\n")
sep_minor <- paste0(strrep("-", 80), "\n")

cat(sep_major)
cat("DADA2 Pipeline - Session Information\n")
cat(sep_major)
cat(sprintf("Date: %s\n\n", Sys.time()))

cat(sep_minor)
cat("INPUT FILES\n")
cat(sep_minor)
cat(sprintf("Input directory: %s\n", input_dir))
cat(sprintf("Number of sample pairs: %d\n", length(raw_r1)))
cat(sprintf("Pattern R1: %s\n", pattern_r1))
cat(sprintf("Pattern R2: %s\n\n", pattern_r2))

cat("Sample files:\n")
for (i in seq_along(sample_names)) {
  cat(sprintf("  %d. %s\n", i, sample_names[i]))
}
cat("\n")

cat(sep_minor)
cat("PARAMETERS USED\n")
cat(sep_minor)
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("Number of threads (nslots): %d\n", nslots))
cat(sprintf("Truncation length R1: %d\n", trunc_r1))
cat(sprintf("Truncation length R2: %d\n", trunc_r2))
cat(sprintf("Minimum overlap for merging: %d\n", min_overlap))
cat(sprintf("Bimeras detection method: %s\n", bimeras_method))
cat(sprintf("Pooled option: %s\n", pool_option))
cat(sprintf("Quality plots generated: %s\n", qual_plot))
cat(sprintf("Error plots generated: %s\n", err_plot))
cat(sprintf("Workspace saved: %s\n\n", save_workspace))

cat(sep_minor)
cat("PACKAGE VERSIONS\n")
cat(sep_minor)
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("dada2: %s\n", packageVersion("dada2")))
cat(sprintf("tidyverse: %s\n", packageVersion("tidyverse")))
cat(sprintf("ShortRead: %s\n\n", packageVersion("ShortRead")))

cat(sep_minor)
cat("FULL SESSION INFO\n")
cat(sep_minor)
print(sessionInfo())

cat("\n")
cat(sep_major)
cat("END OF SESSION INFO\n")
cat(sep_major)

sink()

log_msg(sprintf("Session info saved to: %s", filename_session_info))

###############################################################################
### 21. Save R workspace
###############################################################################

if (save_workspace) {

  filename_r_data <- file.path(results_dir, ".RData")
  save.image(file = filename_r_data)

}

###############################################################################
### 22. Summary report
###############################################################################

log_msg("Generating summary report...")
filename_report <- file.path(output_dir, "summary_report.txt")

all_row <- track_n_seqs[track_n_seqs$samples == "all_samples", ]

sink(filename_report)
cat(sep_major)
cat("DADA2 Pipeline - Summary Report\n")
cat(sep_major)
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Input directory:  %s\n", input_dir))
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("Samples processed: %d\n\n", length(raw_r1)))

cat("Parameters:\n")
cat(sep_minor)
cat(sprintf("  Threads:             %d\n", nslots))
cat(sprintf("  Truncation R1:       %d\n", trunc_r1))
cat(sprintf("  Truncation R2:       %d\n", trunc_r2))
cat(sprintf("  Min merge overlap:   %d\n", min_overlap))
cat(sprintf("  Bimeras method:      %s\n", bimeras_method))
cat(sprintf("  Pooled:              %s\n", pool_option))
cat(sprintf("  Bimeras removed:     %.2f%%\n\n", perc_bim))

cat("Output files:\n")
cat(sep_minor)
cat(sprintf("  ASV table:  %s\n", filename_asv))
cat(sprintf("  Results:    %s\n", results_dir))
cat(sprintf("  Stats:      %s\n", filename_stats))
cat(sprintf("  Logs:       %s\n\n", logs_dir))

cat(sep_major)
cat("\nRead tracking (all samples pooled):\n\n")
cat(sprintf("  raw:        %s\n", all_row$raw))
cat(sprintf("  filtered:   %s\n", all_row$filtered))
cat(sprintf("  denoisedR1: %s\n", all_row$denoisedR1))
cat(sprintf("  denoisedR2: %s\n", all_row$denoisedR2))
cat(sprintf("  merged:     %s\n", all_row$merged))
cat(sprintf("  nobim:      %s\n", all_row$nobim))
sink()

###############################################################################
### 23. Write run log
###############################################################################

log_msg("\033[0;32m2.1-dada2-pipeline.R completed successfully\033[0m")
filename_log <- file.path(logs_dir, "dada2.log")
ansi_re <- "\033\\[[0-9;]*m"
writeLines(gsub(ansi_re, "", log_lines), filename_log)
