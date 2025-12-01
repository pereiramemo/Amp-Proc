#!/usr/bin/env Rscript

# This code is a modified version of the DADA2 tutorial: 
# https://benjjneb.github.io/dada2/tutorial.html

###############################################################################
### 1. Def. env
###############################################################################

suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ShortRead)))

###############################################################################
### 2. Parse command line arguments
###############################################################################

# Function to display help
show_usage <- function() {
  cat("Usage: ./dada2_pipeline.R <options>\n")
  cat("--help                          print this help\n")
  cat("--input_dir CHAR                directory with the input raw fastq files (required)\n")
  cat("--output_dir CHAR               directory to output generated data (required)\n")
  cat("--nslots NUM                    number of threads used (default: 12)\n")
  cat("--trunc_r1 NUM                  number of nuc to remove in R1 from the 3' end (default: 250)\n")
  cat("--trunc_r2 NUM                  number of nuc to remove in R2 from the 3' end (default: 200)\n")
  cat("--pattern_r1 CHAR               pattern of R1 reads to load fastq files (default: _L001_R1_001.fastq.gz)\n")
  cat("--pattern_r2 CHAR               pattern of R2 reads to load fastq files (default: _L001_R2_001.fastq.gz)\n")
  cat("--bimeras_method CHAR           method to check bimeras: pooled, consensus, per-sample (default: consensus)\n")
  cat("--min_overlap NUM               minimum number of nucleotides to overlap in merging (default: 12)\n")
  cat("--pooled                        use pooled option when running dada2 (default: TRUE)\n")
  cat("--no_pooled                     disable pooled option\n")
  cat("--qual_plot                     create quality plots (default: TRUE)\n")
  cat("--no_qual_plot                  disable quality plots\n")
  cat("--err_plot                      create error plots (default: TRUE)\n")
  cat("--no_err_plot                   disable error plots\n")
  cat("--save_workspace                save R workspace image (default: TRUE)\n")
  cat("--no_save_workspace             disable saving workspace\n")
  cat("--overwrite                     overwrite previous directory (default: FALSE)\n")
  quit(status = 0)
}

# Initialize parameters with defaults
INPUT_DIR <- NULL
OUTPUT_DIR <- NULL
NSLOTS <- 12
TRUNC_R1 <- 250
TRUNC_R2 <- 200
PATTERN_R1 <- "_L001_R1_001.fastq.gz"
PATTERN_R2 <- "_L001_R2_001.fastq.gz"
MIN_OVERLAP <- 12
BIMERAS_METHOD <- "consensus"
POOL_OPTION <- TRUE
QUAL_PLOT <- TRUE
ERR_PLOT <- TRUE
SAVE_WORKSPACE <- TRUE
OVERWRITE <- FALSE

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  arg <- args[i]
  
  if (arg == "--help" || arg == "-h") {
    show_usage()
  } else if (arg == "--input_dir") {
    INPUT_DIR <- args[i + 1]
    i <- i + 1
  } else if (arg == "--output_dir") {
    OUTPUT_DIR <- args[i + 1]
    i <- i + 1
  } else if (arg == "--nslots") {
    NSLOTS <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--trunc_r1") {
    TRUNC_R1 <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--trunc_r2") {
    TRUNC_R2 <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--pattern_r1") {
    PATTERN_R1 <- args[i + 1]
    i <- i + 1
  } else if (arg == "--pattern_r2") {
    PATTERN_R2 <- args[i + 1]
    i <- i + 1
  } else if (arg == "--min_overlap") {
    MIN_OVERLAP <- as.numeric(args[i + 1])
    i <- i + 1
  } else if (arg == "--bimeras_method") {
    BIMERAS_METHOD <- args[i + 1]
    i <- i + 1
  } else if (arg == "--pooled") {
    POOL_OPTION <- TRUE
  } else if (arg == "--no_pooled") {
    POOL_OPTION <- FALSE
  } else if (arg == "--qual_plot") {
    QUAL_PLOT <- TRUE
  } else if (arg == "--no_qual_plot") {
    QUAL_PLOT <- FALSE
  } else if (arg == "--err_plot") {
    ERR_PLOT <- TRUE
  } else if (arg == "--no_err_plot") {
    ERR_PLOT <- FALSE
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

# Validate required parameters
if (is.null(INPUT_DIR)) {
  cat("Error: --input_dir is required\n")
  show_usage()
}

if (is.null(OUTPUT_DIR)) {
  cat("Error: --output_dir is required\n")
  show_usage()
}

# Check input directory exists
if (!dir.exists(INPUT_DIR)) {
  cat(sprintf("Error: Input directory '%s' does not exist\n", INPUT_DIR))
  quit(status = 1)
}

# Handle output directory
if (dir.exists(OUTPUT_DIR)) {
  if (OVERWRITE) {
    cat(sprintf("Removing existing output directory: %s\n", OUTPUT_DIR))
    unlink(OUTPUT_DIR, recursive = TRUE)
  } else {
    cat(sprintf("Error: Output directory '%s' already exists. Use --overwrite to overwrite\n", OUTPUT_DIR))
    quit(status = 1)
  }
}

###############################################################################
### 3. Create output dirs
###############################################################################

dir.create(OUTPUT_DIR)
file.path(OUTPUT_DIR,"plots") %>% dir.create 
file.path(OUTPUT_DIR,"filtered") %>% dir.create
file.path(OUTPUT_DIR,"tables") %>% dir.create

###############################################################################
### 4. Load data
###############################################################################

print("Loading data ...")
rawR1 <- sort(list.files(INPUT_DIR, pattern = PATTERN_R1, full.names = T))
rawR2 <- sort(list.files(INPUT_DIR, pattern = PATTERN_R2, full.names = T))

sample.names <- basename(rawR1) %>%
                sub(pattern = PATTERN_R1, replacement = "")

###############################################################################
### 5. Create quality data plots
###############################################################################

if (QUAL_PLOT == TRUE) {
  
  print("Creating quality plots ...")
  pqR1 <- plotQualityProfile(rawR1)
  pqR2 <- plotQualityProfile(rawR2)

  filename_qual_R1 <- file.path(OUTPUT_DIR,"plots/quality_plot_R1.pdf")
  filename_qual_R2 <- file.path(OUTPUT_DIR,"plots/quality_plot_R2.pdf")

  ggsave(pqR1, device = "pdf", filename = filename_qual_R1, width = 20, height = 20)
  ggsave(pqR2, device = "pdf", filename = filename_qual_R2, width = 20, height = 20)

}
 
###############################################################################
### 6. Quatlity check
###############################################################################

print("Quality check ...")
filtR1 <- file.path(OUTPUT_DIR, "filtered", paste(sample.names, "R1_filt.fastq.gz", sep = "_"))
filtR2 <- file.path(OUTPUT_DIR, "filtered", paste(sample.names, "R2_filt.fastq.gz", sep = "_"))

names(filtR1) <- sample.names
names(filtR2) <- sample.names

filterAndTrim_log <- filterAndTrim(fwd = rawR1, filt = filtR1, 
                                   rev = rawR2, filt.rev = filtR2, 
                                   truncLen = c(TRUNC_R1, TRUNC_R2),
                                   maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                                   compress=TRUE, 
                                   multithread = NSLOTS)

###############################################################################
### 7. Learn error rates
###############################################################################

print("Learning errors ...")
errR1 <- learnErrors(filtR1, multithread=NSLOTS)
errR2 <- learnErrors(filtR2, multithread=NSLOTS)

###############################################################################
### 8. Create error plots
###############################################################################

if (ERR_PLOT == TRUE) {

  print("Creating error plots ...")
  peR1 <- plotErrors(errR1, nominalQ=TRUE)
  peR2 <- plotErrors(errR2, nominalQ=TRUE)

  filename_error_R1 <- file.path(OUTPUT_DIR,"plots/error_plot_R1.pdf")
  filename_error_R2 <- file.path(OUTPUT_DIR,"plots/error_plot_R2.pdf")

  ggsave(peR1, device = "pdf", filename = filename_error_R1, width = 10, height = 10)
  ggsave(peR1, device = "pdf", filename = filename_error_R2, width = 10, height = 10)

}

###############################################################################
### 9. Derepliacate
###############################################################################

derepR1 <- derepFastq(filtR1, verbose=TRUE)
derepR2 <- derepFastq(filtR2, verbose=TRUE)
  
###############################################################################
### 10. Apply sample inference algorithms
###############################################################################

print("Finding ASVs ...")
dadaR1 <- dada(derepR1, err=errR1, multithread=NSLOTS, pool = POOL_OPTION)
dadaR2 <- dada(derepR2, err=errR2, multithread=NSLOTS, pool = POOL_OPTION)

###############################################################################
### 11. Merge paired reads
###############################################################################

print("Merging ...")
mergers <- mergePairs(dadaR1, filtR1, 
                      dadaR2, filtR2,
                      minOverlap = MIN_OVERLAP,
                      verbose=TRUE)
# Note:
# The output is a list of data.frames from each sample. 
# Each data.frame contains the merged $sequence, its $abundance, 
# and the indices of the $forward and $reverse sequence variants that were merged.

###############################################################################
### 12. Construct sequence table
###############################################################################

seqtab <- makeSequenceTable(mergers)
x <- dim(seqtab)
print(paste("asv table dim:", x[1], "x", x[2], sep = " "))

###############################################################################
### 13. Remove chimeras
###############################################################################

print("Bimeras check ...")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = BIMERAS_METHOD,
                                    multithread = NSLOTS, verbose=TRUE)

x <- dim(seqtab.nochim)
print(paste("asv table (no bimeras) dim:", x[1], "x", x[2], sep = " "))

perc_bim <- (1 - sum(seqtab.nochim)/sum(seqtab))*100
print(paste("bimeras: ",round(perc_bim,4), "%", sep =""))

###############################################################################
### 14. Save ASV table
###############################################################################

filename_asv <- paste(OUTPUT_DIR,"/asv_table.csv", sep ="")

write.csv(x = t(seqtab.nochim), file = filename_asv)

###############################################################################
### 15. Track number of seqs
###############################################################################

print("Creating n seq and length plots ...")
count_seqs <- function(x){ sum(getUniques(x)) }
                       
track_n_seqs <- data.frame(samples = sample.names,
                           raw = filterAndTrim_log[,1],
                           filtered = filterAndTrim_log[,2],
                           denoisedR1 = sapply(dadaR1, count_seqs),
                           denoisedR2 = sapply(dadaR2, count_seqs), 
                           merged = sapply(mergers, count_seqs),
                           nobim = rowSums(seqtab.nochim),
                           stringsAsFactors = F)

track_n_seqs["all", ] <- data.frame(samples = "all_samples", 
                                    colSums(track_n_seqs %>% select(-samples)) %>% t,
                                    stringsAsFactors = F)
  
track_n_seqs_long <- track_n_seqs %>%
                     gather(key = "var", value = "value", raw:nobim)

track_n_seqs_long$var <- factor(track_n_seqs_long$var, 
                                levels = c("raw","filtered", "denoisedR1", 
                                           "denoisedR2", "merged","nobim"))

###############################################################################
### 16. Track reads length
###############################################################################

track_mergers_length_long <- lapply(mergers, "[[", 1) %>%
                             lapply(., nchar) %>%
                             plyr::ldply(., cbind) 

colnames(track_mergers_length_long) <- c("samples","length")

###############################################################################
### 17. Create plots
###############################################################################

nseq_barplots <- ggplot(track_n_seqs_long, aes(x = var, y = value)) +
                 facet_wrap(~ samples, ncol = 7, scales = "free") +
                 geom_bar(stat = "identity", fill = "gray50", alpha = 0.7) +
                 ylab("Number of seqs") +
                 theme_light() +
                 theme(strip.background = element_blank(),
                       strip.text = element_text(color = "black", face = "bold"),
                       axis.text.x = element_text(angle = 45, hjust = 1),
                       axis.title.x = element_blank())

seqlength_barplots <- ggplot(track_mergers_length_long, aes(y = length)) +
                      facet_wrap(~ samples, ncol = 5, scales = "free") +
                      geom_boxplot(fill = "gray80", alpha = 0.7) +
                      ylab("Read length") +
                      theme_light() +
                      theme(strip.background = element_blank(),
                            strip.text = element_text(color = "black", face = "bold"),
                            axis.title.x = element_blank())

###############################################################################
### 18. Save plots
###############################################################################

filename_nseq <- file.path(OUTPUT_DIR,"plots/nseq_barplot.pdf")
filename_seqlength <- file.path(OUTPUT_DIR,"plots/seq_length_hist.pdf")

ggsave(nseq_barplots, filename = filename_nseq,
       device = "pdf", width = 18, height = 30)

ggsave(seqlength_barplots, filename = filename_seqlength,
       device = "pdf", width = 10, height = 30)

###############################################################################
### 19. Save track stats tables
###############################################################################

# nseq
filename_nseq <- file.path(OUTPUT_DIR,"tables/nseq_counts.csv")
write.csv(file = filename_nseq, track_n_seqs)

# seq length
track_mergers_length_stats <- track_mergers_length_long %>% 
                              group_by(samples) %>%
                              summarize(mean = mean(length), 
                                        sd = sd(length),
                                        max = max(length), 
                                        min = min(length))

filename_seqlength_stats <- file.path(OUTPUT_DIR,"tables/seq_length_stats.csv")
write.csv(file = filename_seqlength_stats, track_mergers_length_stats)

###############################################################################
### 20. Save R workspace
###############################################################################

if (SAVE_WORKSPACE == T) {
  
  filename_RData <-paste(OUTPUT_DIR,"/.RData", sep = "") 
  save.image(file = filename_RData)

}
