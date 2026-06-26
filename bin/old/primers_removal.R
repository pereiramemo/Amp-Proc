#!/usr/bin/env Rscript

# This code is a modified version of the DADA2 tutorial: 
# https://benjjneb.github.io/dada2/ITS_workflow.html

###############################################################################
### 1. Def. env
###############################################################################

suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ShortRead)))
suppressMessages(suppressWarnings(library(Biostrings)))

REPO <- "/home/epereira/workspace/repositories/amplicon_pipelines/"
toolbox <- file.path(REPO, "scripts_v2/toolbox.R")
source(toolbox)

###############################################################################
### 2. Get parameters
###############################################################################

args = commandArgs(trailingOnly=TRUE)

INPUT_DIR <- args[1]
OUTPUT_DIR <- args[2]
PATTERN_R1 <- args[3]
PATTERN_R2 <- args[4]
PRIMER_FWD <- args[5]
PRIMER_REV <- args[6]
NSLOTS <- args[7] %>% as.numeric
SAVE_WORKSPACE <- args[8] %>% as.logical()

# INPUT_DIR <- "/home/epereira/workspace/indicadores_cuencas_2018/data/toy_dataset/"
# OUTPUT_DIR <- "/home/epereira/workspace/indicadores_cuencas_2018/data/toy_dataset/"
# PRIMER_FWD <- "GTGYCAGCMGCCGCGGTAA"
# PRIMER_REV <- "CCGYCAATTYMTTTRAGTTT"
# PATTERN_R1 <- "_L001_R1_001_redu.fastq.gz"
# PATTERN_R2 <- "_L001_R2_001_redu.fastq.gz"
# NSLOTS <- 12

###############################################################################
### 3. Load data
###############################################################################

print("Loading data ...")
rawR1 <- sort(list.files(INPUT_DIR, pattern = PATTERN_R1, full.names = T))
rawR2 <- sort(list.files(INPUT_DIR, pattern = PATTERN_R2, full.names = T))

sample.names <- basename(rawR1) %>%
                sub(pattern = PATTERN_R1, replacement = "")

###############################################################################
### 4. Create output dirs
###############################################################################

dir.create(OUTPUT_DIR)

###############################################################################
### 5. Create output file names
###############################################################################

print("Creating output files ...")
rmprimerR1 <- file.path(OUTPUT_DIR, paste(sample.names, "R1_rmprimer.fastq.gz", sep = "_"))
rmprimerR2 <- file.path(OUTPUT_DIR, paste(sample.names, "R2_rmprimer.fastq.gz", sep = "_"))

##############################################################################
### 6. Check presence of primers
###############################################################################

print("Searching primers brefore removal ...")
FWD_orients <- all_orients(PRIMER_FWD)
REV_orients <- all_orients(PRIMER_REV)

counts_log <- rbind(FWD.ForwardReads = sapply(FWD_orients, primer_hits, INPUT = rawR1[[1]]), 
                    FWD.ReverseReads = sapply(FWD_orients, primer_hits, INPUT = rawR2[[1]]), 
                    REV.ForwardReads = sapply(REV_orients, primer_hits, INPUT = rawR1[[1]]), 
                    REV.ReverseReads = sapply(REV_orients, primer_hits, INPUT = rawR2[[1]]))

filename_counts_log <- file.path(OUTPUT_DIR,"primer_counts_before_rmprimer.log")
write.csv(file = filename_counts_log, counts_log)

###############################################################################
### 7. Remove primers
###############################################################################

print("Removing primres ...")
cutadapt_runner(primer_fwd = PRIMER_FWD,
                primer_rev = PRIMER_REV,
                input_r1 = rawR1,
                input_r2 = rawR2,
                output_r1 = rmprimerR1,
                output_r2 = rmprimerR2,
                nslots = NSLOTS)

###############################################################################
### 8. Check removal
###############################################################################

print("Searching primers after removal ...")

FWD_orients <- all_orients(PRIMER_FWD)
REV_orients <- all_orients(PRIMER_REV)

counts_log <- rbind(FWD.ForwardReads = sapply(FWD_orients, primer_hits, INPUT = rmprimerR1[[1]]), 
                    FWD.ReverseReads = sapply(FWD_orients, primer_hits, INPUT = rmprimerR2[[1]]), 
                    REV.ForwardReads = sapply(REV_orients, primer_hits, INPUT = rmprimerR1[[1]]), 
                    REV.ReverseReads = sapply(REV_orients, primer_hits, INPUT = rmprimerR2[[1]]))

filename_counts_log <- file.path(OUTPUT_DIR,"primer_counts_after_rmprimer.log")
write.csv(file = filename_counts_log, counts_log)


###############################################################################
###  9. Save R workspace
###############################################################################

if (SAVE_WORKSPACE == T) {
  
  filename_RData <-paste(OUTPUT_DIR,"/.RData", sep = "") 
  save.image(file = filename_RData)
  
}
