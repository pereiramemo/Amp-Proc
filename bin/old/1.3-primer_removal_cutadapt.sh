#!/usr/bin/env bash

###############################################################################
# primer_removal_cutadapt.sh
###############################################################################

set -euo pipefail

###############################################################################
# 1. Environment
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/conf.sh"

# Colors
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'

log()       { echo -e "[INFO] $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*" >&2; }
log_error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }

# Helpers
check_cmd() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log_error "Required tool '$1' not found in PATH."
        return 1
    fi
}

###############################################################################
# 2. Define help
###############################################################################

show_usage() {
    cat << EOF
Usage: $(basename "$0") <options>

Options:
    --help
        Print this help message and exit

    --input_dir=CHAR
        Directory containing input FASTQ files (required)

    --output_dir=CHAR
        Directory to output trimmed FASTQ files (required)

    --pattern_r1=CHAR
        Pattern for R1 FASTQ files [default=_R1_001.fastq.gz]

    --pattern_r2=CHAR
        Pattern for R2 FASTQ files [default=_R2_001.fastq.gz]

    --primer_fwd=CHAR
        Forward primer sequence (5' to 3') (required)

    --primer_rev=CHAR
        Reverse primer sequence (5' to 3') (required)

    --nslots=NUM
        Number of threads to use [default=12]

    --error_rate=NUM
        Maximum allowed error rate (mismatches/length) [default=0.1]

    --min_overlap=NUM
        Minimum overlap between primer and read [default=3]

    --min_length=NUM
        Discard reads shorter than this after trimming [default=50]

    --discard_untrimmed=t|f
        Discard reads where no primer was found [default=f]

    --compress=t|f
        Compress output files with gzip [default=t]

    --overwrite=t|f
        Overwrite previous directory [default=f]

Examples:
    # Basic usage
    $(basename "$0") \\
        --input_dir=raw_data/ \\
        --output_dir=results/trimmed \\
        --primer_fwd=GTGYCAGCMGCCGCGGTAA \\
        --primer_rev=CCGYCAATTYMTTTRAGTTT

    # Custom settings with untrimmed read filtering
    $(basename "$0") \\
        --input_dir=raw_data/ \\
        --output_dir=results/trimmed \\
        --primer_fwd=GTGYCAGCMGCCGCGGTAA \\
        --primer_rev=CCGYCAATTYMTTTRAGTTT \\
        --nslots=16 \\
        --min_length=75 \\
        --discard_untrimmed=t

EOF
}

###############################################################################
# 3. Define default parameters
###############################################################################

INPUT_DIR=""
OUTPUT_DIR=""
PATTERN_R1="_R1_001.fastq.gz"
PATTERN_R2="_R2_001.fastq.gz"
PRIMER_FWD=""
PRIMER_REV=""
NSLOTS=12
ERROR_RATE=0.1
MIN_OVERLAP=3
MIN_LENGTH=50
DISCARD_UNTRIMMED="f"
COMPRESS="t"
OVERWRITE="f"

###############################################################################
# 4. Parse arguments
###############################################################################

if [[ $# -eq 0 ]]; then
    show_usage
    exit 1
fi

ARGS=$(
  getopt -o h \
  -l help,input_dir:,output_dir:,pattern_r1:,pattern_r2:,primer_fwd:,primer_rev:,\
nslots:,error_rate:,min_overlap:,min_length:,discard_untrimmed:,compress:,overwrite: \
  -n "$(basename "$0")" -- "$@" \
  ) || {
  log_error "Failed to parse arguments."
  show_usage
  exit 1
}

eval set -- "${ARGS}"

while true; do
  case "$1" in
    -h|--help) show_usage; exit 0 ;;
    --input_dir) INPUT_DIR="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    --pattern_r1) PATTERN_R1="$2"; shift 2 ;;
    --pattern_r2) PATTERN_R2="$2"; shift 2 ;;
    --primer_fwd) PRIMER_FWD="$2"; shift 2 ;;
    --primer_rev) PRIMER_REV="$2"; shift 2 ;;
    --nslots) NSLOTS="$2"; shift 2 ;;
    --error_rate) ERROR_RATE="$2"; shift 2 ;;
    --min_overlap) MIN_OVERLAP="$2"; shift 2 ;;
    --min_length) MIN_LENGTH="$2"; shift 2 ;;
    --discard_untrimmed) DISCARD_UNTRIMMED="$2"; shift 2 ;;
    --compress) COMPRESS="$2"; shift 2 ;;
    --overwrite) OVERWRITE="$2"; shift 2 ;;
    --) shift; break ;;
    *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
  esac
done

###############################################################################
# 5. Validate parameters
###############################################################################

# Required parameters
if [[ -z "${INPUT_DIR}" ]]; then
    log_error "--input_dir is required."
    exit 1
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
    log_error "--output_dir is required."
    exit 1
fi

if [[ -z "${PRIMER_FWD}" ]]; then
    log_error "--primer_fwd is required."
    exit 1
fi

if [[ -z "${PRIMER_REV}" ]]; then
    log_error "--primer_rev is required."
    exit 1
fi

# Check input directory exists
if [[ ! -d "${INPUT_DIR}" ]]; then
    log_error "Input directory does not exist: ${INPUT_DIR}"
    exit 1
fi

# Validate boolean flags
for flag in DISCARD_UNTRIMMED COMPRESS OVERWRITE; do
    if ! [[ "${!flag}" =~ ^[tf]$ ]]; then
        log_error "Flag --$(echo "${flag}" | tr 'A-Z_' 'a-z-') must be 't' or 'f' (got '${!flag}')."
        exit 1
    fi
done

# Validate numeric parameters
for num in NSLOTS MIN_OVERLAP MIN_LENGTH; do
    if ! [[ "${!num}" =~ ^[0-9]+$ ]]; then
        log_error "Parameter ${num} must be numeric (got '${!num}')."
        exit 1
    fi
done

# Validate error rate (must be a number between 0 and 1)
if ! [[ "${ERROR_RATE}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
    log_error "Parameter ERROR_RATE must be numeric (got '${ERROR_RATE}')."
    exit 1
fi

###############################################################################
# 6. Check dependencies
###############################################################################

log "Checking dependencies..."

if ! check_cmd cutadapt; then
    log_error "cutadapt not found. Please install it or add it to PATH."
    exit 1
fi

###############################################################################
# 7. Prepare output directory
###############################################################################

if [[ -d "${OUTPUT_DIR}" ]]; then
    if [[ "${OVERWRITE}" == "t" ]]; then
        log_warn "Overwriting existing directory: ${OUTPUT_DIR}"
        rm -rf "${OUTPUT_DIR}"
    else
        log_error "Output directory exists: ${OUTPUT_DIR}. Use --overwrite=t to overwrite."
        exit 1
    fi
fi

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/trimmed"
mkdir -p "${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/stats"

###############################################################################
# 8. Find input files
###############################################################################

log "Searching for input files..."

# Use arrays to hold file paths
readarray -t R1_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*${PATTERN_R1}" | sort)
readarray -t R2_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*${PATTERN_R2}" | sort)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    log_error "No R1 files found with pattern: *${PATTERN_R1}"
    exit 1
fi

if [[ ${#R2_FILES[@]} -eq 0 ]]; then
    log_error "No R2 files found with pattern: *${PATTERN_R2}"
    exit 1
fi

if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
    log_error "Number of R1 and R2 files do not match (R1: ${#R1_FILES[@]}, R2: ${#R2_FILES[@]})"
    exit 1
fi

log "Found ${#R1_FILES[@]} sample pairs"

###############################################################################
# 9. Prepare primer sequences
###############################################################################

log "Forward primer: ${PRIMER_FWD}"
log "Reverse primer: ${PRIMER_REV}"

###############################################################################
# 10. Process samples
###############################################################################

SUMMARY_FILE="${OUTPUT_DIR}/stats/summary.tsv"
echo -e "sample\ttotal_pairs\ttrimmed_pairs\tpercent_trimmed" > "${SUMMARY_FILE}"

for i in "${!R1_FILES[@]}"; do
    R1="${R1_FILES[$i]}"
    R2="${R2_FILES[$i]}"
    
    # Extract sample name by removing pattern
    SAMPLE_NAME=$(basename "${R1}" "${PATTERN_R1}")
    
    log "Processing sample ${SAMPLE_NAME} ($(( i + 1 ))/${#R1_FILES[@]})..."
    
    # Define output files
    if [[ "${COMPRESS}" == "t" ]]; then
        R1_OUT="${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_R1_trimmed.fastq.gz"
        R2_OUT="${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_R2_trimmed.fastq.gz"
    else
        R1_OUT="${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_R1_trimmed.fastq"
        R2_OUT="${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_R2_trimmed.fastq"
    fi
    
    LOG_FILE="${OUTPUT_DIR}/logs/${SAMPLE_NAME}_cutadapt.log"
    
    # Build cutadapt command
    # -g for 5' adapter on R1 (forward primer)
    # -G for 5' adapter on R2 (reverse primer)
    CUTADAPT_CMD=(
        cutadapt
        -g "${PRIMER_FWD}"
        -G "${PRIMER_REV}"
        -o "${R1_OUT}"
        -p "${R2_OUT}"
        -j "${NSLOTS}"
        -e "${ERROR_RATE}"
        -O "${MIN_OVERLAP}"
        -m "${MIN_LENGTH}"
    )
    
    # Add discard untrimmed option if requested
    if [[ "${DISCARD_UNTRIMMED}" == "t" ]]; then
        CUTADAPT_CMD+=(--discard-untrimmed)
    fi
    
    # Add input files
    CUTADAPT_CMD+=(
        "${R1}"
        "${R2}"
    )
    
    # Run cutadapt
    if ! "${CUTADAPT_CMD[@]}" > "${LOG_FILE}" 2>&1; then
        log_error "cutadapt failed for sample ${SAMPLE_NAME}"
        log_error "Check log file: ${LOG_FILE}"
        exit 1
    fi
    
    # Extract summary statistics from log file
    TOTAL_PAIRS=$(grep "Total read pairs processed:" "${LOG_FILE}" | awk '{print $NF}' | tr -d ',')
    PAIRS_WRITTEN=$(grep "Pairs written (passing filters):" "${LOG_FILE}" | awk '{print $(NF-1)}' | tr -d ',')
    
    # Calculate percent trimmed
    if [[ -n "${TOTAL_PAIRS}" && -n "${PAIRS_WRITTEN}" && ${TOTAL_PAIRS} -gt 0 ]]; then
        PERCENT_TRIMMED=$(awk -v pairs="${PAIRS_WRITTEN}" -v total="${TOTAL_PAIRS}" 'BEGIN {printf "%.2f", (pairs/total)*100}')
    else
        PERCENT_TRIMMED="0.00"
    fi
    
    # Write to summary file
    echo -e "${SAMPLE_NAME}\t${TOTAL_PAIRS:-0}\t${PAIRS_WRITTEN:-0}\t${PERCENT_TRIMMED}" >> "${SUMMARY_FILE}"
    
    log "Completed processing ${SAMPLE_NAME}"
done

###############################################################################
# 11. Generate summary report
###############################################################################

log "Generating summary report..."

SUMMARY_TXT="${OUTPUT_DIR}/summary_report.txt"

cat > "${SUMMARY_TXT}" << EOF
================================================================================
Primer Removal Report - cutadapt
================================================================================
Date: $(date)
Input directory: ${INPUT_DIR}
Output directory: ${OUTPUT_DIR}
Number of samples: ${#R1_FILES[@]}

Primers:
--------
Forward primer (5'-3'): ${PRIMER_FWD}
Reverse primer (5'-3'): ${PRIMER_REV}

Parameters:
-----------
Pattern R1: ${PATTERN_R1}
Pattern R2: ${PATTERN_R2}
Threads: ${NSLOTS}
Error rate: ${ERROR_RATE}
Minimum overlap: ${MIN_OVERLAP}
Minimum length: ${MIN_LENGTH}
Discard untrimmed: $([ "${DISCARD_UNTRIMMED}" == "t" ] && echo "yes" || echo "no")
Compress output: $([ "${COMPRESS}" == "t" ] && echo "yes" || echo "no")

Output files:
-------------
- Trimmed reads: ${OUTPUT_DIR}/trimmed/*_trimmed.fastq$([ "${COMPRESS}" == "t" ] && echo ".gz" || echo "")
- Processing logs: ${OUTPUT_DIR}/logs/*_cutadapt.log
- Summary statistics: ${OUTPUT_DIR}/stats/summary.tsv
- This report: ${OUTPUT_DIR}/summary_report.txt

================================================================================
EOF

if [[ -f "${SUMMARY_FILE}" ]]; then
    echo -e "\nSummary Statistics:\n" >> "${SUMMARY_TXT}"
    column -t -s $'\t' "${SUMMARY_FILE}" >> "${SUMMARY_TXT}"
fi

cat "${SUMMARY_TXT}"

###############################################################################
# 12. End
###############################################################################

log "${GREEN}primer_removal_cutadapt.sh completed successfully${NC}"
log "Trimmed reads available in: ${OUTPUT_DIR}/trimmed/"
log "Summary statistics: ${SUMMARY_FILE}"
