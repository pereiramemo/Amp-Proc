#!/usr/bin/env bash

###############################################################################
# quality_check_fastp.sh
###############################################################################

set -euo pipefail

###############################################################################
# 1. Environment
###############################################################################

# Colors
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'

log()       { echo -e "[INFO] $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*" >&2; }
log_error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }

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

    --reads1 FILE
        R1 FASTQ file (required)

    --reads2 FILE
        R2 FASTQ file (required)

    --output_dir CHAR
        Directory to output generated data (required)

    --nslots NUM
        Number of threads to use [default=12]

    --min_length NUM
        Minimum read length (only for reporting) [default=50]

    --qualified_quality_phred NUM
        Minimum quality value for qualified base (Phred score, only for reporting) [default=20]

    --unqualified_percent_limit NUM
        Maximum percent of unqualified bases allowed (only for reporting) [default=40]

    --disable_adapter_trimming t|f
        Disable adapter trimming in report [default=t]

    --html_report t|f
        Generate HTML report [default=t]

    --json_report t|f
        Generate JSON report [default=t]

    --overwrite t|f
        Overwrite previous directory [default=f]

Examples:
    $(basename "$0") \\
        --reads1 raw_data/sample1_R1.fastq.gz \\
        --reads2 raw_data/sample1_R2.fastq.gz \\
        --output_dir results/qc/sample1

EOF
}

###############################################################################
# 3. Define default parameters
###############################################################################

READS1=""
READS2=""
OUTPUT_DIR=""
NSLOTS=12
MIN_LENGTH=50
QUALIFIED_QUALITY_PHRED=20
UNQUALIFIED_PERCENT_LIMIT=40
DISABLE_ADAPTER_TRIMMING="t"
HTML_REPORT="t"
JSON_REPORT="t"
OVERWRITE="f"

# DEV ONLY — comment out before production use

READS1="/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R1_001_redu.fastq.gz"
READS2="/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R2_001_redu.fastq.gz"
OUTPUT_DIR="/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/01_fastp/sample1"
NSLOTS=4
OVERWRITE="t"

###############################################################################
# 4. Parse arguments
###############################################################################

if [[ $# -eq 0 ]]; then
    show_usage
    exit 1
fi

ARGS=$(
  getopt -o h \
  -l help,reads1:,reads2:,output_dir:,nslots:,min_length:,\
qualified_quality_phred:,unqualified_percent_limit:,disable_adapter_trimming:,\
html_report:,json_report:,overwrite: \
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
    --reads1) READS1="$2"; shift 2 ;;
    --reads2) READS2="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    --nslots) NSLOTS="$2"; shift 2 ;;
    --min_length) MIN_LENGTH="$2"; shift 2 ;;
    --qualified_quality_phred) QUALIFIED_QUALITY_PHRED="$2"; shift 2 ;;
    --unqualified_percent_limit) UNQUALIFIED_PERCENT_LIMIT="$2"; shift 2 ;;
    --disable_adapter_trimming) DISABLE_ADAPTER_TRIMMING="$2"; shift 2 ;;
    --html_report) HTML_REPORT="$2"; shift 2 ;;
    --json_report) JSON_REPORT="$2"; shift 2 ;;
    --overwrite) OVERWRITE="$2"; shift 2 ;;
    --) shift; break ;;
    *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
  esac
done

###############################################################################
# 5. Validate parameters
###############################################################################

if [[ -z "${READS1}" ]]; then
    log_error "--reads1 is required."
    exit 1
fi

if [[ -z "${READS2}" ]]; then
    log_error "--reads2 is required."
    exit 1
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
    log_error "--output_dir is required."
    exit 1
fi

if [[ ! -f "${READS1}" ]]; then
    log_error "R1 file does not exist: ${READS1}"
    exit 1
fi

if [[ ! -f "${READS2}" ]]; then
    log_error "R2 file does not exist: ${READS2}"
    exit 1
fi

for flag in DISABLE_ADAPTER_TRIMMING HTML_REPORT JSON_REPORT OVERWRITE; do
    if ! [[ "${!flag}" =~ ^[tf]$ ]]; then
        log_error "Flag --$(echo "${flag}" | tr 'A-Z_' 'a-z-') must be 't' or 'f' (got '${!flag}')."
        exit 1
    fi
done

for num in NSLOTS MIN_LENGTH QUALIFIED_QUALITY_PHRED UNQUALIFIED_PERCENT_LIMIT; do
    if ! [[ "${!num}" =~ ^[0-9]+$ ]]; then
        log_error "Parameter ${num} must be numeric (got '${!num}')."
        exit 1
    fi
done

###############################################################################
# 6. Check dependencies
###############################################################################

log "Checking dependencies..."

if ! check_cmd fastp; then
    log_error "fastp not found. Please install it or add it to PATH."
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
        log_error "Output directory exists: ${OUTPUT_DIR}. Use --overwrite t to overwrite."
        exit 1
    fi
fi

mkdir -p "${OUTPUT_DIR}/reports"
mkdir -p "${OUTPUT_DIR}/stats"

###############################################################################
# 8. Derive sample name
###############################################################################

SAMPLE_NAME=$(basename "${READS1}")
SAMPLE_NAME="${SAMPLE_NAME%.fastq.gz}"
SAMPLE_NAME="${SAMPLE_NAME%.fq.gz}"
SAMPLE_NAME="${SAMPLE_NAME%.fastq}"
SAMPLE_NAME="${SAMPLE_NAME%.fq}"

log "Processing sample: ${SAMPLE_NAME}"

###############################################################################
# 9. Run fastp
###############################################################################

HTML_OUT="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_fastp.html"
JSON_OUT="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_fastp.json"
LOG_OUT="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_fastp.log"

FASTP_CMD=(
    fastp
    -i "${READS1}"
    -I "${READS2}"
    -w "${NSLOTS}"
    --disable_quality_filtering
    --disable_length_filtering
    --disable_trim_poly_g
    --json "${JSON_OUT}"
)

if [[ "${HTML_REPORT}" == "t" ]]; then
    FASTP_CMD+=(--html "${HTML_OUT}")
else
    FASTP_CMD+=(--html /dev/null)
fi

if [[ "${DISABLE_ADAPTER_TRIMMING}" == "t" ]]; then
    FASTP_CMD+=(--disable_adapter_trimming)
fi

if ! "${FASTP_CMD[@]}" 2>&1 | tee "${LOG_OUT}"; then
    log_error "fastp failed for sample ${SAMPLE_NAME}"
    exit 1
fi

###############################################################################
# 10. Extract summary statistics
###############################################################################

read -r TOTAL_READS TOTAL_BASES Q20_BASES Q20_RATE Q30_BASES Q30_RATE R1_MEAN_LEN R2_MEAN_LEN GC_CONTENT < <(
    python3 - "${JSON_OUT}" << 'EOF'
import json, sys
with open(sys.argv[1]) as f:
    bf = json.load(f)["summary"]["before_filtering"]
print(bf["total_reads"], bf["total_bases"],
      bf["q20_bases"], bf["q20_rate"],
      bf["q30_bases"], bf["q30_rate"],
      bf["read1_mean_length"], bf["read2_mean_length"],
      bf["gc_content"])
EOF
)

SUMMARY_FILE="${OUTPUT_DIR}/stats/summary.tsv"
echo -e "sample\ttotal_reads\ttotal_bases\tq20_bases\tq20_rate\tq30_bases\tq30_rate\tread1_mean_length\tread2_mean_length\tgc_content" > "${SUMMARY_FILE}"
echo -e "${SAMPLE_NAME}\t${TOTAL_READS}\t${TOTAL_BASES}\t${Q20_BASES}\t${Q20_RATE}\t${Q30_BASES}\t${Q30_RATE}\t${R1_MEAN_LEN}\t${R2_MEAN_LEN}\t${GC_CONTENT}" >> "${SUMMARY_FILE}"

if [[ "${JSON_REPORT}" != "t" ]]; then
    rm -f "${JSON_OUT}"
fi

###############################################################################
# 11. Generate summary report
###############################################################################

log "Generating summary report..."

SUMMARY_TXT="${OUTPUT_DIR}/summary_report.txt"

cat > "${SUMMARY_TXT}" << EOF
================================================================================
Quality Check Report - fastp
================================================================================
Date: $(date)
Sample: ${SAMPLE_NAME}
R1: ${READS1}
R2: ${READS2}
Output directory: ${OUTPUT_DIR}

Parameters:
-----------
Threads: ${NSLOTS}
Mode: Report only (no filtering applied)
Adapter trimming: $([ "${DISABLE_ADAPTER_TRIMMING}" == "f" ] && echo "enabled" || echo "disabled")

Output files:
-------------
- HTML report: ${HTML_OUT}
- JSON report: ${JSON_OUT}
- Processing log: ${LOG_OUT}
- Summary statistics: ${SUMMARY_FILE}
- This report: ${SUMMARY_TXT}

================================================================================
EOF

echo -e "\nSummary Statistics:\n" >> "${SUMMARY_TXT}"
column -t -s $'\t' "${SUMMARY_FILE}" >> "${SUMMARY_TXT}"

cat "${SUMMARY_TXT}"

###############################################################################
# 12. End
###############################################################################

log "${GREEN}quality_check_fastp.sh completed successfully${NC}"
log "Reports available in: ${OUTPUT_DIR}/reports/"
log "Summary statistics: ${SUMMARY_FILE}"
