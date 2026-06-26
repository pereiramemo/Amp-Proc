#!/usr/bin/env bash

###############################################################################
# 2.2-vsearch_pipeline.sh
###############################################################################

set -euo pipefail

###############################################################################
# 1. Environment
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/conf.sh"

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

OTU pipeline: fastp quality filtering → vsearch PE merging → chimera
removal → OTU clustering → OTU table.

Options:
    --help
        Print this help message and exit

    --input_dir=CHAR
        Directory containing paired-end FASTQ files (required)

    --output_dir=CHAR
        Directory to write all output (required)

    --pattern_r1=CHAR
        Filename pattern for R1 reads [default=_R1_001.fastq.gz]

    --pattern_r2=CHAR
        Filename pattern for R2 reads [default=_R2_001.fastq.gz]

    --nslots=NUM
        Number of threads [default=12]

    --min_quality=NUM
        Minimum base Phred quality for fastp trimming [default=20]

    --min_length=NUM
        Minimum read length after fastp trimming and after merging [default=50]

    --fastq_minovlen=NUM
        Minimum overlap length for vsearch PE merging [default=12]

    --fastq_maxdiffs=NUM
        Maximum mismatches allowed in the overlap region [default=2]

    --fastq_maxee=NUM
        Maximum expected errors per merged read [default=1.0]

    --min_size=NUM
        Minimum abundance to retain a sequence after dereplication [default=2]

    --chimera_method=CHAR
        Chimera detection method: denovo or ref [default=denovo]

    --ref_db=CHAR
        Reference FASTA for chimera checking (required when chimera_method=ref)

    --identity=NUM
        OTU clustering identity threshold (0–1) [default=0.97]

    --overwrite=t|f
        Overwrite existing output directory [default=f]

Examples:
    # De novo chimera removal (no reference required)
    $(basename "$0") \\
        --input_dir=results/trimmed/ \\
        --output_dir=results/otus/

    # Reference-based chimera removal, custom identity
    $(basename "$0") \\
        --input_dir=results/trimmed/ \\
        --output_dir=results/otus/ \\
        --chimera_method=ref \\
        --ref_db=/path/to/silva_138_99_16S.fasta \\
        --identity=0.97 \\
        --nslots=16

EOF
}

###############################################################################
# 3. Default parameters
###############################################################################

INPUT_DIR=""
OUTPUT_DIR=""
PATTERN_R1="_R1_001.fastq.gz"
PATTERN_R2="_R2_001.fastq.gz"
NSLOTS=12
MIN_QUALITY=20
MIN_LENGTH=50
FASTQ_MINOVLEN=12
FASTQ_MAXDIFFS=2
FASTQ_MAXEE=1.0
MIN_SIZE=2
CHIMERA_METHOD="denovo"
REF_DB=""
IDENTITY=0.97
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
  -l help,input_dir:,output_dir:,pattern_r1:,pattern_r2:,nslots:,\
min_quality:,min_length:,fastq_minovlen:,fastq_maxdiffs:,fastq_maxee:,\
min_size:,chimera_method:,ref_db:,identity:,overwrite: \
  -n "$(basename "$0")" -- "$@"
) || {
  log_error "Failed to parse arguments."
  show_usage
  exit 1
}

eval set -- "${ARGS}"

while true; do
  case "$1" in
    -h|--help)           show_usage; exit 0 ;;
    --input_dir)         INPUT_DIR="$2";       shift 2 ;;
    --output_dir)        OUTPUT_DIR="$2";      shift 2 ;;
    --pattern_r1)        PATTERN_R1="$2";      shift 2 ;;
    --pattern_r2)        PATTERN_R2="$2";      shift 2 ;;
    --nslots)            NSLOTS="$2";          shift 2 ;;
    --min_quality)       MIN_QUALITY="$2";     shift 2 ;;
    --min_length)        MIN_LENGTH="$2";      shift 2 ;;
    --fastq_minovlen)    FASTQ_MINOVLEN="$2";  shift 2 ;;
    --fastq_maxdiffs)    FASTQ_MAXDIFFS="$2";  shift 2 ;;
    --fastq_maxee)       FASTQ_MAXEE="$2";     shift 2 ;;
    --min_size)          MIN_SIZE="$2";        shift 2 ;;
    --chimera_method)    CHIMERA_METHOD="$2";  shift 2 ;;
    --ref_db)            REF_DB="$2";          shift 2 ;;
    --identity)          IDENTITY="$2";        shift 2 ;;
    --overwrite)         OVERWRITE="$2";       shift 2 ;;
    --) shift; break ;;
    *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
  esac
done

###############################################################################
# 5. Validate parameters
###############################################################################

if [[ -z "${INPUT_DIR}" ]]; then
    log_error "--input_dir is required."
    exit 1
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
    log_error "--output_dir is required."
    exit 1
fi

if [[ ! -d "${INPUT_DIR}" ]]; then
    log_error "Input directory does not exist: ${INPUT_DIR}"
    exit 1
fi

if ! [[ "${OVERWRITE}" =~ ^[tf]$ ]]; then
    log_error "--overwrite must be 't' or 'f' (got '${OVERWRITE}')."
    exit 1
fi

for num in NSLOTS MIN_QUALITY MIN_LENGTH FASTQ_MINOVLEN FASTQ_MAXDIFFS MIN_SIZE; do
    if ! [[ "${!num}" =~ ^[0-9]+$ ]]; then
        log_error "--$(echo "${num}" | tr 'A-Z_' 'a-z-') must be a positive integer (got '${!num}')."
        exit 1
    fi
done

for flt in FASTQ_MAXEE IDENTITY; do
    if ! [[ "${!flt}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
        log_error "--$(echo "${flt}" | tr 'A-Z_' 'a-z-') must be numeric (got '${!flt}')."
        exit 1
    fi
done

if [[ "${CHIMERA_METHOD}" != "denovo" && "${CHIMERA_METHOD}" != "ref" ]]; then
    log_error "--chimera_method must be 'denovo' or 'ref' (got '${CHIMERA_METHOD}')."
    exit 1
fi

if [[ "${CHIMERA_METHOD}" == "ref" ]]; then
    if [[ -z "${REF_DB}" ]]; then
        log_error "--ref_db is required when --chimera_method=ref."
        exit 1
    fi
    if [[ ! -f "${REF_DB}" ]]; then
        log_error "Reference database not found: ${REF_DB}"
        exit 1
    fi
fi

###############################################################################
# 6. Check dependencies
###############################################################################

log "Checking dependencies..."

for tool in fastp vsearch; do
    if ! check_cmd "${tool}"; then
        log_error "${tool} not found. Please install it or add it to PATH."
        exit 1
    fi
done

###############################################################################
# 7. Prepare output directories
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

mkdir -p \
    "${OUTPUT_DIR}/01_fastp" \
    "${OUTPUT_DIR}/02_merged" \
    "${OUTPUT_DIR}/03_filtered" \
    "${OUTPUT_DIR}/04_derep" \
    "${OUTPUT_DIR}/05_chimera" \
    "${OUTPUT_DIR}/06_otus" \
    "${OUTPUT_DIR}/logs" \
    "${OUTPUT_DIR}/stats"

###############################################################################
# 8. Find input files
###############################################################################

log "Searching for input files..."

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
    log_error "R1/R2 file count mismatch (R1: ${#R1_FILES[@]}, R2: ${#R2_FILES[@]})"
    exit 1
fi

log "Found ${#R1_FILES[@]} sample pair(s)"

###############################################################################
# 9. Quality filtering with fastp (per sample)
###############################################################################

log "Step 1/6: Quality filtering with fastp..."

FASTP_SUMMARY="${OUTPUT_DIR}/stats/fastp_summary.tsv"
echo -e "sample\treads_in\treads_out\tpercent_passed" > "${FASTP_SUMMARY}"

for i in "${!R1_FILES[@]}"; do
    R1="${R1_FILES[$i]}"
    R2="${R2_FILES[$i]}"
    SAMPLE_NAME=$(basename "${R1}" "${PATTERN_R1}")

    log "  [$(( i + 1 ))/${#R1_FILES[@]}] fastp: ${SAMPLE_NAME}"

    R1_OUT="${OUTPUT_DIR}/01_fastp/${SAMPLE_NAME}_R1.fastq.gz"
    R2_OUT="${OUTPUT_DIR}/01_fastp/${SAMPLE_NAME}_R2.fastq.gz"
    JSON_OUT="${OUTPUT_DIR}/01_fastp/${SAMPLE_NAME}_fastp.json"

    fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${R1_OUT}" \
        -O "${R2_OUT}" \
        -j "${JSON_OUT}" \
        --html /dev/null \
        -w "${NSLOTS}" \
        -q "${MIN_QUALITY}" \
        -l "${MIN_LENGTH}" \
        --detect_adapter_for_pe \
        2> "${OUTPUT_DIR}/logs/${SAMPLE_NAME}_fastp.log"

    if [[ -f "${JSON_OUT}" ]]; then
        READS_IN=$(grep -Po '"total_reads":\s*\K[0-9]+' "${JSON_OUT}" | head -1 || echo "0")
        READS_OUT=$(grep -Po '"total_reads":\s*\K[0-9]+' "${JSON_OUT}" | tail -1 || echo "0")
        if [[ "${READS_IN:-0}" -gt 0 ]]; then
            PCT=$(awk -v a="${READS_OUT:-0}" -v b="${READS_IN}" \
                  'BEGIN {printf "%.2f", (a/b)*100}')
        else
            PCT="0.00"
        fi
        echo -e "${SAMPLE_NAME}\t${READS_IN:-0}\t${READS_OUT:-0}\t${PCT}" \
            >> "${FASTP_SUMMARY}"
    fi
done

###############################################################################
# 10. Merge paired-end reads with vsearch (per sample)
###############################################################################

log "Step 2/6: Merging paired-end reads with vsearch..."

MERGE_SUMMARY="${OUTPUT_DIR}/stats/merge_summary.tsv"
echo -e "sample\tpairs_in\tpairs_merged\tpercent_merged" > "${MERGE_SUMMARY}"

ALL_FILTERED_FASTA="${OUTPUT_DIR}/03_filtered/all_samples.fasta"
> "${ALL_FILTERED_FASTA}"

for i in "${!R1_FILES[@]}"; do
    R1="${R1_FILES[$i]}"
    SAMPLE_NAME=$(basename "${R1}" "${PATTERN_R1}")

    log "  [$(( i + 1 ))/${#R1_FILES[@]}] merge: ${SAMPLE_NAME}"

    R1_TRIMMED="${OUTPUT_DIR}/01_fastp/${SAMPLE_NAME}_R1.fastq.gz"
    R2_TRIMMED="${OUTPUT_DIR}/01_fastp/${SAMPLE_NAME}_R2.fastq.gz"
    MERGED_FASTQ="${OUTPUT_DIR}/02_merged/${SAMPLE_NAME}_merged.fastq"
    MERGE_LOG="${OUTPUT_DIR}/logs/${SAMPLE_NAME}_merge.log"

    vsearch \
        --fastq_mergepairs "${R1_TRIMMED}" \
        --reverse "${R2_TRIMMED}" \
        --fastqout "${MERGED_FASTQ}" \
        --fastq_minovlen "${FASTQ_MINOVLEN}" \
        --fastq_maxdiffs "${FASTQ_MAXDIFFS}" \
        --threads "${NSLOTS}" \
        2> "${MERGE_LOG}"

    PAIRS_IN=$(awk '/Pairs processed:/{print $NF}' "${MERGE_LOG}" \
               | tr -d ',' || echo "0")
    PAIRS_MERGED=$(awk '/Pairs merged:/{print $3}' "${MERGE_LOG}" \
                   | tr -d ',' || echo "0")
    if [[ "${PAIRS_IN:-0}" -gt 0 && -n "${PAIRS_MERGED:-}" ]]; then
        PCT=$(awk -v a="${PAIRS_MERGED}" -v b="${PAIRS_IN}" \
              'BEGIN {printf "%.2f", (a/b)*100}')
    else
        PCT="0.00"
    fi
    echo -e "${SAMPLE_NAME}\t${PAIRS_IN:-0}\t${PAIRS_MERGED:-0}\t${PCT}" \
        >> "${MERGE_SUMMARY}"

    # Filter merged reads by expected errors, convert to FASTA, label by sample
    FILTERED_FASTA="${OUTPUT_DIR}/03_filtered/${SAMPLE_NAME}_filtered.fasta"

    vsearch \
        --fastq_filter "${MERGED_FASTQ}" \
        --fastq_maxee "${FASTQ_MAXEE}" \
        --fastq_minlen "${MIN_LENGTH}" \
        --fastaout "${FILTERED_FASTA}" \
        --relabel "${SAMPLE_NAME}." \
        --fasta_width 0 \
        2> "${OUTPUT_DIR}/logs/${SAMPLE_NAME}_filter.log"

    cat "${FILTERED_FASTA}" >> "${ALL_FILTERED_FASTA}"
done

###############################################################################
# 11. Global dereplication
###############################################################################

log "Step 3/6: Dereplicating sequences globally..."

DEREP_FASTA="${OUTPUT_DIR}/04_derep/derep.fasta"

vsearch \
    --derep_fulllength "${ALL_FILTERED_FASTA}" \
    --output "${DEREP_FASTA}" \
    --uc "${OUTPUT_DIR}/04_derep/derep.uc" \
    --minuniquesize "${MIN_SIZE}" \
    --sizeout \
    --fasta_width 0 \
    2> "${OUTPUT_DIR}/logs/derep.log"

DEREP_SEQS=$(grep -c "^>" "${DEREP_FASTA}" || echo "0")
log "  Unique sequences after dereplication: ${DEREP_SEQS}"

###############################################################################
# 12. Chimera removal
###############################################################################

log "Step 4/6: Chimera removal (method: ${CHIMERA_METHOD})..."

NOCHIMERA_FASTA="${OUTPUT_DIR}/05_chimera/nochimeras.fasta"
CHIMERA_LOG="${OUTPUT_DIR}/logs/chimera.log"

if [[ "${CHIMERA_METHOD}" == "denovo" ]]; then
    vsearch \
        --uchime3_denovo "${DEREP_FASTA}" \
        --nonchimeras "${NOCHIMERA_FASTA}" \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        2> "${CHIMERA_LOG}"
else
    vsearch \
        --uchime_ref "${DEREP_FASTA}" \
        --db "${REF_DB}" \
        --nonchimeras "${NOCHIMERA_FASTA}" \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        2> "${CHIMERA_LOG}"
fi

NOCHIMERA_SEQS=$(grep -c "^>" "${NOCHIMERA_FASTA}" || echo "0")
log "  Sequences after chimera removal: ${NOCHIMERA_SEQS}"

###############################################################################
# 13. OTU clustering
###############################################################################

log "Step 5/6: Clustering OTUs at ${IDENTITY} identity..."

OTUS_FASTA="${OUTPUT_DIR}/06_otus/otus.fasta"

vsearch \
    --cluster_size "${NOCHIMERA_FASTA}" \
    --id "${IDENTITY}" \
    --centroids "${OTUS_FASTA}" \
    --sizein \
    --sizeout \
    --relabel "OTU_" \
    --fasta_width 0 \
    --threads "${NSLOTS}" \
    2> "${OUTPUT_DIR}/logs/cluster.log"

N_OTUS=$(grep -c "^>" "${OTUS_FASTA}" || echo "0")
log "  OTUs generated: ${N_OTUS}"

###############################################################################
# 14. Map reads to OTUs and generate OTU table
###############################################################################

log "Step 6/6: Mapping reads to OTUs..."

OTU_TABLE="${OUTPUT_DIR}/06_otus/otu_table.tsv"

# Map all merged+filtered reads (with sample labels) back to OTU centroids
vsearch \
    --usearch_global "${ALL_FILTERED_FASTA}" \
    --db "${OTUS_FASTA}" \
    --id "${IDENTITY}" \
    --otutabout "${OTU_TABLE}" \
    --strand plus \
    --threads "${NSLOTS}" \
    2> "${OUTPUT_DIR}/logs/mapping.log"

log "  OTU table: ${OTU_TABLE}"

###############################################################################
# 15. Generate summary report
###############################################################################

log "Generating summary report..."

SUMMARY_TXT="${OUTPUT_DIR}/summary_report.txt"

cat > "${SUMMARY_TXT}" << EOF
================================================================================
VSEARCH OTU Pipeline - Summary Report
================================================================================
Date: $(date)
Input directory:  ${INPUT_DIR}
Output directory: ${OUTPUT_DIR}
Samples:          ${#R1_FILES[@]}

Parameters:
-----------
  Pattern R1:           ${PATTERN_R1}
  Pattern R2:           ${PATTERN_R2}
  Threads:              ${NSLOTS}
  Min base quality:     ${MIN_QUALITY}
  Min read length:      ${MIN_LENGTH} bp
  Min merge overlap:    ${FASTQ_MINOVLEN} bp
  Max merge diffs:      ${FASTQ_MAXDIFFS}
  Max expected errors:  ${FASTQ_MAXEE}
  Min unique size:      ${MIN_SIZE}
  Chimera method:       ${CHIMERA_METHOD}
  OTU identity:         ${IDENTITY}
EOF

if [[ "${CHIMERA_METHOD}" == "ref" ]]; then
    echo "  Reference database:   ${REF_DB}" >> "${SUMMARY_TXT}"
fi

cat >> "${SUMMARY_TXT}" << EOF

Results:
--------
  Unique seqs after dereplication:  ${DEREP_SEQS}
  Seqs after chimera removal:       ${NOCHIMERA_SEQS}
  OTUs generated:                   ${N_OTUS}

Output files:
-------------
  Fastp trimmed reads:       ${OUTPUT_DIR}/01_fastp/
  Merged reads (FASTQ):      ${OUTPUT_DIR}/02_merged/
  Filtered FASTA per sample: ${OUTPUT_DIR}/03_filtered/
  Dereplicated sequences:    ${OUTPUT_DIR}/04_derep/derep.fasta
  Non-chimeric sequences:    ${OUTPUT_DIR}/05_chimera/nochimeras.fasta
  OTU representatives:       ${OUTPUT_DIR}/06_otus/otus.fasta
  OTU table:                 ${OUTPUT_DIR}/06_otus/otu_table.tsv
  Processing logs:           ${OUTPUT_DIR}/logs/
  Per-step statistics:       ${OUTPUT_DIR}/stats/

================================================================================
EOF

echo -e "\nFastp statistics:\n" >> "${SUMMARY_TXT}"
column -t -s $'\t' "${FASTP_SUMMARY}" >> "${SUMMARY_TXT}"

echo -e "\nMerging statistics:\n" >> "${SUMMARY_TXT}"
column -t -s $'\t' "${MERGE_SUMMARY}" >> "${SUMMARY_TXT}"

cat "${SUMMARY_TXT}"

###############################################################################
# 16. End
###############################################################################

log "${GREEN}2.2-vsearch_pipeline.sh completed successfully${NC}"
log "OTU table:  ${OTU_TABLE}"
log "OTU FASTA:  ${OTUS_FASTA}"
