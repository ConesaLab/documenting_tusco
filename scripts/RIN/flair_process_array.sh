#!/bin/bash

#SBATCH --job-name=flair_RIN_sequin
#SBATCH --output=flair_RIN_sequin_%A_%a.out
#SBATCH --error=flair_RIN_sequin_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --qos=short
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tianyuan.liu@csic.es
#SBATCH --array=0-17

# Load environment and modules
source /home/tyuan/miniconda3/etc/profile.d/conda.sh
conda activate flair

# Set variables
GENOME="/storage/gge/genomes/sequin_standards/rnasequin_decoychr_2.4.fa"
ANNOTATION="/storage/gge/genomes/sequin_standards/rnasequin_annotation_2.4.gtf"
WORK_DIR="/home/tyuan/RIN_sequin"
FASTQ_LIST="$WORK_DIR/fastq_list.txt"
THREADS=12

# Change to working directory
cd "$WORK_DIR" || { echo "Cannot change to WORK_DIR: $WORK_DIR"; exit 1; }

# Read FASTQ file list into array
mapfile -t FASTQS < "$FASTQ_LIST"

# Get current FASTQ file based on SLURM_ARRAY_TASK_ID
CURRENT_FASTQ=${FASTQS[$SLURM_ARRAY_TASK_ID]}

# Check if CURRENT_FASTQ is set
if [ -z "$CURRENT_FASTQ" ]; then
    echo "No FASTQ file assigned to array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Extract sample ID from FASTQ file name
SAMPLE_ID=$(basename "$CURRENT_FASTQ" .fastq)
SAMPLE_ID=$(basename "$SAMPLE_ID" .fq)
SAMPLE_ID=$(basename "$SAMPLE_ID" .pass)
SAMPLE_ID=$(basename "$SAMPLE_ID" .PCR_cDNA.pass)

# Create output directories
OUT_ALIGN="$WORK_DIR/alignments"
OUT_COLLAPSE="$WORK_DIR/collapse"
mkdir -p "$OUT_ALIGN"
mkdir -p "$OUT_COLLAPSE"

# Define output prefixes
ALIGN_PREFIX="$OUT_ALIGN/$SAMPLE_ID.aligned"
CORRECT_PREFIX="$OUT_ALIGN/$SAMPLE_ID.corrected"
COLLAPSE_PREFIX="$OUT_COLLAPSE/$SAMPLE_ID.collapsed"

# Run FLAIR Align
echo "[FLAIR] Aligning $SAMPLE_ID"
flair align \
    --genome "$GENOME" \
    --reads "$CURRENT_FASTQ" \
    --output "$ALIGN_PREFIX" \
    --threads "$THREADS" 

# Run FLAIR Correct
echo "[FLAIR] Correcting $SAMPLE_ID"
flair correct \
    --query "${ALIGN_PREFIX}.bed" \
    --genome "$GENOME" \
    --gtf "$ANNOTATION" \
    --output "$CORRECT_PREFIX" \
    --threads "$THREADS" 

# Run FLAIR Collapse
echo "[FLAIR] Collapsing $SAMPLE_ID"
flair collapse \
    -g  "$GENOME" \
    --gtf "$ANNOTATION" \
    -q "${CORRECT_PREFIX}_all_corrected.bed" \
    -r "$CURRENT_FASTQ" \
    --output "$COLLAPSE_PREFIX" \
    --check_splice \
    --threads "$THREADS" 

echo "[FLAIR] Completed processing for $SAMPLE_ID"

