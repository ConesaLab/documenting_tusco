#!/bin/bash

#SBATCH --job-name=sqanti3_RIN_sequin
#SBATCH --output=sqanti3_RIN_sequin_%A_%a.out
#SBATCH --error=sqanti3_RIN_sequin_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --qos=short
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tianyuan.liu@csic.es
#SBATCH --array=0-17

# Load conda environment
source /home/tyuan/miniconda3/etc/profile.d/conda.sh
conda activate SQANTI3.env

# Variables
ANNOTATION="/storage/gge/genomes/sequin_standards/rnasequin_annotation_2.4.gtf"
GENOME="/storage/gge/genomes/sequin_standards/rnasequin_decoychr_2.4.fa"
GTF_LIST="/home/tyuan/RIN_sequin/gtf_files.txt"
WORK_DIR="/home/tyuan/RIN_sequin"
THREADS=10

# Change to working directory
cd "$WORK_DIR" || { echo "Cannot change to WORK_DIR: $WORK_DIR"; exit 1; }

# Read GTF file list into array
mapfile -t GTF_ARRAY < "$GTF_LIST"

# Get current GTF file based on SLURM_ARRAY_TASK_ID
CURRENT_GTF=${GTF_ARRAY[$SLURM_ARRAY_TASK_ID]}

# Check if CURRENT_GTF is set
if [ -z "$CURRENT_GTF" ]; then
    echo "No GTF file assigned to array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Extract sample ID from GTF file name
SAMPLE_ID=$(basename "$CURRENT_GTF" .collapsed_isoforms.gtf)

# Define output directory for SQANTI3
OUTPUT_DIR="$WORK_DIR/SQANTI3_outputs/$SAMPLE_ID"
mkdir -p "$OUTPUT_DIR"

# Run SQANTI3
echo "[SQANTI3] Running SQANTI3 for $SAMPLE_ID"
/home/tyuan/GitHub/SQANTI3/sqanti3_qc.py "$CURRENT_GTF" \
             "$ANNOTATION" \
             "$GENOME" \
             -o "$SAMPLE_ID" \
             -d "$OUTPUT_DIR" \
             --skipORF \
             -n 10 \
             -t 10

# Verify SQANTI3 Output
if [ $? -ne 0 ]; then
    echo "[SQANTI3] SQANTI3 failed for $SAMPLE_ID"
    exit 1
fi

echo "[SQANTI3] Completed SQANTI3 for $SAMPLE_ID"

