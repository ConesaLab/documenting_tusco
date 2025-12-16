#!/bin/bash
#SBATCH --job-name=minimap2_mapping
#SBATCH --output=minimap2_mapping_%A_%a.out
#SBATCH --error=minimap2_mapping_%A_%a.err
#SBATCH --time=0-24:00:00  # Adjusted time, can be modified
#SBATCH --mem=64G         # Adjusted memory, can be modified
#SBATCH --qos=short        # Adjusted qos, can be modified
#SBATCH --cpus-per-task=10 # Default CPUs, adjust if minimap2 needs more/less
#SBATCH --array=0-9        # Placeholder, will be updated dynamically
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tianyuan.liu@csic.es

# --- Configuration ---
BASE_DATA_DIR="/storage/gge/Tian/lrgasp_data"
OUTPUT_BASE_DIR="${BASE_DATA_DIR}/map_result"
LOG_DIR="${BASE_DATA_DIR}/log" # For minimap2 logs if needed

# Reference files - Human
HUMAN_GENOME="/storage/gge/Tian/lrgasp_data/reference/lrgasp_grch38_sirvs.fasta.gz"

# Reference files - Mouse
MOUSE_GENOME="/storage/gge/Tian/lrgasp_data/reference/lrgasp_grcm39_sirvs.fasta.gz"

# Ensure output and log directories exist
mkdir -p "${OUTPUT_BASE_DIR}"
mkdir -p "${LOG_DIR}"

# --- Helper Functions ---
source_conda_env() {
    echo "Activating Conda environment..."
    # Source Conda initialization script
    CONDA_BASE_PATH=$(conda info --base)
    if [ -z "${CONDA_BASE_PATH}" ]; then
        echo "Error: Conda base path not found. Is Conda installed and configured?"
        exit 1
    fi
    source "${CONDA_BASE_PATH}/etc/profile.d/conda.sh"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to source conda.sh."
        exit 1
    fi

    CONDA_ENV_PATH="/home/tyuan/.conda/envs/minimap2"
    conda activate "${CONDA_ENV_PATH}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to activate conda environment ${CONDA_ENV_PATH}."
        exit 1
    fi

    command -v minimap2 >/dev/null 2>&1 || { echo "Error: minimap2 not found in Conda environment ${CONDA_ENV_PATH}. Exiting." >&2; exit 1; }
    command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found in Conda environment ${CONDA_ENV_PATH}. Exiting." >&2; exit 1; }
    
    echo "Successfully activated Conda environment: ${CONDA_ENV_PATH}"
    echo "Minimap2 path: $(command -v minimap2)"
    echo "Samtools path: $(command -v samtools)"
}

# --- Sample Directories ---
# Define the list of sample directories to process
# Excludes illumina, r2c2, reference, map_result, and log directories
SAMPLE_DIRS_FULL_PATH=()
while IFS= read -r line; do
    SAMPLE_DIRS_FULL_PATH+=("$line")
done < <(find "${BASE_DATA_DIR}" -mindepth 1 -maxdepth 1 -type d \
    ! -name '*illumina*' \
    ! -name '*r2c2*' \
    ! -name 'reference' \
    ! -name 'map_result' \
    ! -name 'log' \
    -print)

NUM_TOTAL_SAMPLES=${#SAMPLE_DIRS_FULL_PATH[@]}

if [ "${NUM_TOTAL_SAMPLES}" -eq 0 ]; then
    echo "Error: No sample directories found in ${BASE_DATA_DIR} matching the criteria."
    exit 1
fi

echo "Found ${NUM_TOTAL_SAMPLES} sample directories to process."

# Update SBATCH array directive based on the number of samples
# This is a common pattern but might require manual adjustment or a pre-flight script
# to update the script file itself if not submitting with dynamic array size.
# For now, we'll assume the user will set it correctly when submitting.
# Example: #SBATCH --array=0-$((NUM_TOTAL_SAMPLES - 1))
echo "INFO: Ensure your #SBATCH --array directive is set correctly, e.g., --array=0-$((NUM_TOTAL_SAMPLES - 1))"

# --- Main Script ---
echo "Starting minimap2 mapping pipeline script"
echo "Date: $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"

# Check if SLURM_ARRAY_TASK_ID is set and within bounds
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID is not set. This script is designed to be run as a job array."
    exit 1
fi

if [ "$SLURM_ARRAY_TASK_ID" -ge "$NUM_TOTAL_SAMPLES" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is out of bounds. Max index is $(($NUM_TOTAL_SAMPLES - 1))."
    exit 1
fi

# Get the specific sample directory for this array task
CURRENT_SAMPLE_PATH=${SAMPLE_DIRS_FULL_PATH[$SLURM_ARRAY_TASK_ID]}
SAMPLE_DIR_NAME=$(basename "${CURRENT_SAMPLE_PATH}")

echo "--------------------------------------------------"
echo "Processing sample directory: ${SAMPLE_DIR_NAME} (Task ID: $SLURM_ARRAY_TASK_ID)"
echo "Full path: ${CURRENT_SAMPLE_PATH}"
echo "--------------------------------------------------"

# Create output directory for the current sample
CURRENT_OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE_DIR_NAME}"
mkdir -p "${CURRENT_OUTPUT_DIR}"

cd "${CURRENT_OUTPUT_DIR}" || { echo "Failed to cd to ${CURRENT_OUTPUT_DIR}"; exit 1; }

# Load modules
source_conda_env

# Determine species and set reference genome
GENOME_FASTA=""
if [[ "${SAMPLE_DIR_NAME}" == es_* ]]; then
    SPECIES="mouse"
    GENOME_FASTA="${MOUSE_GENOME}"
    echo "Species detected: Mouse. Using reference: ${GENOME_FASTA}"
elif [[ "${SAMPLE_DIR_NAME}" == wtc11_* ]]; then
    SPECIES="human"
    GENOME_FASTA="${HUMAN_GENOME}"
    echo "Species detected: Human. Using reference: ${GENOME_FASTA}"
else
    echo "Error: Could not determine species for ${SAMPLE_DIR_NAME}. Exiting."
    exit 1
fi

if [ ! -f "${GENOME_FASTA}" ]; then
    echo "Error: Reference genome ${GENOME_FASTA} not found for ${SPECIES}. Exiting."
    exit 1
fi

# Determine data type and minimap2 parameters
# IMPORTANT: No reference annotation is provided to minimap2 during the mapping stage.
# This is intentional to avoid annotation-dependent alignment bias, which could impact
# downstream transcript reconstruction and isoform discovery. By using only the genome
# reference (without splice junction annotations), we ensure that read alignments are
# based solely on sequence similarity and splice site signals, not on existing gene models.
# This approach tests the true capability of downstream tools to perform data-driven
# isoform discovery, which is critical for the TUSCO-novel evaluation framework.
MINIMAP2_PARAMS=""
MM2_THREADS=$SLURM_CPUS_PER_TASK # Default to Slurm CPUs, adjust if needed

if [[ "${SAMPLE_DIR_NAME}" == *pacbio* ]]; then
    echo "Data type: PacBio"
    # PacBio, minimap2 -ax splice:hq -uf --MD -t 40
    MM2_THREADS=40 # As per user request
    MINIMAP2_PARAMS="-ax splice:hq -uf --MD -t ${MM2_THREADS}"
elif [[ "${SAMPLE_DIR_NAME}" == *drna_ont* ]]; then
    echo "Data type: ONT direct RNA-seq"
    # ONT direct RNA-seq, minimap2 -ax splice -uf -k14 --MD -t 30
    MM2_THREADS=30 # As per user request
    MINIMAP2_PARAMS="-ax splice -uf -k14 --MD -t ${MM2_THREADS}"
elif [[ "${SAMPLE_DIR_NAME}" == *ont* ]]; then # Covers cdna_ont and captrap_ont
    echo "Data type: ONT cDNA"
    # ONT, minimap2 -ax splice --MD -t 30
    MM2_THREADS=30 # As per user request
    MINIMAP2_PARAMS="-ax splice --MD -t ${MM2_THREADS}"
else
    echo "Error: Could not determine data type for ${SAMPLE_DIR_NAME} (pacbio, drna_ont, or ont not in name). Exiting."
    exit 1
fi

echo "Using minimap2 parameters: ${MINIMAP2_PARAMS}"
echo "Minimap2 threads: ${MM2_THREADS} (Requesting ${SLURM_CPUS_PER_TASK} CPUs from Slurm)"
if [ "${MM2_THREADS}" -gt "${SLURM_CPUS_PER_TASK}" ]; then
    echo "Warning: Minimap2 threads (${MM2_THREADS}) exceed allocated Slurm CPUs (${SLURM_CPUS_PER_TASK}). Performance may be affected or job may fail."
    echo "Consider adjusting --cpus-per-task in SBATCH directives to at least ${MM2_THREADS}."
fi


# Find FASTQ files (handles .fastq.gz, .fq.gz, .fastq, .fq)
# Using an array to properly handle spaces in filenames if any (though unlikely for FASTQ)
mapfile -t FASTQ_FILES < <(find "${CURRENT_SAMPLE_PATH}" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \))

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "Error: No FASTQ files found in ${CURRENT_SAMPLE_PATH} for ${SAMPLE_DIR_NAME}. Exiting."
    exit 1
fi

echo "Found ${#FASTQ_FILES[@]} FASTQ file(s) for ${SAMPLE_DIR_NAME}:"
for fq_file in "${FASTQ_FILES[@]}"; do
    echo "  - ${fq_file}"
done

# Define output file names
OUTPUT_SAM="${CURRENT_OUTPUT_DIR}/${SAMPLE_DIR_NAME}.sam"
OUTPUT_BAM="${CURRENT_OUTPUT_DIR}/${SAMPLE_DIR_NAME}.bam"
OUTPUT_SORTED_BAM="${CURRENT_OUTPUT_DIR}/${SAMPLE_DIR_NAME}.sorted.bam"
MINIMAP2_LOG="${LOG_DIR}/minimap2_${SAMPLE_DIR_NAME}_${SLURM_ARRAY_TASK_ID}.log"

echo "Output SAM: ${OUTPUT_SAM}"
echo "Output BAM (sorted): ${OUTPUT_SORTED_BAM}"
echo "Minimap2 Log: ${MINIMAP2_LOG}"

# Run minimap2
echo "Starting minimap2 alignment for ${SAMPLE_DIR_NAME}..."
# minimap2 outputs to stdout, so redirect to SAM file. stderr for logging.
minimap2 ${MINIMAP2_PARAMS} "${GENOME_FASTA}" "${FASTQ_FILES[@]}" > "${OUTPUT_SAM}" 2> "${MINIMAP2_LOG}"

if [ $? -ne 0 ]; then
    echo "Error: minimap2 failed for ${SAMPLE_DIR_NAME}. Check log: ${MINIMAP2_LOG}"
    exit 1
fi

if [ ! -s "${OUTPUT_SAM}" ]; then
    echo "Error: minimap2 output SAM file ${OUTPUT_SAM} is empty for ${SAMPLE_DIR_NAME}. Check log: ${MINIMAP2_LOG}"
    exit 1
fi
echo "Minimap2 alignment completed successfully for ${SAMPLE_DIR_NAME}."

# Convert SAM to BAM, sort, and index
echo "Converting SAM to BAM, sorting, and indexing for ${SAMPLE_DIR_NAME}..."
SAMTOOLS_THREADS=$((SLURM_CPUS_PER_TASK > 1 ? SLURM_CPUS_PER_TASK - 1 : 1)) # Use N-1 cores for samtools sort

samtools view -@ "${SAMTOOLS_THREADS}" -bS "${OUTPUT_SAM}" -o "${OUTPUT_BAM}"
if [ $? -ne 0 ]; then echo "Error: samtools view failed for ${OUTPUT_SAM}"; exit 1; fi

samtools sort -@ "${SAMTOOLS_THREADS}" -o "${OUTPUT_SORTED_BAM}" "${OUTPUT_BAM}"
if [ $? -ne 0 ]; then echo "Error: samtools sort failed for ${OUTPUT_BAM}"; exit 1; fi

samtools index "${OUTPUT_SORTED_BAM}"
if [ $? -ne 0 ]; then echo "Error: samtools index failed for ${OUTPUT_SORTED_BAM}"; exit 1; fi

echo "BAM conversion, sorting, and indexing completed for ${SAMPLE_DIR_NAME}."

# Optional: Remove intermediate SAM and unsorted BAM to save space
echo "Cleaning up intermediate files..."
rm "${OUTPUT_SAM}"
rm "${OUTPUT_BAM}"
echo "Intermediate files removed."


echo "--------------------------------------------------"
echo "Task $SLURM_ARRAY_TASK_ID for sample ${SAMPLE_DIR_NAME} completed."
echo "Output sorted BAM: ${OUTPUT_SORTED_BAM}"
echo "Pipeline task finished at $(date)"
echo "--------------------------------------------------"

exit 0 