#!/bin/bash
#SBATCH --job-name=flair_sqanti_pipeline
#SBATCH --output=flair_sqanti_pipeline_%A_%a.out
#SBATCH --error=flair_sqanti_pipeline_%A_%a.err
#SBATCH --time=0-24:00:00 # Adjusted time for short QOS (24 hours)
#SBATCH --mem=64G      # Adjusted memory per array task, can be modified
#SBATCH --partition=global   # Submit to the global partition where 'short' QoS is available
#SBATCH --qos short       # Adjusted qos, can be modified
#SBATCH --cpus-per-task=20 # Adjusted cpus, can be modified
#SBATCH --array=0-6 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tianyuan.liu@csic.es # Please replace with your email

# --- Configuration ---
# Removed CONDA_BASE_PATH, using module system instead
FLAIR_ENV="flair" # Corrected FLAIR environment name
SQANTI3_ENV="SQANTI3.env"
SQANTI3_PATH="${SQANTI3_PATH:-}" # optional: directory containing sqanti3_qc.py (leave empty if in PATH)
BASE_DATA_DIR="/storage/gge/Tian/lrgasp_data"
OUTPUT_BASE_DIR="/storage/gge/Tian/lrgasp_analysis/flair_sq3" # Updated output directory

# Reference files - Human
HUMAN_GENOME="/storage/gge/Tian/lrgasp_data/reference/lrgasp_grch38_sirvs.fasta"
HUMAN_TUSCO_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_v39_annotation_sirvs.human.full_with_fake.gtf"
HUMAN_REF_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_v39_annotation_sirvs.human.gtf"

# Reference files - Mouse
MOUSE_GENOME="/storage/gge/Tian/lrgasp_data/reference/lrgasp_grcm39_sirvs.fasta"
MOUSE_TUSCO_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_vM28_sirvs.mouse.tusco_sim_sorted.gtf"
MOUSE_REF_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_vM28_sirvs.mouse.gtf"

# --- SQANTI3 auxiliary annotation files ---
# Splice‑junction coverage (comma‑separated lists expected by SQANTI3 -c/--coverage)
HUMAN_COVERAGE_FILES="/storage/gge/Tian/lrgasp_data/reference/wtc11_rep1SJ.out.tab,/storage/gge/Tian/lrgasp_data/reference/wtc11_rep2SJ.out.tab,/storage/gge/Tian/lrgasp_data/reference/wtc11_rep3SJ.out.tab"
MOUSE_COVERAGE_FILES="/storage/gge/Tian/lrgasp_data/reference/es_rep1SJ.out.tab,/storage/gge/Tian/lrgasp_data/reference/es_rep2SJ.out.tab,/storage/gge/Tian/lrgasp_data/reference/es_rep3SJ.out.tab"

# PolyA‑site peak BEDs
HUMAN_POLYA_PEAK="/storage/gge/Tian/lrgasp_data/reference/QuantSeq_WTC11.all_reps.bed"
MOUSE_POLYA_PEAK="/storage/gge/Tian/lrgasp_data/reference/QuantSeq_ES.all_reps.bed"

# CAGE peak BEDs
HUMAN_CAGE_PEAK="/storage/gge/Tian/lrgasp_data/reference/human.refTSS_v3.1.hg38.bed"
MOUSE_CAGE_PEAK="/storage/gge/Tian/lrgasp_data/reference/mouse.refTSS_v3.1.mm39.bed"

# Short read junction files - Human
HUMAN_SJ_MERGED="/storage/gge/Tian/lrgasp_data/reference/wtc11_merged.sorted.SJ.out.tab"

# Short read junction files - Mouse
MOUSE_SJ_MERGED="/storage/gge/Tian/lrgasp_data/reference/es_merged.sorted.SJ.out.tab"


THREADS=${SLURM_CPUS_PER_TASK}

# --- Helper Functions ---
source_conda() {
    echo "Loading Anaconda module..."
    module load anaconda
    if [ $? -ne 0 ]; then
        echo "Error: Failed to load anaconda module."
        exit 1
    fi
    echo "Successfully loaded anaconda module"
}

activate_env() {
    local env_name=$1
    local tool_to_check=$2 # e.g., "flair" or "sqanti3_qc.py"
    local tool_path=$3     # Full path to the tool if provided
    
    echo "Activating conda environment: ${env_name}"
    
    echo "PATH before activating ${env_name}: ${PATH}"
    conda activate "${env_name}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to activate conda environment ${env_name}."
        exit 1
    fi
    echo "Successfully activated ${env_name}"
    echo "PATH after activating ${env_name}: ${PATH}"
    echo "CONDA_PREFIX: ${CONDA_PREFIX}"

    if [ -n "$tool_to_check" ]; then
        if [ -n "$tool_path" ]; then
            # Check specific path
            echo "Verifying existence of '$tool_to_check' at path: $tool_path/$tool_to_check"
            if [ -f "$tool_path/$tool_to_check" ]; then
                echo "'$tool_to_check' found at: $tool_path/$tool_to_check"
                # Add tool path to PATH for convenience
                export PATH="$tool_path:$PATH"
                echo "Added $tool_path to PATH"
            else
                echo "Error: $tool_path/$tool_to_check not found"
                exit 1
            fi
        else
            # Check in PATH
            echo "Verifying command '$tool_to_check' in PATH..."
            command -v "$tool_to_check"
            if [ $? -ne 0 ]; then
                echo "Error: command '$tool_to_check' not found in PATH after activating ${env_name}."
                if [ -d "${CONDA_PREFIX}/bin" ]; then
                    echo "Contents of ${CONDA_PREFIX}/bin:"
                    ls -l "${CONDA_PREFIX}/bin"
                else
                    echo "${CONDA_PREFIX}/bin directory not found."
                fi
                exit 1
            else
                echo "'$tool_to_check' is available at: $(command -v "$tool_to_check")"
            fi
        fi
    fi
}

# --- FLAIR Pipeline Functions ---
# Split into separate align and process functions for better efficiency
run_flair_align() {
    local sample_dir_name=$1
    local output_dir=$2
    
    echo "[FLAIR ALIGN] Starting for ${sample_dir_name}" >&2
    
    mkdir -p "${output_dir}/shared/flair_tmp"
    
    local flair_align_prefix="${output_dir}/shared/${sample_dir_name}"
    local bed_file="${flair_align_prefix}.bed"
    
    # Step 1: FLAIR align - only run if .bed file doesn't exist
    if [ ! -f "${bed_file}" ]; then
        echo "[FLAIR ALIGN] Running for ${sample_dir_name} (output doesn't exist yet)" >&2
        echo "Checking input reads for FLAIR align:" >&2
        for r_file in "${READS_FOR_FLAIR[@]}"; do
            if [ ! -f "$r_file" ]; then
                echo "Error: Read file $r_file not found for FLAIR align! Exiting task at line $LINENO." >&2; exit 1;
            else
                echo "Read file $r_file found." >&2
            fi
        done
        if [ ! -f "${GENOME_FASTA}" ]; then
            echo "Error: Genome FASTA ${GENOME_FASTA} not found for FLAIR align! Exiting task at line $LINENO." >&2; exit 1;
        fi

        echo "Starting align..." >&2
        # Redirect FLAIR align stdout to stderr to avoid contaminating command substitution output
        flair align \
            -g "${GENOME_FASTA}" \
            -r "${READS_FOR_FLAIR[@]}" \
            -t "${THREADS}" \
            -o "${flair_align_prefix}" \
            1>&2

        # Wait a moment to ensure file system has updated
        sleep 5
        
        # Explicitly check that the bed file exists after alignment
        if [ ! -f "${bed_file}" ]; then
            echo "Error: FLAIR align output (${bed_file}) not found for ${sample_dir_name}. Exiting." >&2
            if [ -d "${output_dir}/shared/flair_tmp" ]; then
                ls -lhtr "${output_dir}/shared/flair_tmp" >&2
            fi
            echo "Error encountered (FLAIR align output not found). Exiting task at line $LINENO." >&2; exit 1;
        fi
        echo "FLAIR align completed successfully, output: ${bed_file}" >&2
    else
        echo "[FLAIR ALIGN] Skipping - output ${bed_file} already exists" >&2
    fi
    
    # Confirm bed file exists before returning
    if [ ! -f "${bed_file}" ]; then
        echo "Critical error: BED file ${bed_file} not found after align step. Exiting task at line $LINENO." >&2
        exit 1
    fi
    
    echo "${bed_file}" # Return the bed file path
}

run_flair_correct_collapse() {
    local sample_dir_name=$1
    local output_dir=$2
    local bed_file=$3
    local gtf_file=$4
    local gtf_type=$5  # "TUSCO" or "REF"
    local species_arg=$6 # "human" or "mouse"
    
    echo "[FLAIR CORRECT/COLLAPSE] Starting for ${sample_dir_name} with ${gtf_type} GTF in ${output_dir}" >&2
    
    mkdir -p "${output_dir}/flair_tmp"
    
    local flair_correct_prefix="${output_dir}/${sample_dir_name}" 
    local flair_collapse_prefix="${output_dir}/${sample_dir_name}"
    
    local corrected_bed="${flair_correct_prefix}_all_corrected.bed"
    local isoforms_gtf="${flair_collapse_prefix}.isoforms.gtf"
    
    # Step 2: FLAIR correct - only run if corrected BED file doesn't exist
    if [ ! -f "${corrected_bed}" ]; then
        echo "[FLAIR CORRECT] Running for ${sample_dir_name} with ${gtf_type} GTF (output doesn't exist yet)" >&2
        
        # Double-check bed file exists and is not empty before continuing
        if [ ! -f "${bed_file}" ]; then
            echo "Error: BED file for FLAIR correct (${bed_file}) not found! Exiting task at line $LINENO." >&2
            exit 1
        fi
        
        if [ ! -s "${bed_file}" ]; then
            echo "Error: BED file ${bed_file} exists but is empty! Exiting task at line $LINENO." >&2
            exit 1
        fi
        
        echo "Using BED file for FLAIR correct: ${bed_file} ($(du -h "${bed_file}" | cut -f1))" >&2
        
        if [ ! -f "${GENOME_FASTA}" ]; then
            echo "Error: Genome FASTA ${GENOME_FASTA} not found for FLAIR correct! Exiting task at line $LINENO." >&2; exit 1;
        fi
        if [ ! -f "${gtf_file}" ]; then
            echo "Error: GTF ${gtf_file} not found for FLAIR correct! Exiting task at line $LINENO." >&2; exit 1;
        fi

        echo "Starting correct..." >&2
        # Redirect FLAIR correct stdout to stderr
        local flair_correct_cmd=("flair" "correct" \
            "-q" "${bed_file}" \
            "-g" "${GENOME_FASTA}" \
            "-f" "${gtf_file}" \
            "-t" "${THREADS}" \
            "-o" "${flair_correct_prefix}")

        if [ "${species_arg}" == "human" ]; then
            flair_correct_cmd+=("--shortread" "${HUMAN_SJ_MERGED}")
            echo "Adding merged human short read SJ file to FLAIR correct: ${HUMAN_SJ_MERGED}" >&2
        elif [ "${species_arg}" == "mouse" ]; then
            flair_correct_cmd+=("--shortread" "${MOUSE_SJ_MERGED}")
            echo "Adding merged mouse short read SJ file to FLAIR correct: ${MOUSE_SJ_MERGED}" >&2
        fi

        echo "Executing FLAIR correct command: ${flair_correct_cmd[@]}" >&2
        "${flair_correct_cmd[@]}" 1>&2

        if [ ! -f "${corrected_bed}" ]; then
            echo "Error: FLAIR correct output (${corrected_bed}) not found for ${sample_dir_name}." >&2
            echo "Error encountered (FLAIR correct output not found). Exiting task at line $LINENO." >&2; exit 1;
        fi
        echo "FLAIR correct completed successfully, output: ${corrected_bed}" >&2
    else
        echo "[FLAIR CORRECT] Skipping - output ${corrected_bed} already exists" >&2
    fi

    # Step 3: FLAIR collapse - only run if isoforms GTF file doesn't exist
    if [ ! -f "${isoforms_gtf}" ]; then
        echo "[FLAIR COLLAPSE] Running for ${sample_dir_name} with ${gtf_type} GTF (output doesn't exist yet)" >&2
        echo "Checking input reads for FLAIR collapse:" >&2
        for r_file in "${READS_FOR_FLAIR[@]}"; do
            if [ ! -f "$r_file" ]; then
                echo "Error: Read file $r_file not found for FLAIR collapse! Exiting task at line $LINENO." >&2; exit 1;
            else
                echo "Read file $r_file found." >&2
            fi
        done
        if [ ! -f "${corrected_bed}" ]; then
            echo "Error: Corrected BED file ${corrected_bed} not found for FLAIR collapse! Exiting task at line $LINENO." >&2; exit 1;
        fi
        if [ ! -f "${GENOME_FASTA}" ]; then
            echo "Error: Genome FASTA ${GENOME_FASTA} not found for FLAIR collapse! Exiting task at line $LINENO." >&2; exit 1;
        fi
        if [ ! -f "${gtf_file}" ]; then
            echo "Error: GTF ${gtf_file} not found for FLAIR collapse! Exiting task at line $LINENO." >&2; exit 1;
        fi

        if [ "${gtf_type}" == "REF" ]; then
            COLLAPSE_OPTS=(--stringent --check_splice --annotation_reliant generate --generate_map)
        else
            COLLAPSE_OPTS=(--stringent --check_splice)
        fi

        # Redirect FLAIR collapse stdout to stderr
        flair collapse \
            -g "${GENOME_FASTA}" \
            -r "${READS_FOR_FLAIR[@]}" \
            -q "${corrected_bed}" \
            -f "${gtf_file}" \
            -t "${THREADS}" \
            "${COLLAPSE_OPTS[@]}" \
            -o "${flair_collapse_prefix}" \
            1>&2

        if [ ! -f "${isoforms_gtf}" ]; then
            echo "Error: FLAIR collapse output GTF (${isoforms_gtf}) not found for ${sample_dir_name}." >&2
            echo "Error encountered (FLAIR collapse output not found). Exiting task at line $LINENO." >&2; exit 1;
        fi
        echo "FLAIR collapse completed successfully, output: ${isoforms_gtf}" >&2
    else
        echo "[FLAIR COLLAPSE] Skipping - output ${isoforms_gtf} already exists" >&2
    fi


    echo "[FLAIR] Completed for ${sample_dir_name} with ${gtf_type} GTF. Output GTF: ${isoforms_gtf}" >&2
    
    echo "${isoforms_gtf}" # Return the isoforms GTF path
}

# --- SQANTI3 QC Function ---
run_sqanti3_qc() {
    local sample_dir_name=$1
    local output_dir=$2
    local isoforms_gtf=$3
    local ref_gtf=$4
    local gtf_type=$5  # "TUSCO" or "REF"
    
    local sqanti3_output_dir="${output_dir}/sqanti3_out"
    local sqanti3_out_prefix="${sample_dir_name}_sqanti"
    local sqanti3_result="${sqanti3_output_dir}/${sqanti3_out_prefix}_classification.txt"
    
    # Always rerun SQANTI3: remove any existing output directory then recreate it
    if [ -d "${sqanti3_output_dir}" ]; then
        echo "[SQANTI3] Removing existing output directory ${sqanti3_output_dir} to force fresh run"
        rm -rf "${sqanti3_output_dir}"
    fi
    mkdir -p "${sqanti3_output_dir}"

    echo "[SQANTI3] Starting for ${sample_dir_name} with ${gtf_type} GTF (output doesn't exist yet)"
    source_conda
    activate_env "${SQANTI3_ENV}" "sqanti3_qc.py" "${SQANTI3_PATH}"

    # Check if reference genome for SQANTI3 needs to be unzipped
    local ref_genome_sqanti3="${GENOME_FASTA}"
    if [[ "${GENOME_FASTA}" == *.gz ]]; then
        echo "Reference genome for SQANTI3 is gzipped. Checking if unzipped version exists or unzipping..."
        local unzipped_genome_fasta="${output_dir}/$(basename "${GENOME_FASTA}" .gz)"
        if [ ! -f "${unzipped_genome_fasta}" ]; then
            echo "Unzipping ${GENOME_FASTA} to ${unzipped_genome_fasta}"
            echo "Permissions for ${output_dir} before unzipping genome:"
            ls -ld "${output_dir}"
            gunzip -c "${GENOME_FASTA}" > "${unzipped_genome_fasta}"
            if [ $? -ne 0 ]; then
                echo "Error: Failed to unzip genome ${GENOME_FASTA}. Exiting task at line $LINENO."; exit 1;
            fi
            ls -lh "${unzipped_genome_fasta}"
        else
            echo "Using existing unzipped genome: ${unzipped_genome_fasta}"
        fi
        ref_genome_sqanti3="${unzipped_genome_fasta}"
    fi

    # Check if reference GTF for SQANTI3 needs to be unzipped
    local ref_gtf_sqanti3="${ref_gtf}"
    if [[ "${ref_gtf}" == *.gz ]]; then
        echo "Reference GTF for SQANTI3 is gzipped. Checking if unzipped version exists or unzipping..."
        local unzipped_ref_gtf="${output_dir}/$(basename "${ref_gtf}" .gz)"
        if [ ! -f "${unzipped_ref_gtf}" ]; then
            echo "Unzipping ${ref_gtf} to ${unzipped_ref_gtf}"
            echo "Permissions for ${output_dir} before unzipping reference GTF:"
            ls -ld "${output_dir}"
            gunzip -c "${ref_gtf}" > "${unzipped_ref_gtf}"
            if [ $? -ne 0 ]; then
                echo "Error: Failed to unzip reference GTF ${ref_gtf}. Exiting task at line $LINENO."; exit 1;
            fi
            ls -lh "${unzipped_ref_gtf}"
        else
            echo "Using existing unzipped reference GTF: ${unzipped_ref_gtf}"
        fi
        ref_gtf_sqanti3="${unzipped_ref_gtf}"
    fi

    # --- Select species‑specific auxiliary files for SQANTI3 ---
    if [ "${SPECIES}" == "human" ]; then
        SQANTI3_COVERAGE="${HUMAN_COVERAGE_FILES}"
        SQANTI3_POLYA="${HUMAN_POLYA_PEAK}"
        SQANTI3_CAGE="${HUMAN_CAGE_PEAK}"
    else
        SQANTI3_COVERAGE="${MOUSE_COVERAGE_FILES}"
        SQANTI3_POLYA="${MOUSE_POLYA_PEAK}"
        SQANTI3_CAGE="${MOUSE_CAGE_PEAK}"
    fi

    # Verify auxiliary annotation files exist
    for sj in $(echo "${SQANTI3_COVERAGE}" | tr ',' ' '); do
        if [ ! -f "${sj}" ]; then
            echo "Error: Coverage file ${sj} not found for SQANTI3! Exiting task at line $LINENO." >&2; exit 1;
        fi
    done
    if [ ! -f "${SQANTI3_POLYA}" ]; then
        echo "Error: polyA peak file ${SQANTI3_POLYA} not found for SQANTI3! Exiting task at line $LINENO." >&2; exit 1;
    fi
    if [ ! -f "${SQANTI3_CAGE}" ]; then
        echo "Error: CAGE peak file ${SQANTI3_CAGE} not found for SQANTI3! Exiting task at line $LINENO." >&2; exit 1;
    fi

    echo "[SQANTI3 QC] Running for ${sample_dir_name} with ${gtf_type} GTF"
    if [ ! -f "${isoforms_gtf}" ]; then
        echo "Error: FLAIR result GTF ${isoforms_gtf} not found for SQANTI3! Exiting task at line $LINENO."; exit 1;
    fi
    if [ ! -f "${ref_gtf_sqanti3}" ]; then
        echo "Error: Reference GTF ${ref_gtf_sqanti3} not found for SQANTI3! Exiting task at line $LINENO."; exit 1;
    fi
    if [ ! -f "${ref_genome_sqanti3}" ]; then
        echo "Error: Reference Genome FASTA ${ref_genome_sqanti3} not found for SQANTI3! Exiting task at line $LINENO."; exit 1;
    fi

    "${SQANTI3_PATH}/sqanti3_qc.py" \
        "${isoforms_gtf}" \
        "${ref_gtf_sqanti3}" \
        "${ref_genome_sqanti3}" \
        -d "${sqanti3_output_dir}" \
        -o "${sqanti3_out_prefix}" \
        -c "${SQANTI3_COVERAGE}" \
        --polyA_peak "${SQANTI3_POLYA}" \
        --CAGE_peak "${SQANTI3_CAGE}" \
        --skipORF \
        -t "${THREADS}" \
        --report pdf

    if [ $? -ne 0 ]; then
        echo "Error: SQANTI3 QC failed for ${sample_dir_name} with ${gtf_type} GTF."
        echo "Error encountered (SQANTI3 QC failed). Exiting task at line $LINENO."; exit 1;
    else
        echo "[SQANTI3] Completed for ${sample_dir_name} with ${gtf_type} GTF. Results in ${sqanti3_output_dir}"
    fi

    return 0
}

# --- Main Script ---
echo "Starting FLAIR and SQANTI3 pipeline script"
echo "Date: $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"

mkdir -p "${OUTPUT_BASE_DIR}"
cd "${OUTPUT_BASE_DIR}" || { echo "Failed to cd to ${OUTPUT_BASE_DIR}"; exit 1; }

# Iterate over relevant directories in lrgasp_data
SAMPLE_DIRS=()
while IFS= read -r line; do
    SAMPLE_DIRS+=("$line")
done < <(ls -d ${BASE_DATA_DIR}/*/ | grep -v 'r2c2' | grep -v 'reference' | grep -v 'log' | grep -v 'illumina')


echo "Found ${#SAMPLE_DIRS[@]} total sample directories eligible for processing."
NUM_TOTAL_SAMPLES=${#SAMPLE_DIRS[@]}
echo "This script is task $SLURM_ARRAY_TASK_ID of a $SLURM_ARRAY_JOB_COUNT task array."

# Check if SLURM_ARRAY_TASK_ID is set and within bounds
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID is not set. This script is designed to be run as a job array."
    echo "Please determine the number of samples (M) and submit with --array=0-(M-1)%<concurrency>."
    echo "Example: If you have 15 samples, use --array=0-14%12"
    # To find M: ls -d ${BASE_DATA_DIR}/*/ | grep -v 'r2c2' | grep -v 'reference' | grep -v 'log' | grep -v 'illumina' | wc -l
    exit 1
fi

if [ "$SLURM_ARRAY_TASK_ID" -ge "$NUM_TOTAL_SAMPLES" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is out of bounds. Max index is $(($NUM_TOTAL_SAMPLES - 1))."
    echo "The array directive in your sbatch script might be too large."
    exit 1
fi

# Get the specific sample directory for this array task
SAMPLE_PATH_FULL=${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}

# The rest of the processing logic is for this single SAMPLE_PATH_FULL
SAMPLE_DIR_NAME=$(basename "${SAMPLE_PATH_FULL}")
echo "--------------------------------------------------"
echo "Processing sample directory: ${SAMPLE_DIR_NAME} (Task ID: $SLURM_ARRAY_TASK_ID)"
echo "--------------------------------------------------"

# Create output directories for the two different analyses
NOVEL_EVL_DIR="${OUTPUT_BASE_DIR}/novel_evl/${SAMPLE_DIR_NAME}"
REF_EVL_DIR="${OUTPUT_BASE_DIR}/ref_evl/${SAMPLE_DIR_NAME}"

mkdir -p "${NOVEL_EVL_DIR}"
mkdir -p "${REF_EVL_DIR}"

# Determine species and set paths
SPECIES=""
if [[ "${SAMPLE_DIR_NAME}" == es_* ]]; then
    SPECIES="mouse"
    GENOME_FASTA="${MOUSE_GENOME}"
    TUSCO_GTF="${MOUSE_TUSCO_GTF}"
    REF_GTF="${MOUSE_REF_GTF}"
    echo "Species detected: Mouse"
elif [[ "${SAMPLE_DIR_NAME}" == wtc11_* ]]; then
    SPECIES="human"
    GENOME_FASTA="${HUMAN_GENOME}"
    TUSCO_GTF="${HUMAN_TUSCO_GTF}"
    REF_GTF="${HUMAN_REF_GTF}"
    echo "Species detected: Human"
else
    echo "Warning: Could not determine species for ${SAMPLE_DIR_NAME}. Skipping."
    exit 1
fi

# --- Ensure reference genome is uncompressed (required by pyfaidx in FLAIR correct/collapse) ---
if [[ "${GENOME_FASTA}" == *.gz ]]; then
    # Decompress into a shared reference folder under the output base to avoid repeated work
    UNZIPPED_GENOME_FASTA="${OUTPUT_BASE_DIR}/reference/$(basename "${GENOME_FASTA}" .gz)"
    if [ ! -f "${UNZIPPED_GENOME_FASTA}" ]; then
        echo "Genome FASTA ${GENOME_FASTA} is gzipped; decompressing to ${UNZIPPED_GENOME_FASTA}"
        mkdir -p "$(dirname "${UNZIPPED_GENOME_FASTA}")"
        gunzip -c "${GENOME_FASTA}" > "${UNZIPPED_GENOME_FASTA}"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to decompress genome FASTA ${GENOME_FASTA}. Exiting."; exit 1;
        fi
    else
        echo "Using existing decompressed genome FASTA: ${UNZIPPED_GENOME_FASTA}"
    fi
    GENOME_FASTA="${UNZIPPED_GENOME_FASTA}"
fi

# Find and concatenate FASTQ files
# Assuming .fastq.gz or .fq.gz, adjust glob pattern if necessary
FASTQ_FILES=($(find "${BASE_DATA_DIR}/${SAMPLE_DIR_NAME}" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \)))

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "Warning: No FASTQ files found for ${SAMPLE_DIR_NAME} in ${BASE_DATA_DIR}/${SAMPLE_DIR_NAME}. Skipping."
    echo "Error encountered (No FASTQ files). Exiting task at line $LINENO."; exit 1;
fi

echo "Found ${#FASTQ_FILES[@]} FASTQ file(s) for ${SAMPLE_DIR_NAME}:"
for fq in "${FASTQ_FILES[@]}"; do
    echo "  - $fq"
done

# Removed concatenation logic. FLAIR will receive multiple FASTQ files directly.
READS_FOR_FLAIR=("${FASTQ_FILES[@]}") # Ensure it's an array for clarity, though direct use is fine

# --- Load conda and set up FLAIR ---
source_conda
activate_env "${FLAIR_ENV}" "flair"

# --- Run FLAIR align once and reuse the results ---
echo "================= Running FLAIR Align (shared step) ================="
# Create shared directory for alignment results
mkdir -p "${OUTPUT_BASE_DIR}/shared"
SHARED_BED_FILE=$(run_flair_align "${SAMPLE_DIR_NAME}" "${OUTPUT_BASE_DIR}")

echo "FLAIR align completed, using BED file: ${SHARED_BED_FILE}"
# Verify the bed file exists and has content before proceeding
if [ ! -f "${SHARED_BED_FILE}" ]; then
    echo "Error: Bed file ${SHARED_BED_FILE} does not exist after align step. Exiting task at line $LINENO."
    exit 1
fi

if [ ! -s "${SHARED_BED_FILE}" ]; then
    echo "Error: Bed file ${SHARED_BED_FILE} exists but is empty. Exiting task at line $LINENO."
    exit 1
fi

# --- Run FLAIR correct and collapse with TUSCO GTF (novel_evl) ---
echo "================= Running FLAIR Correct/Collapse with TUSCO GTF (novel_evl) ================="
if ! NOVEL_EVL_GTF=$(run_flair_correct_collapse "${SAMPLE_DIR_NAME}" "${NOVEL_EVL_DIR}" "${SHARED_BED_FILE}" "${TUSCO_GTF}" "TUSCO" "${SPECIES}"); then
    echo "Error: FLAIR correct/collapse with TUSCO GTF failed. Exiting."
    exit 1
fi

# --- Run FLAIR correct and collapse with REF GTF (ref_evl) ---
echo "================= Running FLAIR Correct/Collapse with REF GTF (ref_evl) ================="
if ! REF_EVL_GTF=$(run_flair_correct_collapse "${SAMPLE_DIR_NAME}" "${REF_EVL_DIR}" "${SHARED_BED_FILE}" "${REF_GTF}" "REF" "${SPECIES}"); then
    echo "Error: FLAIR correct/collapse with REF GTF failed. Exiting."
    exit 1
fi

# --- Run SQANTI3 QC for TUSCO GTF output (novel_evl) ---
echo "================= Running SQANTI3 QC on TUSCO GTF Output (novel_evl) ================="
run_sqanti3_qc "${SAMPLE_DIR_NAME}" "${NOVEL_EVL_DIR}" "${NOVEL_EVL_GTF}" "${REF_GTF}" "TUSCO"

# --- Run SQANTI3 QC for REF GTF output (ref_evl) ---
echo "================= Running SQANTI3 QC on REF GTF Output (ref_evl) ================="
run_sqanti3_qc "${SAMPLE_DIR_NAME}" "${REF_EVL_DIR}" "${REF_EVL_GTF}" "${REF_GTF}" "REF"

echo "--------------------------------------------------"
echo "Task $SLURM_ARRAY_TASK_ID for sample ${SAMPLE_DIR_NAME} completed."
echo "Pipeline task finished at $(date)" 
