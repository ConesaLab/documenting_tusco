#!/bin/bash
#SBATCH --job-name=stringtie_sqanti_pipeline
#SBATCH --output=stringtie_sqanti_pipeline_%A_%a.out
#SBATCH --error=stringtie_sqanti_pipeline_%A_%a.err
#SBATCH --time=0-24:00:00 # Adjusted time for short QOS (24 hours)
#SBATCH --mem=32G      # Adjusted memory per array task, can be modified
#SBATCH --qos=short       # Adjusted qos, can be modified
#SBATCH --cpus-per-task=8 # Adjusted cpus, can be modified
#SBATCH --array=0-9 # Adjust this based on the number of samples
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tianyuan.liu@csic.es # Please replace with your email

# --- Configuration ---
SQANTI3_ENV="SQANTI3.env"
SQANTI3_PATH="/home/tyuan/GitHub/SQANTI3" # Add the specific path to SQANTI3
BASE_MAP_DIR="/storage/gge/Tian/lrgasp_data/map_result" # Directory with BAM files
OUTPUT_BASE_DIR="/storage/gge/Tian/lrgasp_analysis/stringtie_sq3" # Updated output directory

# Reference files - Human
HUMAN_GENOME="/storage/gge/Tian/lrgasp_data/reference/lrgasp_grch38_sirvs.fasta.gz"
HUMAN_TUSCO_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_v39_annotation_sirvs.human.tusco_sim.gtf"
HUMAN_REF_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_v39_annotation_sirvs.human.gtf"

# Reference files - Mouse
MOUSE_GENOME="/storage/gge/Tian/lrgasp_data/reference/lrgasp_grcm39_sirvs.fasta.gz"
MOUSE_TUSCO_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_vM28_sirvs.mouse.tusco_sim.gtf"
MOUSE_REF_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_vM28_sirvs.mouse.gtf"

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
    local tool_to_check=$2 # e.g., "stringtie" or "sqanti3_qc.py"
    local tool_path=$3     # Full path to the tool if provided (for SQANTI3)

    echo "Activating conda environment: ${env_name}"
    conda activate "${env_name}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to activate conda environment ${env_name}."
        exit 1
    fi
    echo "Successfully activated ${env_name}"
    echo "PATH after activating ${env_name}: ${PATH}"
    echo "CONDA_PREFIX: ${CONDA_PREFIX}"

    if [ -n "$tool_to_check" ]; then
        if [ "$tool_to_check" == "sqanti3_qc.py" ] && [ -n "$tool_path" ]; then
            echo "Verifying existence of '$tool_to_check' at path: $tool_path/$tool_to_check"
            if [ -f "$tool_path/$tool_to_check" ]; then
                echo "'$tool_to_check' found at: $tool_path/$tool_to_check"
            else
                echo "Error: $tool_path/$tool_to_check not found"
                exit 1
            fi
        else # For tools like stringtie, expect them to be in PATH after env activation
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

# --- StringTie Assembly Function ---
run_stringtie_assembly() {
    local sample_dir_name=$1
    local bam_file_path=$2
    local base_output_dir_for_sample=$3 # e.g., ${OUTPUT_BASE_DIR}/stringtie_tusco/${SAMPLE_DIR_NAME}
    local assembly_type=$4 # "tusco", "ref", "novel"
    local guide_gtf_file=$5 # Path to the guide GTF, can be empty

    local stringtie_out_dir="${base_output_dir_for_sample}/stringtie_out"
    mkdir -p "${stringtie_out_dir}"
    local output_gtf="${stringtie_out_dir}/${sample_dir_name}_stringtie_${assembly_type}.gtf"

    echo "[STRINGTIE ${assembly_type^^}] Starting for ${sample_dir_name}" >&2

    if [ ! -f "${output_gtf}" ]; then
        echo "Running StringTie (${assembly_type}) for ${sample_dir_name}" >&2
        if [ ! -f "${bam_file_path}" ]; then
            echo "Error: BAM file ${bam_file_path} not found for StringTie! Exiting task at line $LINENO." >&2; exit 1;
        fi
        
        local stringtie_cmd=("stringtie" "${bam_file_path}" "-o" "${output_gtf}" "-p" "${THREADS}" "-L")

        if [ -n "${guide_gtf_file}" ]; then # Check if a guide GTF path was provided
            if [ ! -f "${guide_gtf_file}" ]; then
                 echo "Error: Guide GTF ${guide_gtf_file} not found for StringTie! Exiting task at line $LINENO." >&2; exit 1;
            fi
            stringtie_cmd+=("-G" "${guide_gtf_file}") # Add -G option with the guide GTF
        fi
        
        echo "Executing StringTie command: ${stringtie_cmd[@]}" >&2
        "${stringtie_cmd[@]}"
        
        if [ $? -ne 0 ]; then
            echo "Error: StringTie failed for ${sample_dir_name} (${assembly_type}). Exiting task at line $LINENO." >&2; exit 1;
        fi
        if [ ! -f "${output_gtf}" ]; then
            echo "Error: StringTie output GTF ${output_gtf} not found for ${sample_dir_name} (${assembly_type}). Exiting task at line $LINENO." >&2; exit 1;
        fi
        echo "StringTie (${assembly_type}) completed. Output: ${output_gtf}" >&2
    else
        echo "[STRINGTIE ${assembly_type^^}] Skipping - output ${output_gtf} already exists." >&2
    fi
    echo "${output_gtf}" # Return the assembled GTF path
}

# --- SQANTI3 QC Function ---
run_sqanti3_qc() {
    local sample_dir_name=$1
    local base_output_dir_for_sample=$2 # e.g., ${OUTPUT_BASE_DIR}/stringtie_tusco/${SAMPLE_DIR_NAME}
    local assembled_gtf=$3
    local reference_gtf_for_sqanti3=$4
    local genome_fasta_for_sqanti3=$5
    local assembly_type_label=$6 # "tusco", "ref", "novel" for naming

    local sqanti3_output_dir="${base_output_dir_for_sample}/sqanti3_out"
    local sqanti3_out_prefix="${sample_dir_name}_stringtie_${assembly_type_label}_sqanti"
    local sqanti3_result_classification="${sqanti3_output_dir}/${sqanti3_out_prefix}_classification.txt"
    
    mkdir -p "${sqanti3_output_dir}"
    
    if [ ! -f "${sqanti3_result_classification}" ]; then
        echo "[SQANTI3 ${assembly_type_label^^}] Starting for ${sample_dir_name}"
        # No need to call source_conda and activate_env here if already done for stringtie in main
        # and SQANTI3 is in the same env or stringtie env is a prerequisite.
        # However, SQANTI3 needs its own specific path for the script.
        # We ensure the environment is active and SQANTI3 is findable.
        # activate_env is called in main, check sqanti3_qc.py specifically here
        if ! command -v "${SQANTI3_PATH}/sqanti3_qc.py" &> /dev/null; then
             echo "Error: sqanti3_qc.py not found at ${SQANTI3_PATH}/sqanti3_qc.py. Ensure SQANTI3_PATH is correct and environment is active."
             exit 1
        fi


        # Check if reference genome for SQANTI3 needs to be unzipped
        local current_genome_for_sqanti3="${genome_fasta_for_sqanti3}"
        local unzipped_genome_target_dir="${OUTPUT_BASE_DIR}/reference_unzipped" # Shared unzipped dir
        mkdir -p "${unzipped_genome_target_dir}"

        if [[ "${genome_fasta_for_sqanti3}" == *.gz ]]; then
            echo "Reference genome for SQANTI3 is gzipped. Checking/unzipping..."
            local unzipped_genome_fasta="${unzipped_genome_target_dir}/$(basename "${genome_fasta_for_sqanti3}" .gz)"
            if [ ! -f "${unzipped_genome_fasta}" ]; then
                echo "Unzipping ${genome_fasta_for_sqanti3} to ${unzipped_genome_fasta}"
                gunzip -c "${genome_fasta_for_sqanti3}" > "${unzipped_genome_fasta}"
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to unzip genome ${genome_fasta_for_sqanti3}. Exiting task at line $LINENO."; exit 1;
                fi
            else
                echo "Using existing unzipped genome: ${unzipped_genome_fasta}"
            fi
            current_genome_for_sqanti3="${unzipped_genome_fasta}"
        fi
        
        # Check if reference GTF for SQANTI3 needs to be unzipped
        local current_ref_gtf_for_sqanti3="${reference_gtf_for_sqanti3}"
        if [[ "${reference_gtf_for_sqanti3}" == *.gz ]]; then
            echo "Reference GTF for SQANTI3 is gzipped. Checking/unzipping..."
            local unzipped_ref_gtf="${unzipped_genome_target_dir}/$(basename "${reference_gtf_for_sqanti3}" .gz)" # Use same dir
            if [ ! -f "${unzipped_ref_gtf}" ]; then
                echo "Unzipping ${reference_gtf_for_sqanti3} to ${unzipped_ref_gtf}"
                gunzip -c "${reference_gtf_for_sqanti3}" > "${unzipped_ref_gtf}"
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to unzip reference GTF ${reference_gtf_for_sqanti3}. Exiting task at line $LINENO."; exit 1;
                fi
            else
                echo "Using existing unzipped reference GTF: ${unzipped_ref_gtf}"
            fi
            current_ref_gtf_for_sqanti3="${unzipped_ref_gtf}"
        fi

        echo "[SQANTI3 QC] Running for ${sample_dir_name} (${assembly_type_label})"
        if [ ! -f "${assembled_gtf}" ]; then
            echo "Error: Assembled GTF ${assembled_gtf} not found for SQANTI3! Exiting task at line $LINENO."; exit 1;
        fi
        if [ ! -f "${current_ref_gtf_for_sqanti3}" ]; then
            echo "Error: Reference GTF ${current_ref_gtf_for_sqanti3} not found for SQANTI3! Exiting task at line $LINENO."; exit 1;
        fi
        if [ ! -f "${current_genome_for_sqanti3}" ]; then
            echo "Error: Reference Genome FASTA ${current_genome_for_sqanti3} not found for SQANTI3! Exiting task at line $LINENO."; exit 1;
        fi
        
        "${SQANTI3_PATH}/sqanti3_qc.py" \
            "${assembled_gtf}" \
            "${current_ref_gtf_for_sqanti3}" \
            "${current_genome_for_sqanti3}" \
            -d "${sqanti3_output_dir}" \
            -o "${sqanti3_out_prefix}" \
            --skipORF \
            -t "${THREADS}" \
            --report pdf

        if [ $? -ne 0 ]; then
            echo "Error: SQANTI3 QC failed for ${sample_dir_name} (${assembly_type_label}). Exiting task at line $LINENO." >&2; exit 1;
        else
            echo "[SQANTI3 ${assembly_type_label^^}] Completed for ${sample_dir_name}. Results in ${sqanti3_output_dir}"
        fi
    else
        echo "[SQANTI3 ${assembly_type_label^^}] Skipping - output ${sqanti3_result_classification} already exists."
    fi
    
    return 0
}


# --- Main Script ---
echo "Starting StringTie and SQANTI3 pipeline script"
echo "Date: $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"

mkdir -p "${OUTPUT_BASE_DIR}"
cd "${OUTPUT_BASE_DIR}" || { echo "Failed to cd to ${OUTPUT_BASE_DIR}"; exit 1; }

# Iterate over relevant directories in lrgasp_data/map_result
# These are the directories containing the BAM files
SAMPLE_DIRS_FULL_PATH=()
while IFS= read -r dir_path; do
    # Check if it's a directory, not a file or symlink (though ls -d should handle this)
    if [ -d "$dir_path" ]; then
        SAMPLE_DIRS_FULL_PATH+=("$dir_path")
    fi
done < <(ls -d ${BASE_MAP_DIR}/*/ | grep -v 'reference' | grep -v 'log') 
# Removed r2c2 and illumina as they might not be relevant for direct StringTie long-read assembly
# Adjust grep -v if other patterns need exclusion from map_result

echo "Found ${#SAMPLE_DIRS_FULL_PATH[@]} total sample directories in ${BASE_MAP_DIR} for processing."
NUM_TOTAL_SAMPLES=${#SAMPLE_DIRS_FULL_PATH[@]}

if [ "$NUM_TOTAL_SAMPLES" -eq 0 ]; then
    echo "Error: No sample directories found in ${BASE_MAP_DIR}. Please check the path and grep filters."
    exit 1
fi

echo "This script is task $SLURM_ARRAY_TASK_ID of a $SLURM_JOB_ID (Array Job ID: $SLURM_ARRAY_JOB_ID) array processing $NUM_TOTAL_SAMPLES samples."
echo "The --array directive in your sbatch script should be 0-$(($NUM_TOTAL_SAMPLES - 1))"

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
CURRENT_SAMPLE_PATH_FULL=${SAMPLE_DIRS_FULL_PATH[$SLURM_ARRAY_TASK_ID]}
SAMPLE_DIR_NAME=$(basename "${CURRENT_SAMPLE_PATH_FULL}") # e.g., es_captrap_ont

echo "--------------------------------------------------"
echo "Processing sample: ${SAMPLE_DIR_NAME} (Path: ${CURRENT_SAMPLE_PATH_FULL}, Task ID: $SLURM_ARRAY_TASK_ID)"
echo "--------------------------------------------------"

# Locate the BAM file for this sample
# Assuming BAM file is named ${SAMPLE_DIR_NAME}.sorted.bam inside its directory
BAM_FILE="${CURRENT_SAMPLE_PATH_FULL}/${SAMPLE_DIR_NAME}.sorted.bam"
if [ ! -f "${BAM_FILE}" ]; then
    # Alternative check: sometimes it might just be sorted.bam or similar fixed name if not prefixed
    # For now, stick to the convention seen in `es_captrap_ont`
    echo "Error: BAM file ${BAM_FILE} not found for sample ${SAMPLE_DIR_NAME}. Exiting."
    ls -lhtr "${CURRENT_SAMPLE_PATH_FULL}" # List contents for debugging
    exit 1
fi
echo "Using BAM file: ${BAM_FILE}"


# Determine species and set paths
SPECIES=""
GENOME_FASTA=""
TUSCO_GTF_SPECIES="" # For StringTie with TUSCO
REF_GTF_SPECIES=""   # For StringTie with Ref Annotation AND for SQANTI3 reference

if [[ "${SAMPLE_DIR_NAME}" == es_* ]]; then
    SPECIES="mouse"
    GENOME_FASTA="${MOUSE_GENOME}"
    TUSCO_GTF_SPECIES="${MOUSE_TUSCO_GTF}"
    REF_GTF_SPECIES="${MOUSE_REF_GTF}"
    echo "Species detected: Mouse"
elif [[ "${SAMPLE_DIR_NAME}" == wtc11_* ]]; then
    SPECIES="human"
    GENOME_FASTA="${HUMAN_GENOME}"
    TUSCO_GTF_SPECIES="${HUMAN_TUSCO_GTF}"
    REF_GTF_SPECIES="${HUMAN_REF_GTF}"
    echo "Species detected: Human"
else
    echo "Warning: Could not determine species for ${SAMPLE_DIR_NAME}. Skipping."
    exit 1
fi

# --- Load conda and activate environment (needed for StringTie and SQANTI3) ---
source_conda
# Activate the environment that has StringTie and SQANTI3 (and its dependencies like python)
activate_env "${SQANTI3_ENV}" "stringtie" # Check for stringtie first
activate_env "${SQANTI3_ENV}" "sqanti3_qc.py" "${SQANTI3_PATH}" # Then ensure SQANTI3 script is accessible

# Output directories for this sample for each assembly type
STRINGTIE_TUSCO_SAMPLE_DIR="${OUTPUT_BASE_DIR}/stringtie_tusco/${SAMPLE_DIR_NAME}"
STRINGTIE_REF_SAMPLE_DIR="${OUTPUT_BASE_DIR}/stringtie_ref/${SAMPLE_DIR_NAME}"
STRINGTIE_NOVEL_SAMPLE_DIR="${OUTPUT_BASE_DIR}/stringtie_novel/${SAMPLE_DIR_NAME}"

mkdir -p "${STRINGTIE_TUSCO_SAMPLE_DIR}"
mkdir -p "${STRINGTIE_REF_SAMPLE_DIR}"
mkdir -p "${STRINGTIE_NOVEL_SAMPLE_DIR}"

# --- Run StringTie and SQANTI3 for TUSCO-guided assembly ---
echo "================= Processing: StringTie with TUSCO GTF ================="
ASSEMBLED_GTF_TUSCO=$(run_stringtie_assembly "${SAMPLE_DIR_NAME}" "${BAM_FILE}" "${STRINGTIE_TUSCO_SAMPLE_DIR}" "tusco" "${TUSCO_GTF_SPECIES}")
if [ -f "${ASSEMBLED_GTF_TUSCO}" ]; then
    echo "================= Processing: SQANTI3 for TUSCO-guided assembly ================="
    run_sqanti3_qc "${SAMPLE_DIR_NAME}" "${STRINGTIE_TUSCO_SAMPLE_DIR}" "${ASSEMBLED_GTF_TUSCO}" "${REF_GTF_SPECIES}" "${GENOME_FASTA}" "tusco"
else
    echo "Skipping SQANTI3 for TUSCO as StringTie output GTF was not found: ${ASSEMBLED_GTF_TUSCO}"
fi

# --- Run StringTie and SQANTI3 for Reference-guided assembly ---
echo "================= Processing: StringTie with Reference Annotation GTF ================="
ASSEMBLED_GTF_REF=$(run_stringtie_assembly "${SAMPLE_DIR_NAME}" "${BAM_FILE}" "${STRINGTIE_REF_SAMPLE_DIR}" "ref" "${REF_GTF_SPECIES}")
if [ -f "${ASSEMBLED_GTF_REF}" ]; then
    echo "================= Processing: SQANTI3 for Reference-guided assembly ================="
    run_sqanti3_qc "${SAMPLE_DIR_NAME}" "${STRINGTIE_REF_SAMPLE_DIR}" "${ASSEMBLED_GTF_REF}" "${REF_GTF_SPECIES}" "${GENOME_FASTA}" "ref"
else
    echo "Skipping SQANTI3 for REF as StringTie output GTF was not found: ${ASSEMBLED_GTF_REF}"
fi

# --- Run StringTie and SQANTI3 for De Novo assembly (no guide) ---
echo "================= Processing: StringTie De Novo (no guide) ================="
ASSEMBLED_GTF_NOVEL=$(run_stringtie_assembly "${SAMPLE_DIR_NAME}" "${BAM_FILE}" "${STRINGTIE_NOVEL_SAMPLE_DIR}" "novel" "") # Empty string for no guide
if [ -f "${ASSEMBLED_GTF_NOVEL}" ]; then
    echo "================= Processing: SQANTI3 for De Novo assembly ================="
    run_sqanti3_qc "${SAMPLE_DIR_NAME}" "${STRINGTIE_NOVEL_SAMPLE_DIR}" "${ASSEMBLED_GTF_NOVEL}" "${REF_GTF_SPECIES}" "${GENOME_FASTA}" "novel"
else
    echo "Skipping SQANTI3 for NOVEL as StringTie output GTF was not found: ${ASSEMBLED_GTF_NOVEL}"
fi

echo "--------------------------------------------------"
echo "Task $SLURM_ARRAY_TASK_ID for sample ${SAMPLE_DIR_NAME} completed."
echo "Pipeline task finished at $(date)" 