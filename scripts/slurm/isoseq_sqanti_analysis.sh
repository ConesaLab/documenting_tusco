#!/bin/bash
#SBATCH --job-name=tusco_isoseq_pipeline
#SBATCH --output=tusco_isoseq_pipeline_%A_%a.out
#SBATCH --error=tusco_isoseq_pipeline_%A_%a.err
#SBATCH --time=0-24:00:00
#SBATCH --mem=120G
#SBATCH --partition=global
#SBATCH --qos short
#SBATCH --cpus-per-task=16
#SBATCH --array=0-9
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tianyuan.liu@csic.es

# --- Configuration ---
SQANTI3_ENV="SQANTI3.env"
SQANTI3_PATH="/storage/gge/Tian/SQANTI3"
TAMA_PATH="/storage/gge/Tian/tama"
BASE_DATA_DIR="/storage/gge/Tian/lrgasp_data"
MAP_RESULT_DIR="${BASE_DATA_DIR}/map_result"
OUTPUT_BASE_DIR="/storage/gge/Tian/lrgasp_analysis/isoseq_sq3"

# Reference files - Human
HUMAN_GENOME="/storage/gge/Tian/lrgasp_data/reference/lrgasp_grch38_sirvs.fasta.gz"
HUMAN_TUSCO_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_v39_annotation_sirvs.human.tusco_sim.gtf"
HUMAN_REF_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_v39_annotation_sirvs.human.gtf"

# Reference files - Mouse
MOUSE_GENOME="/storage/gge/Tian/lrgasp_data/reference/lrgasp_grcm39_sirvs.fasta.gz"
MOUSE_TUSCO_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_vM28_sirvs.mouse.tusco_sim.gtf"
MOUSE_REF_GTF="/storage/gge/Tian/lrgasp_data/reference/lrgasp_gencode_vM28_sirvs.mouse.gtf"

# --- SQANTI3 auxiliary annotation files ---
HUMAN_COVERAGE_FILES="/storage/gge/Tian/lrgasp_data/reference/wtc11_rep1SJ.out.tab,/storage/gge/Tian/lrgasp_data/reference/wtc11_rep2SJ.out.tab,/storage/gge/Tian/lrgasp_data/reference/wtc11_rep3SJ.out.tab"
MOUSE_COVERAGE_FILES="/storage/gge/Tian/lrgasp_data/reference/es_rep1SJ.out.tab,/storage/gge/Tian/lrgasp_data/reference/es_rep2SJ.out.tab,/storage/gge/Tian/lrgasp_data/reference/es_rep3SJ.out.tab"

HUMAN_POLYA_PEAK="/storage/gge/Tian/lrgasp_data/reference/QuantSeq_WTC11.all_reps.bed"
MOUSE_POLYA_PEAK="/storage/gge/Tian/lrgasp_data/reference/QuantSeq_ES.all_reps.bed"

HUMAN_CAGE_PEAK="/storage/gge/Tian/lrgasp_data/reference/human.refTSS_v3.1.hg38.bed"
MOUSE_CAGE_PEAK="/storage/gge/Tian/lrgasp_data/reference/mouse.refTSS_v3.1.mm39.bed"

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
    local tool_to_check=$2
    local tool_path=$3
    
    echo "Activating conda environment: ${env_name}"
    conda activate "${env_name}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to activate conda environment ${env_name}."
        exit 1
    fi
    echo "Successfully activated ${env_name}"

    if [ -n "$tool_to_check" ]; then
        if [ -n "$tool_path" ]; then
            echo "Verifying existence of '$tool_to_check' at path: $tool_path/$tool_to_check"
            if [ -f "$tool_path/$tool_to_check" ]; then
                echo "'$tool_to_check' found at: $tool_path/$tool_to_check"
                export PATH="$tool_path:$PATH"
                echo "Added $tool_path to PATH"
            else
                echo "Error: $tool_path/$tool_to_check not found"
                exit 1
            fi
        else
            echo "Verifying command '$tool_to_check' in PATH..."
            command -v "$tool_to_check"
            if [ $? -ne 0 ]; then
                echo "Error: command '$tool_to_check' not found in PATH after activating ${env_name}."
                exit 1
            else
                echo "'$tool_to_check' is available at: $(command -v "$tool_to_check")"
            fi
        fi
    fi
}

# --- TAMA Collapse Function ---
run_tama_collapse() {
    local sample_dir_name=$1
    local output_dir=$2
    local mapped_bam=$3
    
    echo "[TAMA COLLAPSE] Starting for ${sample_dir_name}" >&2
    
    mkdir -p "${output_dir}"
    
    local collapse_prefix="${output_dir}/${sample_dir_name}"
    local collapse_bed="${collapse_prefix}.bed"
    local collapse_gtf="${collapse_prefix}.collapsed.gtf"
    
    # Check if final GTF output already exists
    if [ -f "${collapse_gtf}" ]; then
        echo "[TAMA COLLAPSE] Skipping - output ${collapse_gtf} already exists" >&2
        echo "${collapse_gtf}"
        return 0
    fi
    
    # Verify input BAM file exists
    if [ ! -f "${mapped_bam}" ]; then
        echo "Error: Mapped BAM file ${mapped_bam} not found for TAMA collapse! Exiting." >&2
        exit 1
    fi
    
    echo "[TAMA COLLAPSE] Running for ${sample_dir_name} (output doesn't exist yet)" >&2
    echo "Input BAM: ${mapped_bam}" >&2
    echo "Output BED: ${collapse_bed}" >&2
    echo "Final GTF: ${collapse_gtf}" >&2
    
    # Activate tama environment for TAMA collapse
    echo "Activating tama environment for TAMA collapse..."
    conda activate tama
    if [ $? -ne 0 ]; then
        echo "Error: Failed to activate tama conda environment."
        exit 1
    fi
    
    # Verify tama_collapse.py is available
    command -v python
    if [ $? -ne 0 ]; then
        echo "Error: python not found in tama environment"
        exit 1
    fi
    
    # Check if reference genome needs to be unzipped for TAMA
    local ref_genome_tama="${GENOME_FASTA}"
    if [[ "${GENOME_FASTA}" == *.gz ]]; then
        echo "Reference genome for TAMA is gzipped. Checking if unzipped version exists or unzipping..."
        local unzipped_genome_fasta="${OUTPUT_BASE_DIR}/reference/$(basename "${GENOME_FASTA}" .gz)"
        if [ ! -f "${unzipped_genome_fasta}" ]; then
            echo "Unzipping ${GENOME_FASTA} to ${unzipped_genome_fasta}"
            mkdir -p "$(dirname "${unzipped_genome_fasta}")"
            gunzip -c "${GENOME_FASTA}" > "${unzipped_genome_fasta}"
            if [ $? -ne 0 ]; then
                echo "Error: Failed to unzip genome ${GENOME_FASTA}. Exiting."; exit 1;
            fi
        else
            echo "Using existing unzipped genome: ${unzipped_genome_fasta}"
        fi
        ref_genome_tama="${unzipped_genome_fasta}"
    fi
    
    # Run TAMA collapse
    python "${TAMA_PATH}/tama_collapse.py" \
        -s "${mapped_bam}" \
        -f "${ref_genome_tama}" \
        -p "${collapse_prefix}" \
        -x no_cap \
        -b BAM \
        -c 99 \
        -i 85 \
        -a 10 \
        -m 10 \
        -z 10 \
        -d merge_dup
    
    if [ $? -ne 0 ]; then
        echo "Error: TAMA collapse failed for ${sample_dir_name}" >&2
        exit 1
    fi
    
    if [ ! -f "${collapse_bed}" ]; then
        echo "Error: TAMA collapse BED output (${collapse_bed}) not found for ${sample_dir_name}." >&2
        exit 1
    fi
    
    echo "TAMA collapse completed successfully, output: ${collapse_bed}" >&2
    
    # Now convert BED to GTF using bed2gtf
    echo "[BED2GTF] Converting BED to GTF for ${sample_dir_name}" >&2
    
    # Activate SQANTI3 environment for bed2gtf
    echo "Activating SQANTI3 environment for bed2gtf..."
    conda activate "${SQANTI3_ENV}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to activate SQANTI3 conda environment."
        exit 1
    fi
    
    # Verify bed2gtf is available
    command -v bed2gtf
    if [ $? -ne 0 ]; then
        echo "Error: bed2gtf not found in SQANTI3 environment"
        exit 1
    fi
    
    # Check if TAMA read file exists (needed for bed2gtf)
    local tama_read_file="${collapse_prefix}_read.txt"
    if [ ! -f "${tama_read_file}" ]; then
        echo "Error: TAMA read file (${tama_read_file}) not found for bed2gtf conversion." >&2
        exit 1
    fi
    
    # Run bed2gtf conversion
    bed2gtf \
        --bed "${collapse_bed}" \
        --isoforms "${tama_read_file}" \
        --output "${collapse_gtf}"
    
    if [ $? -ne 0 ]; then
        echo "Error: bed2gtf conversion failed for ${sample_dir_name}" >&2
        exit 1
    fi
    
    if [ ! -f "${collapse_gtf}" ]; then
        echo "Error: GTF output (${collapse_gtf}) not found after bed2gtf conversion." >&2
        exit 1
    fi
    
    echo "BED to GTF conversion completed successfully, final output: ${collapse_gtf}" >&2
    echo "${collapse_gtf}"
}

# --- SQANTI3 QC Function ---
run_sqanti3_qc() {
    local sample_dir_name=$1
    local output_dir=$2
    local input_gff=$3
    local ref_gtf=$4
    local analysis_type=$5  # "TUSCO", "REF", or "EVAL"
    
    local sqanti3_output_dir="${output_dir}/sqanti3_out_${analysis_type}"
    local sqanti3_out_prefix="${sample_dir_name}_sqanti_${analysis_type}"
    local sqanti3_result="${sqanti3_output_dir}/${sqanti3_out_prefix}_classification.txt"
    
    # Check if SQANTI3 output already exists
    if [ -f "${sqanti3_result}" ]; then
        echo "[SQANTI3] Skipping - output ${sqanti3_result} already exists for ${analysis_type}"
        return 0
    fi
    
    # Create output directory
    mkdir -p "${sqanti3_output_dir}"

    echo "[SQANTI3] Starting for ${sample_dir_name} with ${analysis_type} analysis"
    source_conda
    activate_env "${SQANTI3_ENV}" "sqanti3_qc.py" "${SQANTI3_PATH}"

    # Check if reference genome needs to be unzipped
    local ref_genome_sqanti3="${GENOME_FASTA}"
    if [[ "${GENOME_FASTA}" == *.gz ]]; then
        echo "Reference genome for SQANTI3 is gzipped. Checking if unzipped version exists or unzipping..."
        local unzipped_genome_fasta="${OUTPUT_BASE_DIR}/reference/$(basename "${GENOME_FASTA}" .gz)"
        if [ ! -f "${unzipped_genome_fasta}" ]; then
            echo "Unzipping ${GENOME_FASTA} to ${unzipped_genome_fasta}"
            mkdir -p "$(dirname "${unzipped_genome_fasta}")"
            gunzip -c "${GENOME_FASTA}" > "${unzipped_genome_fasta}"
            if [ $? -ne 0 ]; then
                echo "Error: Failed to unzip genome ${GENOME_FASTA}. Exiting."; exit 1;
            fi
        else
            echo "Using existing unzipped genome: ${unzipped_genome_fasta}"
        fi
        ref_genome_sqanti3="${unzipped_genome_fasta}"
    fi

    # Check if reference GTF needs to be unzipped
    local ref_gtf_sqanti3="${ref_gtf}"
    if [[ "${ref_gtf}" == *.gz ]]; then
        echo "Reference GTF for SQANTI3 is gzipped. Checking if unzipped version exists or unzipping..."
        local unzipped_ref_gtf="${OUTPUT_BASE_DIR}/reference/$(basename "${ref_gtf}" .gz)"
        if [ ! -f "${unzipped_ref_gtf}" ]; then
            echo "Unzipping ${ref_gtf} to ${unzipped_ref_gtf}"
            mkdir -p "$(dirname "${unzipped_ref_gtf}")"
            gunzip -c "${ref_gtf}" > "${unzipped_ref_gtf}"
            if [ $? -ne 0 ]; then
                echo "Error: Failed to unzip reference GTF ${ref_gtf}. Exiting."; exit 1;
            fi
        else
            echo "Using existing unzipped reference GTF: ${unzipped_ref_gtf}"
        fi
        ref_gtf_sqanti3="${unzipped_ref_gtf}"
    fi

    # Convert GFF to GTF if needed (SQANTI3 expects GTF input)
    local input_gtf="${input_gff}"
    if [[ "${input_gff}" == *.gff ]]; then
        input_gtf="${input_gff%.gff}.gtf"
        if [ ! -f "${input_gtf}" ]; then
            echo "Converting GFF to GTF: ${input_gff} -> ${input_gtf}"
            # Simple conversion assuming GFF3 format - you might need a more sophisticated conversion
            awk 'BEGIN{OFS="\t"} !/^#/ {if($3=="exon") print $0}' "${input_gff}" | \
            sed 's/ID=[^;]*;//g; s/Parent=[^;]*;//g' > "${input_gtf}"
            
            if [ ! -s "${input_gtf}" ]; then
                echo "Warning: GTF conversion resulted in empty file, using original GFF"
                input_gtf="${input_gff}"
            fi
        fi
    fi

    # Select species-specific auxiliary files for SQANTI3
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
            echo "Error: Coverage file ${sj} not found for SQANTI3! Exiting." >&2; exit 1;
        fi
    done
    if [ ! -f "${SQANTI3_POLYA}" ]; then
        echo "Error: polyA peak file ${SQANTI3_POLYA} not found for SQANTI3! Exiting." >&2; exit 1;
    fi
    if [ ! -f "${SQANTI3_CAGE}" ]; then
        echo "Error: CAGE peak file ${SQANTI3_CAGE} not found for SQANTI3! Exiting." >&2; exit 1;
    fi

    echo "[SQANTI3 QC] Running for ${sample_dir_name} with ${analysis_type} analysis"
    if [ ! -f "${input_gtf}" ]; then
        echo "Error: Input GTF/GFF ${input_gtf} not found for SQANTI3! Exiting."; exit 1;
    fi
    if [ ! -f "${ref_gtf_sqanti3}" ]; then
        echo "Error: Reference GTF ${ref_gtf_sqanti3} not found for SQANTI3! Exiting."; exit 1;
    fi
    if [ ! -f "${ref_genome_sqanti3}" ]; then
        echo "Error: Reference Genome FASTA ${ref_genome_sqanti3} not found for SQANTI3! Exiting."; exit 1;
    fi

    "${SQANTI3_PATH}/sqanti3_qc.py" \
        "${input_gtf}" \
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
        echo "Error: SQANTI3 QC failed for ${sample_dir_name} with ${analysis_type} analysis."
        exit 1
    else
        echo "[SQANTI3] Completed for ${sample_dir_name} with ${analysis_type} analysis. Results in ${sqanti3_output_dir}"
    fi

    return 0
}

# --- SQANTI3 ML Sets and Filter Function ---
run_sqanti3_ml_filter() {
    local sample_dir_name=$1
    local sqanti3_output_dir=$2
    local analysis_type=$3  # "TUSCO" or "REF"
    
    local sqanti3_out_prefix="${sample_dir_name}_sqanti_${analysis_type}"
    local class_file="${sqanti3_output_dir}/${sqanti3_out_prefix}_classification.txt"
    local tp_file="${sqanti3_output_dir}/TP_transcripts.txt"
    local tn_file="${sqanti3_output_dir}/TN_transcripts.txt"
    local ml_output_prefix="${sqanti3_out_prefix}_mlfiltered"
    local filtered_gtf="${sqanti3_output_dir}/${ml_output_prefix}_corrected.gtf"
    
    # Check if ML filtering already completed
    if [ -f "${filtered_gtf}" ]; then
        echo "[SQANTI3 ML] Skipping - ML filtered output ${filtered_gtf} already exists for ${analysis_type}"
        echo "${filtered_gtf}"
        return 0
    fi
    
    if [ ! -f "${class_file}" ]; then
        echo "Error: SQANTI3 classification file ${class_file} not found for ML filtering!"
        exit 1
    fi
    
    echo "[SQANTI3 ML] Starting ML filtering for ${sample_dir_name} with ${analysis_type} analysis"
    
    # Generate TP/TN sets if they don't exist
    if [ ! -f "${tp_file}" ] || [ ! -f "${tn_file}" ]; then
        echo "[SQANTI3 ML] Generating TP/TN transcript sets"
        bash "/storage/gge/Tian/lrgasp_analysis/sqanti_ml_sets.sh" -d "${sqanti3_output_dir}" -p "${tp_file}" -n "${tn_file}"
        
        if [ $? -ne 0 ]; then
            echo "Error: Failed to generate TP/TN sets for ${sample_dir_name}"
            exit 1
        fi
    else
        echo "[SQANTI3 ML] Using existing TP/TN transcript sets"
    fi
    
    # Verify TP/TN files exist and are not empty
    if [ ! -s "${tp_file}" ]; then
        echo "Error: TP set ${tp_file} not found or empty"
        exit 1
    fi
    if [ ! -s "${tn_file}" ]; then
        echo "Error: TN set ${tn_file} not found or empty"
        exit 1
    fi
    
    echo "TP set: ${tp_file} ($(wc -l < "${tp_file}") transcripts)"
    echo "TN set: ${tn_file} ($(wc -l < "${tn_file}") transcripts)"
    
    # Run SQANTI3 ML filter
    cd "${sqanti3_output_dir}" || exit 1
    
    python3 "${SQANTI3_PATH}/sqanti3_filter.py" ml \
        --sqanti_class "${class_file}" \
        -p "${tp_file}" \
        -n "${tn_file}" \
        -o "${ml_output_prefix}" \
        -d "${sqanti3_output_dir}" \
        -c "${THREADS}" \
        --skip_report \
        --max_class_size 3000
    
    if [ $? -ne 0 ]; then
        echo "Error: SQANTI3 ML filtering failed for ${sample_dir_name} with ${analysis_type} analysis"
        exit 1
    else
        echo "[SQANTI3 ML] Completed ML filtering for ${sample_dir_name} with ${analysis_type} analysis"
    fi
    
    # Return the filtered GTF path
    if [ -f "${filtered_gtf}" ]; then
        echo "${filtered_gtf}"
    else
        echo "Error: ML filtered GTF ${filtered_gtf} not found after filtering"
        exit 1
    fi
}

# --- Main Script ---
echo "Starting TAMA and SQANTI3 pipeline script"
echo "Date: $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"

mkdir -p "${OUTPUT_BASE_DIR}"
cd "${OUTPUT_BASE_DIR}" || { echo "Failed to cd to ${OUTPUT_BASE_DIR}"; exit 1; }

# Get list of sample directories
SAMPLE_DIRS=()
while IFS= read -r line; do
    SAMPLE_DIRS+=("$line")
done < <(ls -d ${MAP_RESULT_DIR}/*/ | sort)

echo "Found ${#SAMPLE_DIRS[@]} total sample directories eligible for processing."
NUM_TOTAL_SAMPLES=${#SAMPLE_DIRS[@]}

# Check if SLURM_ARRAY_TASK_ID is set and within bounds
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID is not set. This script is designed to be run as a job array."
    echo "Please submit with --array=0-$((NUM_TOTAL_SAMPLES-1))"
    exit 1
fi

if [ "$SLURM_ARRAY_TASK_ID" -ge "$NUM_TOTAL_SAMPLES" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is out of bounds. Max index is $(($NUM_TOTAL_SAMPLES - 1))."
    exit 1
fi

# Get the specific sample directory for this array task
SAMPLE_PATH_FULL=${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}
SAMPLE_DIR_NAME=$(basename "${SAMPLE_PATH_FULL}")

echo "--------------------------------------------------"
echo "Processing sample directory: ${SAMPLE_DIR_NAME} (Task ID: $SLURM_ARRAY_TASK_ID)"
echo "--------------------------------------------------"

# Create output directories
BASE_OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE_DIR_NAME}"
TUSCO_BUILD_DIR="${BASE_OUTPUT_DIR}/tusco_build"
REF_BUILD_DIR="${BASE_OUTPUT_DIR}/ref_build"
TUSCO_EVAL_DIR="${BASE_OUTPUT_DIR}/tusco_eval"
REF_EVAL_DIR="${BASE_OUTPUT_DIR}/ref_eval"

mkdir -p "${BASE_OUTPUT_DIR}" "${TUSCO_BUILD_DIR}" "${REF_BUILD_DIR}" "${TUSCO_EVAL_DIR}" "${REF_EVAL_DIR}"

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

# Find the mapped BAM file
MAPPED_BAM="${SAMPLE_PATH_FULL}/${SAMPLE_DIR_NAME}.sorted.bam"
if [ ! -f "${MAPPED_BAM}" ]; then
    echo "Error: Mapped BAM file ${MAPPED_BAM} not found for ${SAMPLE_DIR_NAME}. Skipping."
    exit 1
fi

echo "Found mapped BAM file: ${MAPPED_BAM}"

# Load conda and set up environment
source_conda
activate_env "${SQANTI3_ENV}"

# ================= STEP 0: TAMA Collapse (shared) =================
echo "================= STEP 0: TAMA Collapse ================="

COLLAPSED_GTF=$(run_tama_collapse "${SAMPLE_DIR_NAME}" "${BASE_OUTPUT_DIR}" "${MAPPED_BAM}")

if [ ! -f "${COLLAPSED_GTF}" ]; then
    echo "Error: TAMA collapse failed for ${SAMPLE_DIR_NAME}"
    exit 1
fi

echo "TAMA collapse completed. Output: ${COLLAPSED_GTF}"

# ================= BUILD TRANSCRIPTOME 1: TUSCO-filtered =================
echo "================= BUILD TRANSCRIPTOME 1: TUSCO-filtered ================="

# Step 1: SQANTI3 QC with TUSCO reference
echo "Running SQANTI3 QC with TUSCO reference..."
run_sqanti3_qc "${SAMPLE_DIR_NAME}" "${TUSCO_BUILD_DIR}" "${COLLAPSED_GTF}" "${TUSCO_GTF}" "TUSCO"

# Step 2: SQANTI3 ML filtering based on TUSCO QC
echo "Running SQANTI3 ML filtering based on TUSCO QC..."
TUSCO_FILTERED_GTF=$(run_sqanti3_ml_filter "${SAMPLE_DIR_NAME}" "${TUSCO_BUILD_DIR}/sqanti3_out_TUSCO" "TUSCO")

echo "TUSCO-filtered transcriptome created: ${TUSCO_FILTERED_GTF}"

# ================= BUILD TRANSCRIPTOME 2: Reference-filtered =================
echo "================= BUILD TRANSCRIPTOME 2: Reference-filtered ================="

# Step 1: SQANTI3 QC with reference annotation
echo "Running SQANTI3 QC with reference annotation..."
run_sqanti3_qc "${SAMPLE_DIR_NAME}" "${REF_BUILD_DIR}" "${COLLAPSED_GTF}" "${REF_GTF}" "REF"

# Step 2: SQANTI3 ML filtering based on reference QC
echo "Running SQANTI3 ML filtering based on reference QC..."
REF_FILTERED_GTF=$(run_sqanti3_ml_filter "${SAMPLE_DIR_NAME}" "${REF_BUILD_DIR}/sqanti3_out_REF" "REF")

echo "Reference-filtered transcriptome created: ${REF_FILTERED_GTF}"

# ================= EVALUATE BOTH TRANSCRIPTOMES using Reference Annotation =================
echo "================= EVALUATE TRANSCRIPTOME 1 (TUSCO-filtered) with Reference ================="

# Evaluate TUSCO-filtered transcriptome with reference annotation
run_sqanti3_qc "${SAMPLE_DIR_NAME}_tusco_filtered" "${TUSCO_EVAL_DIR}" "${TUSCO_FILTERED_GTF}" "${REF_GTF}" "EVAL"

echo "================= EVALUATE TRANSCRIPTOME 2 (Reference-filtered) with Reference ================="

# Evaluate Reference-filtered transcriptome with reference annotation
run_sqanti3_qc "${SAMPLE_DIR_NAME}_ref_filtered" "${REF_EVAL_DIR}" "${REF_FILTERED_GTF}" "${REF_GTF}" "EVAL"

echo "--------------------------------------------------"
echo "Task $SLURM_ARRAY_TASK_ID for sample ${SAMPLE_DIR_NAME} completed."
echo "Results:"
echo "  Collapsed GTF: ${COLLAPSED_GTF}"
echo "  TUSCO-filtered transcriptome: ${TUSCO_FILTERED_GTF}"
echo "  Reference-filtered transcriptome: ${REF_FILTERED_GTF}"
echo "  TUSCO-filtered evaluation: ${TUSCO_EVAL_DIR}"
echo "  Reference-filtered evaluation: ${REF_EVAL_DIR}"
echo "Pipeline task finished at $(date)" 