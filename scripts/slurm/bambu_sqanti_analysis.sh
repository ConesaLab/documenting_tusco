#!/bin/bash
#SBATCH --job-name=bambu_sqanti_pipeline
#SBATCH --output=bambu_sqanti_pipeline_%A_%a.out
#SBATCH --error=bambu_sqanti_pipeline_%A_%a.err
#SBATCH --time=0-24:00:00  # Adjust as necessary
#SBATCH --mem=64G
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --array=0-9%10   # Update after counting samples
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tianyuan.liu@csic.es

###############################################################################
# Script: bambu_sqanti_analysis.sh
# Description: For each LRGASP sample, build transcriptomes with Bambu using
#              (1) TUSCO annotations and (2) GENCODE reference annotations, then
#              evaluate both assemblies with SQANTI3 against the reference GTF.
#              Inspired by flair_sqanti_analysis.sh and isoquant_sqanti_analysis.sh
###############################################################################

# ----------------------------- CONFIGURATION -------------------------------- #
BAMBU_ENV="bambu_env"                 # conda env that has R + bambu + dependencies
SQANTI3_ENV="SQANTI3.env"         # conda env with SQANTI3 + python deps
SQANTI3_PATH="${SQANTI3_PATH:-}"  # optional: directory containing sqanti3_qc.py (leave empty if in PATH)

BASE_DATA_DIR="/storage/gge/Tian/lrgasp_data"        # FASTQ root
OUTPUT_BASE_DIR="/storage/gge/Tian/lrgasp_analysis/bambu_sq3"  # main out dir

# ---- Reference files ----
HUMAN_GENOME="/storage/gge/Tian/lrgasp_data/reference/GRCh38.primary_assembly.genome.fa"
HUMAN_TUSCO_GTF="/storage/gge/Tian/lrgasp_data/reference/tusco_human_full.gtf"
HUMAN_REF_GTF="/storage/gge/Tian/lrgasp_data/reference/gencode.v47.annotation.gtf"

MOUSE_GENOME="/storage/gge/Tian/lrgasp_data/reference/GRCm39.primary_assembly.genome.fa"
MOUSE_TUSCO_GTF="/storage/gge/Tian/lrgasp_data/reference/tusco_mouse_full.gtf"
MOUSE_REF_GTF="/storage/gge/Tian/lrgasp_data/reference/gencode.vM36.annotation.gtf"

THREADS=${SLURM_CPUS_PER_TASK}

# -------------------------- HELPER FUNCTIONS -------------------------------- #
source_conda() {
    echo "[INFO] Loading Anaconda module..."
    module load anaconda
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Failed to load Anaconda module."; exit 1;
    fi
}

activate_env() {
    local env_name=$1; local exe_to_check=$2
    echo "[INFO] Activating conda env $env_name"
    conda activate "$env_name"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Cannot activate env $env_name"; exit 1;
    fi
    if [[ -n "$exe_to_check" ]]; then
        command -v "$exe_to_check" >/dev/null 2>&1 || { echo "[ERROR] $exe_to_check not found in PATH"; exit 1; }
    fi
}

# ------------- Alignment step (minimap2 -> sorted BAM) ---------------------- #
run_alignment() {
    local sample_dir_name=$1
    local bam_out_dir="$OUTPUT_BASE_DIR/shared"
    mkdir -p "$bam_out_dir"
    local bam_path="$bam_out_dir/${sample_dir_name}.sorted.bam"

    if [[ -f "$bam_path" ]]; then
        echo "[ALIGN] Skipping – $bam_path already exists" >&2
        echo "$bam_path"; return 0
    fi

    echo "[ALIGN] Building BAM for $sample_dir_name" >&2

    # Choose minimap2 preset depending on platform inferred from dir name
    local preset="splice"  # generic splice preset
    if [[ "$sample_dir_name" == *ont* ]]; then
        preset="splice"          # ONT direct RNA / cDNA
    elif [[ "$sample_dir_name" == *pacbio* ]]; then
        preset="splice"          # PacBio CCS/CLR – splice mode still fine
    fi

    minimap2 -t "$THREADS" -ax "$preset" --secondary=no "$GENOME_FASTA" "${READS_FOR_ALIGN[@]}" | \
        samtools sort -@ "$THREADS" -o "$bam_path" -
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] minimap2 or samtools sort failed" >&2; exit 1;
    fi
    samtools index "$bam_path"
    echo "[ALIGN] Created $bam_path" >&2
    echo "$bam_path"
}

# --------------------------- Bambu runner ---------------------------------- #
run_bambu() {
    local sample_dir_name=$1   # wtc11_cdna_ont
    local output_dir=$2        # novel_evl/sample or ref_evl/sample
    local gtf_file=$3          # annotation to use
    local gtf_type=$4          # label: TUSCO or REF
    local bam_file=$5          # aligned BAM

    mkdir -p "$output_dir"
    local bambu_out_dir="$output_dir/bambu_out"
    local extended_gtf="$bambu_out_dir/extended_annotations.gtf"  # produced by writeBambuOutput()

    if [[ -f "$extended_gtf" ]]; then
        echo "[BAMBU] Skipping – $extended_gtf exists (sample $sample_dir_name, $gtf_type)" >&2
        echo "$extended_gtf"; return 0
    fi

    echo "[BAMBU] Running Bambu for $sample_dir_name ($gtf_type)" >&2

    # Build Rscript snippet and execute (redirect stdout to stderr so that only the
    # final echo statement below is captured by the caller)
    Rscript --vanilla - << EOF 1>&2
        suppressPackageStartupMessages({ library(bambu) })
        bam_file <- "$bam_file"
        ann_file <- "$gtf_file"
        genome_fa <- "$GENOME_FASTA"
        out_dir <- "$bambu_out_dir"
        threads <- $THREADS
        message("Bambu parameters:")
        message("  bam    : ", bam_file)
        message("  ann    : ", ann_file)
        message("  genome : ", genome_fa)
        
        # Configure BiocManager repositories properly
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        
        # Set up repositories correctly
        options(repos = BiocManager::repositories())
        
        # Use parameters according to Bambu 3.8.3
        se <- bambu(reads=bam_file,
                   annotations=ann_file,
                   genome=genome_fa,
                   ncore=threads,
                   NDR=0.2,
                   opt.discovery=list(
                     min.txScore.multiExon=0,
                     min.txScore.singleExon=1,
                     min.sampleNumber=1,
                     remove.subsetTx=FALSE,
                     min.readFractionByGene=0
                   ))
        
        dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
        writeBambuOutput(se, path=out_dir)
EOF

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Bambu failed for $sample_dir_name ($gtf_type)" >&2; exit 1;
    fi

    if [[ ! -f "$extended_gtf" ]]; then
        echo "[ERROR] Expected Bambu output $extended_gtf missing." >&2; exit 1;
    fi
    echo "[BAMBU] Completed – output $extended_gtf" >&2
    echo "$extended_gtf"
}

# --------------------------- SQANTI3 runner -------------------------------- #
run_sqanti3_qc() {
    local sample_dir_name=$1
    local output_dir=$2
    local isoforms_gtf=$3
    local ref_gtf=$4

    local sqanti3_output_dir="$output_dir/sqanti3_out"
    local sqanti3_prefix="${sample_dir_name}_sqanti"
    local classification_txt="$sqanti3_output_dir/${sqanti3_prefix}_classification.txt"

    mkdir -p "$sqanti3_output_dir"

    if [[ -f "$classification_txt" ]]; then
        echo "[SQANTI3] Skipping – $classification_txt exists" >&2; return 0; fi

    echo "[SQANTI3] Running on $isoforms_gtf (sample $sample_dir_name)" >&2

    source_conda
    activate_env "$SQANTI3_ENV" ""

    local genome_fa="$GENOME_FASTA"
    if [[ "$genome_fa" == *.gz ]]; then
        local unzipped_genome="$OUTPUT_BASE_DIR/reference/$(basename "$genome_fa" .gz)"
        mkdir -p "$(dirname "$unzipped_genome")"
        if [[ ! -f "$unzipped_genome" ]]; then
            gunzip -c "$genome_fa" > "$unzipped_genome" || { echo "[ERROR] gunzip failed"; exit 1; }
        fi
        genome_fa="$unzipped_genome"
    fi

    local ref_gtf_sqanti3="$ref_gtf"
    local sqanti3_stderr_log="$sqanti3_output_dir/${sqanti3_prefix}_stderr.log"
    {
        "$SQANTI3_PATH/sqanti3_qc.py" \
            "$isoforms_gtf" \
            "$ref_gtf_sqanti3" \
            "$genome_fa" \
            -d "$sqanti3_output_dir" \
            -o "$sqanti3_prefix" \
            --skipORF \
            -t "$THREADS" \
            --report pdf
    } 2> "$sqanti3_stderr_log"

    # Check exit status of sqanti3_qc.py. If failed, print its stderr and exit.
    if [[ $? -ne 0 ]]; then 
        cat "$sqanti3_stderr_log" >&2
        echo "[ERROR] SQANTI3 failed for $sample_dir_name on $isoforms_gtf" >&2; 
        exit 1; 
    fi

    echo "[SQANTI3] Completed – results in $sqanti3_output_dir" >&2
}

# ------------------------------- MAIN -------------------------------------- #

echo "================ Bambu + SQANTI3 Pipeline ================"
echo "Date: $(date) | SLURM Job ID: $SLURM_JOB_ID | Array task: $SLURM_ARRAY_TASK_ID/$SLURM_ARRAY_JOB_COUNT"

mkdir -p "$OUTPUT_BASE_DIR" || { echo "[ERROR] Cannot create output dir $OUTPUT_BASE_DIR"; exit 1; }
cd "$OUTPUT_BASE_DIR" || { echo "[ERROR] cd to $OUTPUT_BASE_DIR failed"; exit 1; }

# Enumerate sample dirs (exclude r2c2, reference, log, illumina)
SAMPLE_DIRS=()
while IFS= read -r line; do
    SAMPLE_DIRS+=("$line")
done < <(ls -d ${BASE_DATA_DIR}/*/ | grep -v 'r2c2' | grep -v 'reference' | grep -v 'log' | grep -v 'illumina')

NUM_SAMPLES=${#SAMPLE_DIRS[@]}

if [[ -z "$SLURM_ARRAY_TASK_ID" ]] || (( SLURM_ARRAY_TASK_ID >= NUM_SAMPLES )); then
    echo "[ERROR] Invalid SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID (max index $((NUM_SAMPLES-1)))"; exit 1; fi

SAMPLE_PATH_FULL="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"
SAMPLE_DIR_NAME=$(basename "$SAMPLE_PATH_FULL")

echo "Processing sample $SAMPLE_DIR_NAME (index $SLURM_ARRAY_TASK_ID)"

# Output directories per sample
NOVEL_EVL_DIR="$OUTPUT_BASE_DIR/novel_evl/$SAMPLE_DIR_NAME"
REF_EVL_DIR="$OUTPUT_BASE_DIR/ref_evl/$SAMPLE_DIR_NAME"
mkdir -p "$NOVEL_EVL_DIR" "$REF_EVL_DIR"

# Detect species and set references
if [[ "$SAMPLE_DIR_NAME" == es_* ]]; then
    SPECIES="mouse"
    GENOME_FASTA="$MOUSE_GENOME"
    TUSCO_GTF="$MOUSE_TUSCO_GTF"
    REF_GTF="$MOUSE_REF_GTF"
elif [[ "$SAMPLE_DIR_NAME" == wtc11_* ]]; then
    SPECIES="human"
    GENOME_FASTA="$HUMAN_GENOME"
    TUSCO_GTF="$HUMAN_TUSCO_GTF"
    REF_GTF="$HUMAN_REF_GTF"
else
    echo "[WARN] Cannot determine species from $SAMPLE_DIR_NAME. Skipping."; exit 0;
fi

echo "[INFO] Species: $SPECIES"

# Gather FASTQ files
FASTQ_FILES=($(find "$SAMPLE_PATH_FULL" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \)))
if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No FASTQ files found in $SAMPLE_PATH_FULL"; exit 1; fi

echo "[INFO] Found ${#FASTQ_FILES[@]} FASTQ files"
for fq in "${FASTQ_FILES[@]}"; do echo "  - $fq"; done

READS_FOR_ALIGN=("${FASTQ_FILES[@]}")

# ------------- Alignment (shared) ----------------------------------------- #
source_conda
activate_env "$BAMBU_ENV" "Rscript"  # Rscript used for bambu

BAM_FILE=$(run_alignment "$SAMPLE_DIR_NAME")

# ------------- Bambu analyses --------------------------------------------- #

echo "================ Bambu with TUSCO (novel_evl) ================"
NOVEL_EVL_GTF=$(run_bambu "$SAMPLE_DIR_NAME" "$NOVEL_EVL_DIR" "$TUSCO_GTF" "TUSCO" "$BAM_FILE")

echo "================ Bambu with REF GTF (ref_evl) ================"
REF_EVL_GTF=$(run_bambu "$SAMPLE_DIR_NAME" "$REF_EVL_DIR" "$REF_GTF" "REF" "$BAM_FILE")

# ------------- SQANTI3 QC -------------------------------------------------- #

echo "================ SQANTI3 QC on TUSCO assembly =================="
run_sqanti3_qc "$SAMPLE_DIR_NAME" "$NOVEL_EVL_DIR" "$NOVEL_EVL_GTF" "$REF_GTF"

echo "================ SQANTI3 QC on REF assembly ===================="
run_sqanti3_qc "$SAMPLE_DIR_NAME" "$REF_EVL_DIR" "$REF_EVL_GTF" "$REF_GTF"

echo "[DONE] Sample $SAMPLE_DIR_NAME finished at $(date)" 
