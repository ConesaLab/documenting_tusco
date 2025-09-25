#!/bin/bash

#SBATCH --job-name=sqanti3_flair_mouse_isoseq
#SBATCH --output=sqanti3_flair_mouse_%A_%a.out
#SBATCH --error=sqanti3_flair_mouse_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tianyuan.liu@csic.es
#SBATCH --array=0-17

# ~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Load conda environment
# ~~~~~~~~~~~~~~~~~~~~~~~~
module load anaconda
conda activate SQANTI3.env

# ~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Define references
# ~~~~~~~~~~~~~~~~~~~~~~~~
ANNOTATION="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
GENOME="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"
CAGE_peak="/home/tyuan/BUGSI_PAPER/JoinCall/CAGE_peak/refTSS_v4.1_mouse_coordinate.mm39.bed"

OUTPUT_REF_DIR="/home/tyuan/BUGSI_PAPER/REFERENCE/MOUSE/"
REF_CLASSIF="${OUTPUT_REF_DIR}/$(basename "${ANNOTATION}" .gtf)_classification.txt"

SQANTI_rules="/home/tyuan/GitHub/SQANTI3/utilities/filter/filter_default.json"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.1) [ADDED] Check if reference classification exists
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if [ ! -f "$REF_CLASSIF" ]; then
    echo "ERROR: Reference classification file is missing: $REF_CLASSIF"
    echo "       Please run 'sqanti3_reference_classification.slurm' first!"
    exit 1
fi
echo "[Reference QC] Found reference classification file: $REF_CLASSIF"

# ~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Define coverage file groups
# ~~~~~~~~~~~~~~~~~~~~~~~~
B_COVERAGE=(
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/B31/B31_EKRN230014690-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/B32/B32_EKRN230014691-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/B33/B33_EKRN230014692-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/B34/B34_EKRN230014693-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/B35/B35_EKRN230014694-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/B151/B151_EKRN230014700-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/B152/B152_EKRN230014701-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/B153/B153_EKRN230014702-1A_HFHGVDSX7_L3SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/B154/B154_EKRN230014703-1A_HFHGVDSX7_L3SJ.out.tab"
)

K_COVERAGE=(
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/K31/K31_EKRN230014695-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/K32/K32_EKRN230014696-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/K33/K33_EKRN230014697-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/K34/K34_EKRN230014698-1A_HFHGVDSX7_L4SJ.out.tab"
    "/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/K35/K35_EKRN230014699-1A_HFHGVDSX7_L4SJ.out.tab"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~
# Helper function to pick coverage
# ~~~~~~~~~~~~~~~~~~~~~~~~
get_coverage_files() {
    local gtf_file="$1"
    local coverage=""
    
    # Get the basename without the .gtf extension
    local basename_gtf
    basename_gtf=$(basename "$gtf_file" .gtf)
    
    # Check for the special cases by looking for the substrings in the filename.
    if [[ "$basename_gtf" == *B0K100* ]]; then
        # Whether the filename is B0K100.join.isoforms or B0K100_TAMA, use the full B_COVERAGE.
        coverage=$(IFS=, ; echo "${B_COVERAGE[*]}")
    elif [[ "$basename_gtf" == *B100K0* ]]; then
        # Use the full K_COVERAGE for any file that contains B100K0.
        coverage=$(IFS=, ; echo "${K_COVERAGE[*]}")
    elif [[ "$basename_gtf" =~ ^B[0-9]+ ]]; then
        # For other B* files, extract the sample prefix and choose the matching coverage file.
        sample_prefix="${basename_gtf%%.*}"
        for cov_file in "${B_COVERAGE[@]}"; do
            if [[ "$cov_file" == */"$sample_prefix"_* ]]; then
                coverage="$cov_file"
                break
            fi
        done
        if [[ -z "$coverage" ]]; then
            echo "ERROR: Coverage file for sample prefix $sample_prefix not found in B_COVERAGE."
            exit 1
        fi
    elif [[ "$basename_gtf" =~ ^K[0-9]+ ]]; then
        # For K* files, do the analogous search.
        sample_prefix="${basename_gtf%%.*}"
        for cov_file in "${K_COVERAGE[@]}"; do
            if [[ "$cov_file" == */"$sample_prefix"_* ]]; then
                coverage="$cov_file"
                break
            fi
        done
        if [[ -z "$coverage" ]]; then
            echo "ERROR: Coverage file for sample prefix $sample_prefix not found in K_COVERAGE."
            exit 1
        fi
    else
        echo "ERROR: Unrecognized GTF name: $gtf_file"
        exit 1
    fi
    
    echo "$coverage"
}

# ~~~~~~~~~~~~~~~~~~~~~~~~
# 4) GTF list / array setup
# ~~~~~~~~~~~~~~~~~~~~~~~~
GTF_LIST="/home/tyuan/BUGSI_PAPER/JoinCall/result/my_flair_isoseq_gtfs_ar.txt"
mapfile -t GTF_ARRAY < "$GTF_LIST"

if [ "$SLURM_ARRAY_TASK_ID" -ge "${#GTF_ARRAY[@]}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID is out of range"
    exit 1
fi

CURRENT_GTF=${GTF_ARRAY[$SLURM_ARRAY_TASK_ID]}
if [ -z "$CURRENT_GTF" ]; then
  echo "ERROR: No GTF for index $SLURM_ARRAY_TASK_ID"
  exit 1
fi

SAMPLE_ID=$(basename "$CURRENT_GTF" .gtf)

# ~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Determine coverage
# ~~~~~~~~~~~~~~~~~~~~~~~~
COVERAGE_FILES=$(get_coverage_files "$CURRENT_GTF")
IFS=',' read -ra COVERAGE_ARRAY <<< "$COVERAGE_FILES"
for file in "${COVERAGE_ARRAY[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Coverage file does not exist: $file"
        exit 1
    fi
done

# ~~~~~~~~~~~~~~~~~~~~~~~~
# 6) Output folder
# ~~~~~~~~~~~~~~~~~~~~~~~~
OUTPUT_DIR="/home/tyuan/BUGSI_PAPER/JoinCall/result/isoseq_ar/${SAMPLE_ID}"
mkdir -p "$OUTPUT_DIR"

# ~~~~~~~~~~~~~~~~~~~~~~~~
# 7) Run SQANTI3 QC
# ~~~~~~~~~~~~~~~~~~~~~~~~
echo "[SQANTI3 QC] Starting for sample: $SAMPLE_ID"

/home/tyuan/GitHub/SQANTI3/sqanti3_qc.py \
    "$CURRENT_GTF" \
    "$ANNOTATION" \
    "$GENOME" \
    -o "$SAMPLE_ID" \
    -d "$OUTPUT_DIR" \
    --skipORF \
    -c "$COVERAGE_FILES" \
    --report skip \
    --CAGE_peak "$CAGE_peak"

if [ $? -ne 0 ]; then
    echo "[SQANTI3 QC] ERROR: SQANTI3 QC failed on $SAMPLE_ID"
    exit 1
fi

echo "[SQANTI3 QC] Completed for sample: $SAMPLE_ID"

# ~~~~~~~~~~~~~~~~~~~~~~~~
# 8) Run SQANTI3 Filter
# ~~~~~~~~~~~~~~~~~~~~~~~~
CLASSIFICATION_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_classification.txt"
FILTER_OUTPUT_DIR="${OUTPUT_DIR}/filtered"
mkdir -p "$FILTER_OUTPUT_DIR"

# The corrected GTF from the QC step (SQANTI3 typically outputs *_corrected.gtf):
CORRECTED_GTF="${OUTPUT_DIR}/${SAMPLE_ID}_corrected.gtf"

echo "[SQANTI3 Filter] Starting filter for sample: $SAMPLE_ID"
/home/tyuan/GitHub/SQANTI3/sqanti3_filter.py rules \
    "$CLASSIFICATION_FILE" \
    --gtf "$CORRECTED_GTF" \
    -o "${SAMPLE_ID}_filtered" \
    -d "$FILTER_OUTPUT_DIR" \
    --skip_report \
    -j "$SQANTI_rules"

if [ $? -ne 0 ]; then
    echo "[SQANTI3 Filter] ERROR: SQANTI3 Filter failed on $SAMPLE_ID"
    exit 1
fi

echo "[SQANTI3 Filter] Completed for sample: $SAMPLE_ID"

# ~~~~~~~~~~~~~~~~~~~~~~~~
# 9) Run SQANTI3 Rescue
# ~~~~~~~~~~~~~~~~~~~~~~~~
RESCUE_OUTPUT_DIR="${OUTPUT_DIR}/rescued"
mkdir -p "$RESCUE_OUTPUT_DIR"

# Classification from the rules filter step
RESCUE_CLASSIFICATION="${FILTER_OUTPUT_DIR}/${SAMPLE_ID}_filtered_RulesFilter_result_classification.txt"
# GTF from the rules filter step
RESCUE_GTF="${FILTER_OUTPUT_DIR}/${SAMPLE_ID}_filtered.filtered.gtf"
# The corrected FASTA from the initial SQANTI3 QC
CORRECTED_FASTA="${OUTPUT_DIR}/${SAMPLE_ID}_corrected.fasta"

echo "[SQANTI3 Rescue] Starting rescue for sample: $SAMPLE_ID"

/home/tyuan/GitHub/SQANTI3/sqanti3_rescue.py rules \
    "$RESCUE_CLASSIFICATION" \
    -f "$GENOME" \
    --isoforms "$CORRECTED_FASTA" \
    --gtf "$RESCUE_GTF" \
    -g "$ANNOTATION" \
    -k "$REF_CLASSIF" \
    -j "$SQANTI_rules" \
    -o "${SAMPLE_ID}_rescued" \
    -d "$RESCUE_OUTPUT_DIR" \
    --mode full

if [ $? -ne 0 ]; then
    echo "[SQANTI3 Rescue] ERROR: SQANTI3 Rescue failed on $SAMPLE_ID"
    exit 1
fi

echo "[SQANTI3 Rescue] Completed for sample: $SAMPLE_ID"