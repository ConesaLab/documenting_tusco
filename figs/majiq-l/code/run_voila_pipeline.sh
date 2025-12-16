#!/bin/bash
###############################################################################
# run_voila_pipeline.sh
#
# Purpose: Run the complete MAJIQ-L pipeline for junction comparison
#          1. Extract isoform counts from SQANTI3 classification files
#          2. Run voila lr to map long-reads to short-read splicegraph
#          3. Run MAJIQ-L comparison (majiql) using zarr splicegraph
#
# Prerequisites:
#   - MAJIQ/VOILA installed in ../../../env/
#   - Splicegraph built: data/majiq-l/majiq_out/splicegraph.zarr
#   - SQANTI3 classification files in data/raw/lrgasp/human/
#
# Run from: figs/majiq-l/code/
###############################################################################

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== MAJIQ-L Pipeline ===${NC}"

# Define paths (relative to figs/majiq-l/code/)
BASE_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
LICENSE_FILE="${BASE_DIR}/licenses/majiq/majiq_license_academic_official.lic"
SPLICEGRAPH="${BASE_DIR}/data/majiq-l/majiq_out/splicegraph.zarr"
LRGASP_DIR="${BASE_DIR}/data/raw/lrgasp/human"
OUTPUT_DIR="${BASE_DIR}/figs/majiq-l/tables/raw"
VOILA_OUTPUT_DIR="${BASE_DIR}/figs/majiq-l/tables/voila_outputs"
ENV_DIR="${BASE_DIR}/env"

# Pipelines to process
PIPELINES=(
    "WTC11_drna_ont"
    "WTC11_drna_ont_ls"
    "WTC11_cdna_ont"
    "WTC11_cdna_ont_ls"
    "WTC11_cdna_pacbio"
    "WTC11_cdna_pacbio_ls"
)

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$VOILA_OUTPUT_DIR"

# Check prerequisites
echo -e "\n${YELLOW}Checking prerequisites...${NC}"
if [ ! -f "$LICENSE_FILE" ]; then
    echo -e "${RED}ERROR: License file not found: $LICENSE_FILE${NC}"
    exit 1
fi
if [ ! -d "$SPLICEGRAPH" ]; then
    echo -e "${RED}ERROR: Splicegraph not found: $SPLICEGRAPH${NC}"
    exit 1
fi
echo -e "${GREEN}Prerequisites OK${NC}"

###############################################################################
# Step 1: Extract isoform counts from SQANTI3 classification
###############################################################################
echo -e "\n${YELLOW}Step 1: Extracting isoform counts from SQANTI3...${NC}"

extract_counts() {
    local PIPELINE=$1
    local SQANTI_FILE="${LRGASP_DIR}/${PIPELINE}/${PIPELINE}_classification.txt"
    local OUTPUT_FILE="${VOILA_OUTPUT_DIR}/${PIPELINE}_counts.tsv"

    if [ ! -f "$SQANTI_FILE" ]; then
        echo -e "${RED}  WARNING: SQANTI3 file not found for $PIPELINE${NC}"
        return 1
    fi

    # Extract isoform ID and FL count (columns 1 and 4)
    # Skip header, filter for FL > 0
    awk -F'\t' 'NR>1 && $4>0 {print $1"\t"$4}' "$SQANTI_FILE" > "$OUTPUT_FILE"
    echo "  $PIPELINE: $(wc -l < "$OUTPUT_FILE") isoforms"
}

for PIPELINE in "${PIPELINES[@]}"; do
    extract_counts "$PIPELINE" || true
done

###############################################################################
# Step 2: Run voila lr
###############################################################################
echo -e "\n${YELLOW}Step 2: Running voila lr...${NC}"

run_voila_lr() {
    local PIPELINE=$1
    local GTF_FILE="${LRGASP_DIR}/${PIPELINE}/${PIPELINE}_corrected_fixed.gtf"
    local TSV_FILE="${VOILA_OUTPUT_DIR}/${PIPELINE}_counts.tsv"
    local OUTPUT_FILE="${VOILA_OUTPUT_DIR}/${PIPELINE}.lr.voila"

    if [ -f "$OUTPUT_FILE" ]; then
        echo -e "  ${YELLOW}$PIPELINE: SKIP (already exists)${NC}"
        return 0
    fi

    if [ ! -f "$GTF_FILE" ]; then
        echo -e "${RED}  ERROR: GTF not found: $GTF_FILE${NC}"
        return 1
    fi

    if [ ! -f "$TSV_FILE" ]; then
        echo -e "${RED}  ERROR: Counts file not found: $TSV_FILE${NC}"
        return 1
    fi

    echo "  Processing $PIPELINE..."
    /bin/bash -c "source ${ENV_DIR}/bin/activate && voila --license '$LICENSE_FILE' lr \
        --lr-gtf-file '$GTF_FILE' \
        --lr-tsv-file '$TSV_FILE' \
        -sg '$SPLICEGRAPH' \
        -o '$OUTPUT_FILE' \
        -j 4 --silent" 2>&1 | grep -v "^$" || true

    if [ -f "$OUTPUT_FILE" ]; then
        echo -e "  ${GREEN}$PIPELINE: SUCCESS${NC}"
    else
        echo -e "  ${RED}$PIPELINE: FAILED${NC}"
        return 1
    fi
}

for PIPELINE in "${PIPELINES[@]}"; do
    run_voila_lr "$PIPELINE" || true
done

###############################################################################
# Step 3: Run MAJIQ-L comparison (majiql)
###############################################################################
echo -e "\n${YELLOW}Step 3: Running MAJIQ-L comparison...${NC}"

# Inline Python script for zarr-compatible majiql
run_majiql() {
    local PIPELINE=$1
    local VOILA_FILE="${VOILA_OUTPUT_DIR}/${PIPELINE}.lr.voila"
    local OUTPUT_FILE="${OUTPUT_DIR}/${PIPELINE}_comparison.tsv"

    if [ -f "$OUTPUT_FILE" ]; then
        echo -e "  ${YELLOW}$PIPELINE: SKIP (already exists)${NC}"
        return 0
    fi

    if [ ! -f "$VOILA_FILE" ]; then
        echo -e "${RED}  ERROR: Voila file not found: $VOILA_FILE${NC}"
        return 1
    fi

    echo "  Processing $PIPELINE..."

    /bin/bash -c "source ${ENV_DIR}/bin/activate && python3.12 << 'PYTHON_EOF'
import pickle
import zarr
import numpy as np
import pandas as pd
from csv import DictWriter
from collections import defaultdict

def strip_gene_version(gene_id):
    '''Remove version number from gene ID (e.g., ENSG00000123456.10 -> ENSG00000123456)'''
    if '.' in gene_id and gene_id.split('.')[-1].isdigit():
        return gene_id.rsplit('.', 1)[0]
    return gene_id

# Load splicegraph
sg = zarr.open('$SPLICEGRAPH', 'r')
gene_ids = sg['genes/gene_id'][:]
junc_starts = sg['junctions/start'][:]
junc_ends = sg['junctions/end'][:]
junc_gene_idx = sg['junctions/gene_idx'][:]
junc_denovo = sg['junctions/denovo'][:]
junc_passed = sg['junctions/passed_build'][:]

# Organize SR junctions by gene
gene_junctions = defaultdict(lambda: {'annot_reads': set(), 'annot_no_reads': set(), 'denovo': set()})
for i in range(len(junc_starts)):
    gid = gene_ids[junc_gene_idx[i]]
    junc = (int(junc_starts[i]), int(junc_ends[i]))
    is_annot = junc_denovo[i] == 0
    has_reads = junc_passed[i] == 1
    if is_annot and has_reads:
        gene_junctions[gid]['annot_reads'].add(junc)
    elif is_annot:
        gene_junctions[gid]['annot_no_reads'].add(junc)
    elif has_reads:
        gene_junctions[gid]['denovo'].add(junc)

# Load LR junctions
with open('$VOILA_FILE', 'rb') as f:
    lr = pickle.load(f)
lr_juncs = {}
lr_juncs_by_clean = {}  # Map clean gene ID -> LR junctions
for gid in lr:
    lr_juncs[gid] = set()
    for tx in lr[gid].get('transcripts', []):
        for j in tx.get('junctions', []):
            lr_juncs[gid].add((int(j[0]), int(j[1])))
    # Also index by clean gene ID for version-tolerant matching
    clean_gid = strip_gene_version(gid)
    lr_juncs_by_clean[clean_gid] = lr_juncs[gid]

print(f'Loaded {len(lr_juncs)} genes from LR, {len(lr_juncs_by_clean)} unique clean IDs')

# Compare
with open('$OUTPUT_FILE', 'w') as fw:
    writer = DictWriter(fw, fieldnames=['gene_id', 'All', 'Both de novo', 'MAJIQ & Annotation',
                                        'MAJIQ de novo', 'LR & Annotaion', 'LR de novo', 'Annotation'], delimiter='\t')
    writer.writeheader()
    for gid in gene_ids:
        counts = {'gene_id': gid, 'All': 0, 'Both de novo': 0, 'MAJIQ & Annotation': 0,
                  'MAJIQ de novo': 0, 'LR & Annotaion': 0, 'LR de novo': 0, 'Annotation': 0}
        # Try exact match first, then version-tolerant match
        if gid in lr_juncs:
            lr_g = set(lr_juncs[gid])
        else:
            clean_gid = strip_gene_version(gid)
            lr_g = set(lr_juncs_by_clean.get(clean_gid, []))
        sr = gene_junctions[gid]
        for j in sr['annot_no_reads']:
            if j in lr_g: counts['LR & Annotaion'] += 1; lr_g.discard(j)
            else: counts['Annotation'] += 1
        for j in sr['annot_reads']:
            if j in lr_g: counts['All'] += 1; lr_g.discard(j)
            else: counts['MAJIQ & Annotation'] += 1
        for j in sr['denovo']:
            if j in lr_g: counts['Both de novo'] += 1; lr_g.discard(j)
            else: counts['MAJIQ de novo'] += 1
        counts['LR de novo'] = len(lr_g)
        writer.writerow(counts)

print('Done: $PIPELINE')
PYTHON_EOF"

    if [ -f "$OUTPUT_FILE" ]; then
        echo -e "  ${GREEN}$PIPELINE: SUCCESS${NC}"
    else
        echo -e "  ${RED}$PIPELINE: FAILED${NC}"
        return 1
    fi
}

for PIPELINE in "${PIPELINES[@]}"; do
    run_majiql "$PIPELINE" || true
done

###############################################################################
# Summary
###############################################################################
echo -e "\n${GREEN}=== Pipeline Complete ===${NC}"
echo -e "Raw outputs: $OUTPUT_DIR"
ls -la "$OUTPUT_DIR"/*.tsv 2>/dev/null || echo "No output files found"

echo -e "\n${YELLOW}Next step: Run 'Rscript run_majiq_analysis.R' to generate figures${NC}"
