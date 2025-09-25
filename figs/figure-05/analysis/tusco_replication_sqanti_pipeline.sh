#!/usr/bin/env bash

# TUSCO replication test — final pipeline (merge and intersect via select_transcripts.py → single SQANTI3 run per experiment)
#
# This script performs:
#   1) Per-replicate merges of per-sample GTFs using select_transcripts.py --mode merge (annotation-aware tie-breaking)
#   2) Intersection across replicates using the same select_transcripts.py --mode intersect (match by splice chain, annotation-aware tie-breaking)
#   3) One SQANTI3 run per experiment on the final intersection GTF (no CAGE, skip ORF)
#
# Expected outputs:
#   - 10M_2reps_kidney/rep1_10M.merge.gtf, rep2_10M.merge.gtf, 10M_2reps_kidney.intersection.gtf, sqanti3_out/
#   - 7M_3reps_mixed/repA_7M.merge.gtf, repB_7M.merge.gtf, repC_7M.merge.gtf, 7M_3reps_mixed.intersection.gtf, sqanti3_out/
#
# Usage example:
#   bash tusco_replication_sqanti_pipeline.sh \
#     --genome /path/to/genome.fa \
#     --annotation /path/to/annotation.gtf
#
# Optional overrides:
#   --ref-base, --single, --sel, --sqanti, --out10, --out7

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
# Defaults per repo layout (new data structure)
REF_BASE_DEFAULT="$REPO_ROOT/data/raw/reference/mouse"
SINGLE_DEFAULT="$REPO_ROOT/data/raw/nih/single_sample"
SEL_DEFAULT="$REPO_ROOT/data/raw/nih/select_transcripts.py"
OUT10_DEFAULT="$REPO_ROOT/data/raw/nih/10M_2reps_kidney"
OUT7_DEFAULT="$REPO_ROOT/data/raw/nih/7M_3reps_mixed"
SQANTI_DEFAULT="/Users/tianyuan/Desktop/GitHub/SQANTI3/sqanti3_qc.py"

# User-configurable vars (can override via flags)
REF_BASE="$REF_BASE_DEFAULT"
SINGLE="$SINGLE_DEFAULT"
SEL="$SEL_DEFAULT"
OUT10="$OUT10_DEFAULT"
OUT7="$OUT7_DEFAULT"
SQANTI="$SQANTI_DEFAULT"
GENOME=""
ANNOTATION=""

PYTHON_BIN="python3"
# Conda activation (best-effort)
CONDA_ENV_DEFAULT="SQANTI3.env"
CONDA_ENV="$CONDA_ENV_DEFAULT"
CONDA_SH=""

usage() {
  cat <<USAGE
TUSCO replication pipeline

Required:
  --genome PATH        Genome FASTA file (e.g., GRCm39.fa)
  --annotation PATH    Annotation GTF file (e.g., gencode.vM37.annotation.gtf)

Optional overrides:
  --ref-base PATH      Reference base dir (default: $REF_BASE_DEFAULT)
  --single PATH        Per-sample GTF dir (default: $SINGLE_DEFAULT)
  --sel PATH           select_transcripts.py for merge/intersect (default: $SEL_DEFAULT)
  --sqanti PATH        sqanti3_qc.py path or command (default: $SQANTI_DEFAULT)
  --out10 PATH         Output root for 10M_2reps_kidney (default: $OUT10_DEFAULT)
  --out7 PATH          Output root for 7M_3reps_mixed (default: $OUT7_DEFAULT)
  --python BIN         Python executable (default: python3)
  --conda-env NAME     Conda env name to activate if needed (default: $CONDA_ENV_DEFAULT)
  --conda-sh PATH      Path to conda.sh (auto-detected if possible)
  --all-10m-combos     Also run all 3 replicate-pair combinations for 10M into subfolders
  --no-sqanti          Skip running SQANTI3 (useful for quick debug)
  --force              Force regeneration even if outputs exist
  -h, --help           Show this help

Notes:
  - Activate your SQANTI3 environment before running:
      conda activate SQANTI3.env
  - This script performs no file splitting and no CAGE integration.
USAGE
}

# Parse args
ALL_10M_COMBOS=0
NO_SQANTI=0
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome) GENOME="$2"; shift 2;;
    --annotation) ANNOTATION="$2"; shift 2;;
    --ref-base) REF_BASE="$2"; shift 2;;
    --single) SINGLE="$2"; shift 2;;
    --sel) SEL="$2"; shift 2;;
    --sqanti) SQANTI="$2"; shift 2;;
    --out10) OUT10="$2"; shift 2;;
    --out7) OUT7="$2"; shift 2;;
    --python) PYTHON_BIN="$2"; shift 2;;
    --conda-env) CONDA_ENV="$2"; shift 2;;
    --conda-sh) CONDA_SH="$2"; shift 2;;
    --all-10m-combos) ALL_10M_COMBOS=1; shift;;
    --no-sqanti) NO_SQANTI=1; shift;;
    --force) FORCE=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown argument: $1"; usage; exit 1;;
  esac
done

# Validate dependencies

if [[ ! -x "$SQANTI" ]] && ! command -v "$SQANTI" >/dev/null 2>&1; then
  echo "Error: sqanti3_qc.py not found or not executable at '$SQANTI'." >&2
  echo "       Set --sqanti to the full path or ensure it is in PATH." >&2
  exit 1
fi

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Error: $PYTHON_BIN not found." >&2
  exit 1
fi

if [[ -z "$GENOME" || -z "$ANNOTATION" ]]; then
  echo "Error: --genome and --annotation are required." >&2
  usage
  exit 1
fi

# Validate inputs
for f in "$GENOME" "$ANNOTATION"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: file not found: $f" >&2
    exit 1
  fi
done

if [[ ! -d "$SINGLE" ]]; then
  echo "Error: single-sample GTF dir not found: $SINGLE" >&2
  exit 1
fi

if [[ ! -f "$SEL" ]]; then
  echo "Error: selection script not found: $SEL" >&2
  exit 1
fi

mkdir -p "$OUT10" "$OUT7"

echo "Using configuration:"
echo "  GENOME      = $GENOME"
echo "  ANNOTATION  = $ANNOTATION"
echo "  SINGLE      = $SINGLE"
echo "  SEL         = $SEL"
echo "  SQANTI      = $SQANTI"
echo "  OUT10       = $OUT10"
echo "  OUT7        = $OUT7"
echo "  CONDA_ENV   = $CONDA_ENV"
echo "  ALL_10M_COMBOS = $ALL_10M_COMBOS"

# Helper: check if output is up-to-date vs inputs (return 0 if up-to-date)
is_up_to_date() {
  local out_file="$1"; shift
  local refs=("$@")
  if [[ ! -f "$out_file" ]]; then
    return 1
  fi
  for r in "${refs[@]}"; do
    if [[ -f "$r" && "$r" -nt "$out_file" ]]; then
      return 1
    fi
  done
  return 0
}

# Helper to run per-replicate merges using select_transcripts.py
run_merge() {
  local out_gtf="$1"; shift
  local inputs=("$@")
  for f in "${inputs[@]}"; do
    if [[ ! -f "$f" ]]; then
      echo "Error: input GTF not found: $f" >&2
      exit 1
    fi
  done
  if [[ "$FORCE" -eq 0 ]] && is_up_to_date "$out_gtf" "${inputs[@]}" "$SEL" "$ANNOTATION"; then
    echo "[merge] Reusing existing up-to-date: ${out_gtf##*/}"
  else
    echo "[merge] Building representative merge: ${out_gtf##*/}"
    "$PYTHON_BIN" "$SEL" -f "${inputs[@]}" --mode merge --annotation "$ANNOTATION" -o "$out_gtf"
  fi
  # Report basic counts for transparency
  for f in "${inputs[@]}"; do
    local cnt
    cnt=$(grep -c $'\ttranscript\t' "$f" || true)
    echo "  input: ${f##*/}: $cnt transcripts"
  done
  local out_cnt
  out_cnt=$(grep -c $'\ttranscript\t' "$out_gtf" || true)
  echo "  merge: ${out_gtf##*/}: $out_cnt transcripts"
}

# Replicate-level merges
rep1_10M="$OUT10/rep1_10M.merge.gtf"
rep2_10M="$OUT10/rep2_10M.merge.gtf"
repA_7M="$OUT7/repA_7M.merge.gtf"
repB_7M="$OUT7/repB_7M.merge.gtf"
repC_7M="$OUT7/repC_7M.merge.gtf"

echo "Step 1/3: Representative merges per simulated replicate"
run_merge "$rep1_10M" "$SINGLE/K31.isoforms.gtf" "$SINGLE/K34.isoforms.gtf"
run_merge "$rep2_10M" "$SINGLE/K32.isoforms.gtf" "$SINGLE/K35.isoforms.gtf"
run_merge "$repA_7M" "$SINGLE/K33.isoforms.gtf" "$SINGLE/B31.isoforms.gtf"
run_merge "$repB_7M" "$SINGLE/B32.isoforms.gtf" "$SINGLE/B33.isoforms.gtf"
run_merge "$repC_7M" "$SINGLE/B34.isoforms.gtf" "$SINGLE/B35.isoforms.gtf"

# Quick sanity checks (non-fatal)
echo "Sanity (transcript lines count):"
grep -H -c $'\ttranscript\t' "$OUT10"/rep*_10M.merge.gtf || true
grep -H -c $'\ttranscript\t' "$OUT7"/rep*_7M.merge.gtf || true

echo "Step 2/3: Intersection across replicates (annotation-aware)"
final10="$OUT10/10M_2reps_kidney.intersection.gtf"
final7="$OUT7/7M_3reps_mixed.intersection.gtf"

if [[ "$FORCE" -eq 0 ]] && is_up_to_date "$final10" "$rep1_10M" "$rep2_10M" "$SEL" "$ANNOTATION"; then
  echo "[intersect] Reusing existing up-to-date: ${final10##*/}"
else
  "$PYTHON_BIN" "$SEL" -f "$rep1_10M" "$rep2_10M" --mode intersect --annotation "$ANNOTATION" -o "$final10"
fi
if [[ "$FORCE" -eq 0 ]] && is_up_to_date "$final7" "$repA_7M" "$repB_7M" "$repC_7M" "$SEL" "$ANNOTATION"; then
  echo "[intersect] Reusing existing up-to-date: ${final7##*/}"
else
  "$PYTHON_BIN" "$SEL" -f "$repA_7M" "$repB_7M" "$repC_7M" --mode intersect --annotation "$ANNOTATION" -o "$final7"
fi

echo "Sanity (intersection transcript counts):"
cnt10=$(grep -c $'\ttranscript\t' "$final10" || true)
cnt7=$(grep -c $'\ttranscript\t' "$final7" || true)
echo "$final10: $cnt10"
echo "$final7: $cnt7"

if [[ "${cnt10:-0}" -eq 0 ]]; then
  echo "Error: intersection for 10M_2reps_kidney is empty (no transcript features)." >&2
  echo "       Confirm inputs and splice-chain matching in: $SEL" >&2
  exit 1
fi
if [[ "${cnt7:-0}" -eq 0 ]]; then
  echo "Error: intersection for 7M_3reps_mixed is empty (no transcript features)." >&2
  echo "       Confirm inputs and splice-chain matching in: $SEL" >&2
  exit 1
fi

#############################################
# If requested, prepare 10M combos first, then run SQANTI on them
#############################################

# SQANTI runner (defined early so it can be reused below)
run_sqanti() {
  local input_gtf="$1"; shift
  local out_dir="$1"; shift
  mkdir -p "$out_dir"
  local class_file="$out_dir/sqanti3_classification.txt"
  if [[ "$FORCE" -eq 0 ]] && is_up_to_date "$class_file" "$input_gtf" "$ANNOTATION" "$GENOME" "$SQANTI"; then
    echo "[SQANTI3] Reusing existing results in $out_dir"
    return
  fi
  echo "[SQANTI3] Running on ${input_gtf##*/} → $out_dir"
  if [[ -x "$SQANTI" ]]; then
    "$SQANTI" --isoforms "$input_gtf" --refGTF "$ANNOTATION" --refFasta "$GENOME" -o sqanti3 -d "$out_dir" --skipORF --report skip
  else
    # Assume callable from PATH
    "$SQANTI" --isoforms "$input_gtf" --refGTF "$ANNOTATION" --refFasta "$GENOME" -o sqanti3 -d "$out_dir" --skipORF --report skip
  fi
}

# Helper: prepare a 10M combo (create unions and intersection only)
prepare_10m_combo() {
  local dest_dir="$1"; shift
  local p1a="$1"; local p1b="$2"; local p2a="$3"; local p2b="$4"
  mkdir -p "$dest_dir"
  local r1="$dest_dir/rep1_10M.merge.gtf"
  local r2="$dest_dir/rep2_10M.merge.gtf"
  echo "[10M combo PREP] ${p1a}+${p1b} vs ${p2a}+${p2b} → $dest_dir"
  run_merge "$r1" "$SINGLE/${p1a}.isoforms.gtf" "$SINGLE/${p1b}.isoforms.gtf"
  run_merge "$r2" "$SINGLE/${p2a}.isoforms.gtf" "$SINGLE/${p2b}.isoforms.gtf"
  local inter="$dest_dir/10M_2reps_kidney.intersection.gtf"
  if [[ "$FORCE" -eq 0 ]] && is_up_to_date "$inter" "$r1" "$r2" "$SEL" "$ANNOTATION"; then
    echo "  [intersect] Reusing existing up-to-date: ${inter##*/}"
  else
    "$PYTHON_BIN" "$SEL" -f "$r1" "$r2" --mode intersect --annotation "$ANNOTATION" -o "$inter"
  fi
  local cnt
  cnt=$(grep -c $'\ttranscript\t' "$inter" || true)
  echo "  intersection transcripts: $cnt"
}

# Helper: run SQANTI on a prepared 10M combo
sqanti_10m_combo() {
  local dest_dir="$1"; shift
  local inter="$dest_dir/10M_2reps_kidney.intersection.gtf"
  if [[ ! -f "$inter" ]]; then
    echo "Warning: missing intersection at $inter; skipping SQANTI" >&2
    return
  fi
  local cnt
  cnt=$(grep -c $'\ttranscript\t' "$inter" || true)
  if [[ "${cnt:-0}" -eq 0 ]]; then
    echo "Warning: empty intersection at $inter; skipping SQANTI" >&2
    return
  fi
  run_sqanti "$inter" "$dest_dir/sqanti3_out"
}

if [[ "$ALL_10M_COMBOS" -eq 1 ]]; then
  echo "Preparing all 10M replicate-pair combinations under: $OUT10"
  prepare_10m_combo "$OUT10/K31_K34K32_K35" K31 K34 K32 K35
  prepare_10m_combo "$OUT10/K31_K32K34_K35" K31 K32 K34 K35
  prepare_10m_combo "$OUT10/K31_K35K32_K34" K31 K35 K32 K34

  if [[ "$NO_SQANTI" -eq 0 ]]; then
    echo "Running SQANTI3 on prepared 10M combinations"
  else
    echo "Skipping SQANTI3 on combinations due to --no-sqanti"
  fi
fi

echo "Step 3/3: SQANTI3 on final intersection per experiment"

# If combos were requested, run SQANTI on them first (unless skipped)
if [[ "$ALL_10M_COMBOS" -eq 1 && "$NO_SQANTI" -eq 0 ]]; then
  sqanti_10m_combo "$OUT10/K31_K34K32_K35"
  sqanti_10m_combo "$OUT10/K31_K32K34_K35"
  sqanti_10m_combo "$OUT10/K31_K35K32_K34"
fi

# Then run SQANTI for the default intersections (unless skipped)
if [[ "$NO_SQANTI" -eq 0 ]]; then
  run_sqanti "$final10" "$OUT10/sqanti3_out"
  run_sqanti "$final7" "$OUT7/sqanti3_out"
else
  echo "Skipping SQANTI3 due to --no-sqanti"
fi

echo "Done. Expected layout:"
cat <<LAYOUT
$OUT10/
  rep1_10M.merge.gtf
  rep2_10M.merge.gtf
  10M_2reps_kidney.intersection.gtf
  sqanti3_out/
    sqanti3_classification.txt
    ...

$OUT7/
  repA_7M.merge.gtf
  repB_7M.merge.gtf
  repC_7M.merge.gtf
  7M_3reps_mixed.intersection.gtf
  sqanti3_out/
    sqanti3_classification.txt
    ...
LAYOUT

echo "Notes: No splitting; no CAGE. If intersections look small, confirm splice-chain comparison in: $SEL"

