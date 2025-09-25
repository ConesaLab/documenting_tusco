#!/usr/bin/env bash
set -euo pipefail

# Rebuild all intersection "*_common.gtf" sets under the intersection folder
# using the new annotation and re-run SQANTI3 in-place, preserving layout.
#
# It parses folder names like K31_K32, B31_B32_B33, etc., and intersects the
# corresponding single-sample isoform GTFs from SINGLE_DIR to produce
#   <folder>/<folder>_common.gtf
# Then runs SQANTI3 with --skipORF --report skip, writing results into the same folder.
#
# Usage:
#   bash rerun_all_intersection_sqanti.sh [--annotation PATH] [--genome PATH] [--select PATH] [--sqanti PATH] [--single PATH] [--intersection PATH] [--python BIN]
#
# Defaults match this repo layout.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

ANNOTATION_DEFAULT="$REPO_ROOT/data/processed/tusco/hsa/tusco_human.gtf"
GENOME_DEFAULT="$REPO_ROOT/data/raw/reference/mouse/mm39_SIRV.fa"
SELECT_DEFAULT="$REPO_ROOT/data/raw/nih/select_transcripts.py"
SQANTI_DEFAULT="/Users/tianyuan/Desktop/GitHub/SQANTI3/sqanti3_qc.py"
SINGLE_DEFAULT="$REPO_ROOT/data/raw/nih/single_sample"
INTERSECTION_DEFAULT="$REPO_ROOT/data/raw/nih/intersection"
PYTHON_BIN_DEFAULT="python3"

ANNOTATION="$ANNOTATION_DEFAULT"
GENOME="$GENOME_DEFAULT"
SELECT="$SELECT_DEFAULT"
SQANTI="$SQANTI_DEFAULT"
SINGLE_DIR="$SINGLE_DEFAULT"
INTERSECTION_DIR="$INTERSECTION_DEFAULT"
PYTHON_BIN="$PYTHON_BIN_DEFAULT"
FORCE=0
NO_SQANTI=0

usage() {
  cat <<USAGE
Rebuild all intersection common sets and re-run SQANTI3.

Options:
  --annotation PATH      Annotation GTF (default: $ANNOTATION_DEFAULT)
  --genome PATH          Genome FASTA (default: $GENOME_DEFAULT)
  --select PATH          select_transcripts.py (default: $SELECT_DEFAULT)
  --sqanti PATH          sqanti3_qc.py (default: $SQANTI_DEFAULT)
  --single PATH          Directory of single-sample isoforms (default: $SINGLE_DEFAULT)
  --intersection PATH    Intersection root directory (default: $INTERSECTION_DEFAULT)
  --python BIN           Python executable (default: $PYTHON_BIN_DEFAULT)
  --force                Force regenerate even if outputs exist
  --no-sqanti            Skip running SQANTI3 (generate *_common.gtf only)
  -h, --help             Show this help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --annotation) ANNOTATION="$2"; shift 2;;
    --genome) GENOME="$2"; shift 2;;
    --select) SELECT="$2"; shift 2;;
    --sqanti) SQANTI="$2"; shift 2;;
    --single) SINGLE_DIR="$2"; shift 2;;
    --intersection) INTERSECTION_DIR="$2"; shift 2;;
    --python) PYTHON_BIN="$2"; shift 2;;
    --force) FORCE=1; shift;;
    --no-sqanti) NO_SQANTI=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1;;
  esac
done

# Validate deps
for f in "$ANNOTATION" "$GENOME"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: file not found: $f" >&2; exit 1
  fi
done
if [[ ! -f "$SELECT" ]]; then
  echo "Error: select_transcripts.py not found: $SELECT" >&2; exit 1
fi
if [[ ! -d "$SINGLE_DIR" ]]; then
  echo "Error: single-sample dir not found: $SINGLE_DIR" >&2; exit 1
fi
if [[ ! -d "$INTERSECTION_DIR" ]]; then
  echo "Error: intersection dir not found: $INTERSECTION_DIR" >&2; exit 1
fi

if [[ ! -x "$SQANTI" ]] && ! command -v "$SQANTI" >/dev/null 2>&1; then
  echo "Error: sqanti3_qc.py not found or not executable: $SQANTI" >&2; exit 1
fi
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Error: Python not found: $PYTHON_BIN" >&2; exit 1
fi

# Helpers
is_up_to_date() {
  local out_file="$1"; shift
  local refs=("$@")
  if [[ ! -f "$out_file" ]]; then return 1; fi
  for r in "${refs[@]}"; do
    if [[ -f "$r" && "$r" -nt "$out_file" ]]; then return 1; fi
  done
  return 0
}

run_sqanti() {
  local input_gtf="$1"; shift
  local out_dir="$1"; shift
  if [[ "$NO_SQANTI" -eq 1 ]]; then
    echo "[SQANTI3] Skipping due to --no-sqanti"
    return
  fi
  mkdir -p "$out_dir"
  local class_file="$out_dir/sqanti3_classification.txt"
  if [[ "$FORCE" -eq 0 ]] && is_up_to_date "$class_file" "$input_gtf" "$ANNOTATION" "$GENOME" "$SQANTI"; then
    echo "[SQANTI3] Reusing up-to-date results in $out_dir"
    return
  fi
  echo "[SQANTI3] Running on ${input_gtf##*/} → $out_dir"
  if [[ -x "$SQANTI" ]]; then
    "$SQANTI" --isoforms "$input_gtf" --refGTF "$ANNOTATION" --refFasta "$GENOME" -o sqanti3 -d "$out_dir" --skipORF --report skip
  else
    "$SQANTI" --isoforms "$input_gtf" --refGTF "$ANNOTATION" --refFasta "$GENOME" -o sqanti3 -d "$out_dir" --skipORF --report skip
  fi
}

# Main loop: visit only directories whose basename matches BK groups like B31, K31_K32, etc.
shopt -s nullglob
for dir in "$INTERSECTION_DIR"/*/; do
  base="$(basename "$dir")"
  if [[ ! "$base" =~ ^[BK][0-9]+(_[BK][0-9]+)*$ ]]; then
    # Skip non-sample folders (e.g., Figure5b_*.tsv/pdf/png live in root, not in subdirs)
    continue
  fi

  echo "\n=== Processing $base ==="

  # Build list of sample isoform paths from name parts
  IFS='_' read -r -a parts <<< "$base"
  inputs=()
  for p in "${parts[@]}"; do
    src="$SINGLE_DIR/$p.isoforms.gtf"
    if [[ ! -f "$src" ]]; then
      echo "Warning: missing input for $p at $src — skipping $base" >&2
      inputs=()
      break
    fi
    inputs+=("$src")
  done
  if (( ${#inputs[@]} == 0 )); then
    continue
  fi

  out_common="$dir/$base"_common.gtf
  if [[ -f "$out_common" ]]; then
    echo "Removing existing: $out_common"
    rm -f "$out_common"
  fi

  # Create new common GTF by intersection across inputs
  if [[ "$FORCE" -eq 1 ]] || ! is_up_to_date "$out_common" "${inputs[@]}" "$SELECT" "$ANNOTATION"; then
    mode="intersect"
    if (( ${#inputs[@]} == 1 )); then
      mode="merge"
    fi
    echo "[$mode] $base ← ${inputs[*]}"
    "$PYTHON_BIN" "$SELECT" -f "${inputs[@]}" --mode "$mode" --annotation "$ANNOTATION" -o "$out_common"
  else
    echo "[intersect] Reusing up-to-date: ${out_common##*/}"
  fi

  # Quick sanity
  if [[ ! -s "$out_common" ]]; then
    echo "Error: produced empty $out_common — skipping SQANTI" >&2
    continue
  fi

  # Run SQANTI3 into the same folder
  run_sqanti "$out_common" "$dir"

done

echo "\nAll done."
