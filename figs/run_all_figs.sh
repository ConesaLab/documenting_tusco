#!/usr/bin/env bash
set -euo pipefail

# Minimal runner to execute all figure R scripts and report success/failure.
# - Runs each R script from its own code directory so relative paths work.
# - Limits BLAS/OMP threads to avoid OpenMP SHM issues on macOS/sandboxes.
# - Captures stderr to .err next to each script and shows a concise summary.

ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"

# Ensure required data is present; download and extract if missing
DATA_DIR="$ROOT_DIR/data"
DATA_URL="https://tusco-paper-data.s3.eu-north-1.amazonaws.com/data.zip"
DATA_ARCHIVE="$ROOT_DIR/tmp/tusco_data.zip"

if [[ ! -d "$DATA_DIR" ]]; then
  echo "[setup] Data directory not found; downloading from $DATA_URL"
  mkdir -p "$ROOT_DIR/tmp"

  if command -v curl >/dev/null 2>&1; then
    curl --fail --location --show-error --output "$DATA_ARCHIVE" "$DATA_URL"
  elif command -v wget >/dev/null 2>&1; then
    wget -nv -O "$DATA_ARCHIVE" "$DATA_URL"
  else
    echo "Neither curl nor wget is available to download data." >&2
    exit 1
  fi

  if ! command -v unzip >/dev/null 2>&1; then
    echo "The unzip command is required to extract downloaded data." >&2
    exit 1
  fi

  unzip -q -o "$DATA_ARCHIVE" -d "$ROOT_DIR"
  rm -f "$DATA_ARCHIVE"

  if [[ ! -d "$DATA_DIR" ]]; then
    echo "Data directory still missing after extraction; please check archive contents." >&2
    exit 1
  fi
fi

# Harden OpenMP/BLAS environment for stability in constrained environments
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OMP_PROC_BIND=${OMP_PROC_BIND:-FALSE}
export OMP_WAIT_POLICY=${OMP_WAIT_POLICY:-PASSIVE}
export OMP_DYNAMIC=${OMP_DYNAMIC:-FALSE}
export OMP_THREAD_LIMIT=${OMP_THREAD_LIMIT:-1}
export KMP_INIT_AT_FORK=${KMP_INIT_AT_FORK:-0}
export KMP_AFFINITY=${KMP_AFFINITY:-disabled}
export KMP_BLOCKTIME=${KMP_BLOCKTIME:-0}
export KMP_SETTINGS=${KMP_SETTINGS:-0}
export KMP_DUPLICATE_LIB_OK=${KMP_DUPLICATE_LIB_OK:-TRUE}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}
export R_MAX_NUM_DLLS=${R_MAX_NUM_DLLS:-150}

SCRIPTS=(
  #figs/figure-01/code/figure1d.R
  #figs/figure-01/code/figure1c_figure-s1.R
  figs/figure-03/code/figure3a-human.R
  figs/figure-03/code/figure3a-mouse.R
  figs/figure-03/code/figure3b-human.R
  figs/figure-03/code/figure3b-mouse.R
  figs/figure-03/code/figure3c.R
  figs/figure-03/code/figure3d-human.R
  figs/figure-03/code/figure3d-mouse.R
  figs/figure-03/code/table-s1.R
  figs/figure-04/code/fig-4b.R
  figs/figure-04/code/fig-s5.R
  figs/figure-05/code/fig-5b-5c-s6.R
  figs/supp-fig-02/code/fig-s2.R
  figs/supp-fig-03/code/FigureS3_correlation.R
  figs/supp-fig-04/code/fig-s4.R
)

STATUS_KEYS=()
STATUS_VALS=()
DURATIONS_S=()

# Helper to run one R script from its directory with proper logs
run_rscript() {
  local script_rel="$1"
  local script_dir
  script_dir=$(dirname "$script_rel")
  local script_base
  script_base=$(basename "$script_rel")
  local err_log="$ROOT_DIR/${script_rel%.R}.err"

  local start_ts end_ts
  start_ts=$(date +%s)

  # Decide LOCAL_ONLY for fig-5b-5c-s6 if precomputed TSVs exist
  local use_local_only=0
  if [[ "$script_base" == "fig-5b-5c-s6.R" ]]; then
    local plot_dir="$ROOT_DIR/figs/figure-05/plots"
    local p1="$plot_dir/fig-5b_points.tsv"
    local p2="$plot_dir/fig-5b_bars.tsv"
    local p3="$plot_dir/fig-5c_points.tsv"
    if [[ -f "$p1" && -f "$p2" && -f "$p3" ]]; then
      use_local_only=1
    fi
  fi

  # Run inside the script directory so relative paths in R match expectations
  local rc=0
  # If FIGS_CONDA_ENV is set, prefer running via conda to pick up R libs
  local runner=(Rscript)
  if [[ -n "${FIGS_CONDA_ENV:-}" ]]; then
    runner=(conda run -n "$FIGS_CONDA_ENV" Rscript)
  fi
  (
    cd "$script_dir"
    if [[ "$use_local_only" == "1" ]]; then
      LOCAL_ONLY=1 "${runner[@]}" "$script_base" >/dev/null 2>"$err_log"
    else
      "${runner[@]}" "$script_base" >/dev/null 2>"$err_log"
    fi
  ) || rc=$?

  end_ts=$(date +%s)
  DURATIONS_S+=("$((end_ts - start_ts))")
  return $rc
}

for s in "${SCRIPTS[@]}"; do
  if [[ ! -f "$s" ]]; then
    STATUS_KEYS+=("$s"); STATUS_VALS+=("SKIP (missing script)"); DURATIONS_S+=("0")
    echo "[skip] $s (not found)"
    continue
  fi
  echo "[run] $s"
  if run_rscript "$s"; then
    STATUS_KEYS+=("$s"); STATUS_VALS+=("OK")
  else
    STATUS_KEYS+=("$s"); STATUS_VALS+=("FAIL")
  fi
done

echo
echo "=== Run Summary ==="
pad() { printf "%-60s" "$1"; }
for i in "${!STATUS_KEYS[@]}"; do
  s=${STATUS_KEYS[$i]}
  st=${STATUS_VALS[$i]}
  dt=${DURATIONS_S[$i]:-0}
  pad "$s"; printf " %4ss  %s\n" "$dt" "$st"
done

echo
echo "Logs (stderr) for failures are at figs/.../*.err next to scripts."
