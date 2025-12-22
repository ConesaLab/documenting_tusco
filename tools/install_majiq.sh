#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

# Create virtual environment
PYTHON_BIN="${PYTHON_BIN:-python3.12}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  PYTHON_BIN="python3"
fi

"$PYTHON_BIN" -m venv env
source env/bin/activate

# Set environment variables for HTSlib
export HTSLIB_LIBRARY_DIR=/opt/homebrew/Cellar/htslib/1.22.1/lib
export HTSLIB_INCLUDE_DIR=/opt/homebrew/Cellar/htslib/1.22.1/include

# Set environment variables for libomp and zlib (Homebrew keg-only)
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/zlib/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/zlib/include"

# Set license file
export MAJIQ_LICENSE_FILE="$REPO_ROOT/licenses/majiq/majiq_license_academic_official.lic"
if [ ! -f "$MAJIQ_LICENSE_FILE" ]; then
  echo "ERROR: MAJIQ license file not found: $MAJIQ_LICENSE_FILE" >&2
  exit 1
fi

# Install MAJIQ
pip install --upgrade pip
pip install ./tools/majiq/academic
