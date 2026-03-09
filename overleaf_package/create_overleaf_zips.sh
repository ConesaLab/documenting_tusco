#!/bin/bash
# Script to create Overleaf-ready zip files for TUSCO manuscript

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="$BASE_DIR/overleaf_package"

rewrite_path() {
  local file="$1"
  local from="$2"
  local to="$3"
  sed -i '' "s|$from|$to|g" "$file"
}

python "$BASE_DIR/scripts/publish_assets.py" --target all

# Clean up any previous attempts
rm -rf "$OUT_DIR"/{main_paper,supplement,reviewer_response}
mkdir -p "$OUT_DIR/main_paper/fig" "$OUT_DIR/supplement/fig" "$OUT_DIR/supplement/pdf" "$OUT_DIR/reviewer_response/fig"

# ===== MAIN PAPER =====
echo "Creating main paper package..."
# Main tex file
cp "$BASE_DIR/manuscript/paper/src/tusco_paper.tex" "$OUT_DIR/main_paper/"
# BBL file (compiled bibliography)
cp "$BASE_DIR/manuscript/paper/build/tusco_paper.bbl" "$OUT_DIR/main_paper/"
# Class file
cp "$BASE_DIR/manuscript/paper/src/sn-jnl.cls" "$OUT_DIR/main_paper/"
# Style files (BST)
cp "$BASE_DIR/manuscript/paper/styles/sn-nature.bst" "$OUT_DIR/main_paper/"
# Figures
for figure in fig-1.pdf fig-2.pdf fig-3.pdf fig-4.pdf fig-5.pdf; do
  cp "$BASE_DIR/manuscript/paper/assets/fig/$figure" "$OUT_DIR/main_paper/fig/"
done
# Bibliography source
cp "$BASE_DIR/ref-tusco/ref-tusco.bib" "$OUT_DIR/main_paper/"

cd "$OUT_DIR/main_paper"
# Adjust paths in tex file for Overleaf flat structure
rewrite_path tusco_paper.tex "../assets/fig/" "fig/"
rewrite_path tusco_paper.tex "../../../ref-tusco/ref-tusco" "ref-tusco"
cd "$OUT_DIR"
zip -r tusco_main_paper.zip main_paper/

# ===== SUPPLEMENT =====
echo "Creating supplement package..."
# Supplement master file
cp "$BASE_DIR/manuscript/paper/src/tusco_paper_supplement.tex" "$OUT_DIR/supplement/"
# Supplement sub-files
cp "$BASE_DIR/manuscript/paper/src/tusco_paper_supplement_figures.tex" "$OUT_DIR/supplement/"
cp "$BASE_DIR/manuscript/paper/src/tusco_paper_supplement_tables.tex" "$OUT_DIR/supplement/"
# Supplementary note
cp "$BASE_DIR/manuscript/supplementary-note/src/supplementary_note.tex" "$OUT_DIR/supplement/"
# Class file
cp "$BASE_DIR/manuscript/paper/src/sn-jnl.cls" "$OUT_DIR/supplement/"
# Supplementary figures
for figure in fig-s1.pdf fig-s2.pdf fig-s3.pdf fig-s4.pdf fig-s5.pdf fig-s6.pdf fig-s7.pdf; do
  cp "$BASE_DIR/manuscript/paper/assets/fig/$figure" "$OUT_DIR/supplement/fig/"
done
# Supplementary note PDFs
cp "$BASE_DIR/manuscript/supplementary-note/assets/pdf/ENSMUSG00000020671.pdf" "$OUT_DIR/supplement/pdf/"
cp "$BASE_DIR/manuscript/supplementary-note/assets/pdf/ENSMUSG00000031983.pdf" "$OUT_DIR/supplement/pdf/"
cp "$BASE_DIR/manuscript/supplementary-note/assets/pdf/ENSMUSG00000036376.pdf" "$OUT_DIR/supplement/pdf/"
cp "$BASE_DIR/manuscript/supplementary-note/assets/pdf/ENSMUSG00000045210.pdf" "$OUT_DIR/supplement/pdf/"

cd "$OUT_DIR/supplement"
# Adjust paths in tex files for Overleaf flat structure
rewrite_path tusco_paper_supplement_figures.tex "../assets/fig/" "fig/"
rewrite_path tusco_paper_supplement.tex "\\input{../../supplementary-note/src/supplementary_note.tex}" "\\input{supplementary_note.tex}"
# Fix the supplementary note asset paths
rewrite_path supplementary_note.tex "\\newcommand{\\SupplementNoteAsset}[1]{../../supplementary-note/assets/pdf/#1}" "\\newcommand{\\SupplementNoteAsset}[1]{pdf/#1}"
rewrite_path supplementary_note.tex "\\newcommand{\\SupplementNoteAsset}[1]{../assets/pdf/#1}" "\\newcommand{\\SupplementNoteAsset}[1]{pdf/#1}"
cd "$OUT_DIR"
zip -r tusco_supplement.zip supplement/

# ===== REVIEWER RESPONSE =====
echo "Creating reviewer response package..."
cp "$BASE_DIR/reviewer_response/round_1/src/reviewer_response.tex" "$OUT_DIR/reviewer_response/"
cp "$BASE_DIR/reviewer_response/round_1/assets/fig/fig-s4.pdf" "$OUT_DIR/reviewer_response/fig/"

cd "$OUT_DIR/reviewer_response"
rewrite_path reviewer_response.tex "../assets/fig/" "fig/"

cd "$OUT_DIR"
zip -r tusco_reviewer_response.zip reviewer_response/

echo ""
echo "===== Overleaf packages created successfully! ====="
echo "Location: $OUT_DIR"
echo ""
echo "Files created:"
ls -lh "$OUT_DIR"/*.zip
echo ""
echo "To use on Overleaf:"
echo "1. Upload tusco_main_paper.zip for the main paper"
echo "2. Upload tusco_supplement.zip for the supplementary materials"
echo "3. Upload tusco_reviewer_response.zip for the reviewer response"
