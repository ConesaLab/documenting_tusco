#!/bin/bash
# Create a single combined Overleaf-ready zip file

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="$BASE_DIR/overleaf_package/tusco_overleaf"

rewrite_path() {
  local file="$1"
  local from="$2"
  local to="$3"
  sed -i '' "s|$from|$to|g" "$file"
}

python "$BASE_DIR/scripts/publish_assets.py" --target all

# Clean and create combined directory
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/fig" "$OUT_DIR/pdf"

# ===== MAIN PAPER FILES =====
cp "$BASE_DIR/manuscript/paper/src/tusco_paper.tex" "$OUT_DIR/"
cp "$BASE_DIR/manuscript/paper/build/tusco_paper.bbl" "$OUT_DIR/"
cp "$BASE_DIR/manuscript/paper/src/sn-jnl.cls" "$OUT_DIR/"
cp "$BASE_DIR/manuscript/paper/styles/sn-nature.bst" "$OUT_DIR/"
cp "$BASE_DIR/ref-tusco/ref-tusco.bib" "$OUT_DIR/"

# Main figures
for figure in fig-1.pdf fig-2.pdf fig-3.pdf fig-4.pdf fig-5.pdf; do
  cp "$BASE_DIR/manuscript/paper/assets/fig/$figure" "$OUT_DIR/fig/"
done

# ===== SUPPLEMENT FILES =====
cp "$BASE_DIR/manuscript/paper/src/tusco_paper_supplement.tex" "$OUT_DIR/"
cp "$BASE_DIR/manuscript/paper/src/tusco_paper_supplement_figures.tex" "$OUT_DIR/"
cp "$BASE_DIR/manuscript/paper/src/tusco_paper_supplement_tables.tex" "$OUT_DIR/"
cp "$BASE_DIR/manuscript/supplementary-note/src/supplementary_note.tex" "$OUT_DIR/"

# Supplementary figures
for figure in fig-s1.pdf fig-s2.pdf fig-s3.pdf fig-s4.pdf fig-s5.pdf fig-s6.pdf fig-s7.pdf; do
  cp "$BASE_DIR/manuscript/paper/assets/fig/$figure" "$OUT_DIR/fig/"
done

# Supplementary note PDFs
cp "$BASE_DIR/manuscript/supplementary-note/assets/pdf/"*.pdf "$OUT_DIR/pdf/"

# ===== REVIEWER RESPONSE FILES =====
cp "$BASE_DIR/reviewer_response/round_1/src/reviewer_response.tex" "$OUT_DIR/"
cp "$BASE_DIR/reviewer_response/round_1/assets/fig/fig-s4.pdf" "$OUT_DIR/fig/fig-s4-reviewer.pdf"

# ===== ADJUST PATHS FOR FLAT STRUCTURE =====
cd "$OUT_DIR"

# Main paper paths
rewrite_path tusco_paper.tex "../assets/fig/" "fig/"
rewrite_path tusco_paper.tex "../../../ref-tusco/ref-tusco" "ref-tusco"

# Supplement paths
rewrite_path tusco_paper_supplement_figures.tex "../assets/fig/" "fig/"
rewrite_path tusco_paper_supplement.tex "\\input{../../supplementary-note/src/supplementary_note.tex}" "\\input{supplementary_note.tex}"
rewrite_path supplementary_note.tex "\\newcommand{\\SupplementNoteAsset}[1]{../../supplementary-note/assets/pdf/#1}" "\\newcommand{\\SupplementNoteAsset}[1]{pdf/#1}"
rewrite_path supplementary_note.tex "\\newcommand{\\SupplementNoteAsset}[1]{../assets/pdf/#1}" "\\newcommand{\\SupplementNoteAsset}[1]{pdf/#1}"

# Reviewer response paths
rewrite_path reviewer_response.tex "../assets/fig/fig-s4.pdf" "fig/fig-s4-reviewer.pdf"

# Create the combined zip
cd "$BASE_DIR/overleaf_package"
rm -f tusco_complete_overleaf.zip
zip -r tusco_complete_overleaf.zip tusco_overleaf/

echo ""
echo "===== Combined Overleaf package created! ====="
echo ""
ls -lh "$BASE_DIR/overleaf_package/tusco_complete_overleaf.zip"
echo ""
echo "Contents:"
unzip -l "$BASE_DIR/overleaf_package/tusco_complete_overleaf.zip" | head -40
echo ""
echo "To use on Overleaf:"
echo "1. Upload tusco_complete_overleaf.zip"
echo "2. Select main document to compile:"
echo "   - tusco_paper.tex (main paper)"
echo "   - tusco_paper_supplement.tex (supplement)"
echo "   - reviewer_response.tex (reviewer response)"
