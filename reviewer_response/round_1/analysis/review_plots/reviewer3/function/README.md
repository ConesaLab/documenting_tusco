This folder contains a small, self-contained functional enrichment analysis of the
TUSCO human/mouse gene sets (GO over-representation; offline once Bioconductor
packages are installed).

Inputs
- `data/processed/tusco/tusco_human.tsv`
- `data/processed/tusco/tusco_mouse.tsv`

Run
- From repo root:
  - `Rscript reviewer_response/round_1/analysis/review_plots/reviewer3/function/tusco_go_enrichment.R`

Outputs (written to this folder)
- `human_tusco_gene_annotations.tsv`, `mouse_tusco_gene_annotations.tsv`
- `human_enrichGO_*.tsv`, `mouse_enrichGO_*.tsv` (BP/CC/MF)
- `human_GO_*_dotplot.pdf`, `mouse_GO_*_dotplot.pdf`
- `*_keyword_summary.tsv` (quick check for “ribosome / metabolic / structural” keywords)
