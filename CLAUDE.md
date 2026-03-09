# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains code for generating figures in the TUSCO paper, including the `tusco-selector` and `tusco-novel-simulator` Python modules. TUSCO is a gene selection pipeline for transcriptome analysis.

## Common Commands

### Run All Figure Scripts
```bash
./scripts/run_all_figs.sh
```
This downloads the dataset automatically if missing (from S3) and executes all R figure scripts. Stderr logs go to `.err` files next to each script.

### Run a Single Figure Script
```bash
cd figs/figure-03/code && Rscript figure3a-human.R
```
Scripts must be run from their own `code/` directory for relative paths to resolve correctly.

### Run Python Modules
```bash
export PYTHONPATH=src
python -m tusco_selector hsa -o output_gtfs
python -m tusco_novel_simulator --gtf input.gtf --target_genes genes.txt --genome ref.fa --output_fake fake.gtf --output_full full.gtf --log log.txt
```
Species codes: `hsa` (human), `mmu` (mouse), `dre` (zebrafish).

### Build LaTeX Manuscript
From `manuscript/paper/src`:
```bash
TEXINPUTS=../styles//: BIBINPUTS=../assets/reference//: BSTINPUTS=../styles//: latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error -outdir=../build tusco_paper.tex
```

### Setup Conda Environment
```bash
conda env create -f envs/tusco_selector.yml
conda activate tusco_selector-env
```

## Architecture

### Directory Structure
- `src/tusco_selector/` - Python pipeline for TUSCO gene selection (7-step pipeline with expression filtering, splice junction analysis)
- `src/tusco_novel_simulator/` - Generates synthetic GTF annotations for target genes
- `figs/figure-0N/` and `figs/supp-fig-0N/` - Each contains `code/`, `plots/`, `tables/`
- `data/raw/` and `data/processed/` - Input data and derived outputs (download via run script)
- `manuscript/paper/` - LaTeX sources, styles, and assets for the publication

### R Figure Framework
All figure scripts use `scripts/figure_utils.R` which provides:
- `figure_context()` - Creates unified context with path resolution, output directories, and logging
- `ctx$resolve_input(...)` - Searches figure-local `data/`, shared `figs/data/`, and repository `data/` trees
- `ctx$save_plot()` and `ctx$write_table()` - Standardized output to `plots/` and `tables/`
- `theme_tusco()` - Consistent ggplot2 styling (7pt Helvetica base)
- `TUSCO_COLORS` - Standard color palette for human/mouse TUSCO, GENCODE, MANE, SIRVs, ERCCs

### Python Pipeline (tusco_selector)
Entry point: `src/tusco_selector/cli.py:main()`

The CLI accepts filtering thresholds for:
- Splice junction analysis (`--novel-threshold`, `--min-novel-length`, `--tss-scope`)
- Expression filtering (`--expression-source`, GTEx/Bgee/ENCODE parameters)
- Universal vs tissue-specific gene classification

### Data Dependencies
Download the dataset archive and extract to repository root:
```
https://tusco-paper-data.s3.eu-north-1.amazonaws.com/data.zip
```
The `run_all_figs.sh` script handles this automatically.

### Reviewer Response Guidelines

#### Manuscript Formatting Requirements
1. All changes must be highlighted in **red font** in the main manuscript
2. Line numbers must be added to the manuscript
3. In all reviewer responses, indicate:
   - The changed text (if applicable)
   - The specific line number(s) where changes appear

#### Quality Control Before Submission
1. Cross-reference every supplementary figure citation in manuscript
2. Verify sequential ordering of supplementary materials
3. Ensure line numbers are accurate in all reviewer response references
4. Final read-through for oversight errors

#### Key File Locations

| Component | Path |
|-----------|------|
| Main paper | `manuscript/paper/src/tusco_paper.tex` |
| Supplementary figures | `manuscript/paper/src/tusco_paper_supplement_figures.tex` |
| Supplementary tables | `manuscript/paper/src/tusco_paper_supplement_tables.tex` |
| Reviewer response | `reviewer_response/reviewer_response.tex` |
| Figure source code | `figs/figure-0N/code/` and `figs/supp-fig-0N/code/` |
| Figure outputs (plots) | `figs/*/plots/` |
| Manuscript assets | `manuscript/paper/assets/` |

#### Revision Workflow

When making manuscript edits:
1. Update the main paper or supplementary `.tex` files as needed
2. Copy updated figures from `figs/*/plots/` to `manuscript/paper/assets/`
3. Update the reviewer response to reflect changes (include line numbers)
4. Compile all modified LaTeX files to verify changes