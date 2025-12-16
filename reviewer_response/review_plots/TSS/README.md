# TSS Window Validation Analysis

This directory contains the complete TSS/CAGE window validation analysis for reviewer response question 3.a, demonstrating that the 300bp window threshold is optimal for TSS filtering.

## Overview

The analysis validates the 300bp CAGE window size used in TUSCO's single-exon gene selection by:

1. **Distance distributions (ECDF)**: Comparing within-gene TSS variation (gene TSS to refTSS peaks) vs. between-gene proximity (gene TSS to other GENCODE genes' TSSs)
2. **Peak count analysis**: Showing how many refTSS peaks are captured as window size increases

## Files

### Scripts
- **`calculate_tss_distances_corrected.py`** - Python script that calculates all distance distributions and peak counts
  - Loads step2 single-isoform genes (1,617 human, 911 mouse)
  - Loads full GENCODE annotation (78,691 human genes, 78,275 mouse genes)
  - Loads refTSS v4.1 peaks (241,049 human, 172,324 mouse)
  - Computes:
    - Within-gene: Distance from each step2 gene's TSS to nearest refTSS peak within ±3000bp
    - Between-gene: Distance from each step2 gene's TSS to nearest other GENCODE gene's TSS
    - Peak counts at window sizes: [0, 100, 200, 300, 500, 1000, 1200, 1500, 2000, 3000]

- **`plot_tss_validation_combined.R`** - R script that creates the two-panel validation figure
  - Panel a: ECDF curves showing within-gene vs between-gene distance distributions
  - Panel b: RefTSS peak count vs window size (median with 25-75% percentile ribbons)
  - Uses standard TUSCO styling:
    - theme_classic, 7pt base font, Helvetica
    - Panel labels via cowplot (lowercase, bold, 7pt)
    - Nature double column width (7.09" × 3.55")
  - Both panels show 0-1500bp x-axis range

### Data Files
- **`tss_within_gene_distances_corrected_human.tsv`** (1,161 genes with nearby refTSS)
- **`tss_within_gene_distances_corrected_mouse.tsv`** (611 genes with nearby refTSS)
- **`tss_between_gene_distances_corrected_human.tsv`** (1,617 genes)
- **`tss_between_gene_distances_corrected_mouse.tsv`** (911 genes)
- **`reftss_peak_counts_by_window_human.tsv`** (10 window sizes)
- **`reftss_peak_counts_by_window_mouse.tsv`** (10 window sizes)

### Output
- **`tss_validation_combined.pdf`** - Final two-panel validation figure

## Key Results

At the 300bp threshold:
- **Human**: 96.6% sensitivity (within-gene), 2.1% collision rate (between-gene), median 1.0 peaks/gene
- **Mouse**: 93.8% sensitivity (within-gene), 4.3% collision rate (between-gene), median 1.0 peaks/gene

**Interpretation**:
- High sensitivity (>93%): Most step2 genes have refTSS within 300bp → captures canonical TSSs
- Low collision rate (<5%): Few step2 genes have other genes' TSSs within 300bp → avoids false filtering
- Optimal peak count: ~1 refTSS peak per gene at 300bp → single TSS capture without over-inclusion

## Usage

### Regenerate analysis from scratch:

```bash
cd /Users/tianyuan/Desktop/github_dev/documenting_tusco/reviewer_response/review_plots/TSS

# 1. Calculate distance distributions and peak counts (takes ~2-3 minutes)
python calculate_tss_distances_corrected.py

# 2. Generate validation plot
Rscript plot_tss_validation_combined.R
```

### Dependencies

**Python**:
- gzip, sys, collections, pathlib, typing, numpy

**R**:
- ggplot2, dplyr, scales, cowplot

## Analysis Details

### Within-gene Distance Calculation
```python
for each step2 single-isoform gene:
    gene_tss = gene's TSS position

    # Find refTSS peaks within ±3000bp window
    nearby_peaks = find_reftss_peaks(
        chrom=gene.chrom,
        strand=gene.strand,
        center=gene_tss,
        window=3000
    )

    if nearby_peaks:
        distance = min(abs(gene_tss - peak.center) for peak in nearby_peaks)
    else:
        # No refTSS within 3000bp - exclude from analysis
        skip
```

**Result**: 1,161 human genes and 611 mouse genes have refTSS within 3000bp

### Between-gene Distance Calculation
```python
for each step2 single-isoform gene:
    gene_tss = gene's TSS position

    # Find all OTHER GENCODE genes on same chrom/strand
    other_genes = get_all_gencode_genes(
        chrom=gene.chrom,
        strand=gene.strand,
        exclude=gene.id
    )

    # Distance to nearest other gene's TSS
    distance = min(abs(gene_tss - other.TSS) for other in other_genes)
```

**Result**: Compares step2 genes to ALL 78K GENCODE genes (not just step2)

### Peak Count Analysis
```python
for window_size in [0, 100, 200, 300, 500, 1000, 1200, 1500, 2000, 3000]:
    for each step2 gene:
        count = number of refTSS peaks within ±window_size of gene_tss

    calculate median, Q25, Q75 across all genes
```

**Result**: Shows peak count plateau around 1.0 at 100-500bp, then increases as neighboring genes' peaks are captured

## Color Scheme

Following standard TUSCO palette:
- **Within-gene (TUSCO)**: Human = `#a8d5a0` (light green), Mouse = `#1b9e77` (dark green)
- **Between-gene (GENCODE)**: Human = `#fdbf6f` (orange), Mouse = `#e66101` (dark orange)

## References

- Used in reviewer response document: `reviewer_response/reviewer_response.tex` (lines 569-589)
- RefTSS v4.1: Abugessaisa et al., 2019 (CAGE-supported transcription start sites)
- GENCODE v49 (human) / vM38 (mouse) annotations
