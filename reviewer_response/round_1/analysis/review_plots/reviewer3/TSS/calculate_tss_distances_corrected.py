#!/usr/bin/env python3
"""
Calculate CORRECTED TSS distance distributions for ECDF plot validation.

This script uses step2 single-isoform genes and calculates:
1. Within-gene: Distance from each gene's TSS to nearest refTSS peak (within ±3000bp)
2. Between-gene: Distance from each gene's TSS to nearest OTHER gene's TSS
3. Peak count: Number of refTSS peaks at different window sizes

Both curves start from the SAME gene set (step2 single-isoform genes).
"""

import gzip
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]


def parse_gtf_attributes(attr_string: str) -> dict:
    """Parse GTF attribute string into dictionary."""
    attrs = {}
    for item in attr_string.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(" ", 1)
        if len(parts) == 2:
            key = parts[0]
            value = parts[1].strip('"')
            attrs[key] = value
    return attrs


def load_step2_genes(gtf_path: Path, species_name: str) -> Dict[str, Tuple[str, str, int]]:
    """Load step2 single-isoform genes.

    Returns:
        {gene_id: (chrom, strand, tss_pos)}
    """
    print(f"Loading {species_name} step2 single-isoform genes from {gtf_path}...")

    gene_tss = {}

    open_fn = gzip.open if str(gtf_path).endswith('.gz') else open

    with open_fn(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'transcript':
                continue

            chrom = fields[0]
            start = int(fields[3])  # 1-based
            end = int(fields[4])    # 1-based inclusive
            strand = fields[6]
            attrs = parse_gtf_attributes(fields[8])

            gene_id = attrs.get('gene_id', '').split('.')[0]
            if not gene_id:
                continue

            # Determine TSS based on strand
            if strand == '+':
                tss = start
            elif strand == '-':
                tss = end
            else:
                continue

            # Should be only one transcript per gene in step2
            if gene_id in gene_tss:
                print(f"  Warning: Gene {gene_id} has multiple transcripts!")

            gene_tss[gene_id] = (chrom, strand, tss)

    print(f"  Loaded {len(gene_tss):,} single-isoform genes")
    return gene_tss


def load_all_gencode_genes(gtf_path: Path, species_name: str) -> Dict[str, Tuple[str, str, int]]:
    """Load ALL genes from GENCODE annotation.

    For genes with multiple transcripts, uses the TSS from the canonical/longest transcript.

    Returns:
        {gene_id: (chrom, strand, tss_pos)}
    """
    print(f"Loading ALL {species_name} GENCODE genes from {gtf_path}...")

    gene_transcripts = defaultdict(list)

    open_fn = gzip.open if str(gtf_path).endswith('.gz') else open

    with open_fn(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'transcript':
                continue

            chrom = fields[0]
            start = int(fields[3])  # 1-based
            end = int(fields[4])    # 1-based inclusive
            strand = fields[6]
            attrs = parse_gtf_attributes(fields[8])

            gene_id = attrs.get('gene_id', '').split('.')[0]
            if not gene_id:
                continue

            # Determine TSS based on strand
            if strand == '+':
                tss = start
            elif strand == '-':
                tss = end
            else:
                continue

            # Store transcript length for later selection
            transcript_length = end - start + 1
            gene_transcripts[gene_id].append((chrom, strand, tss, transcript_length))

    # For each gene, select TSS from longest transcript
    gene_tss = {}
    for gene_id, transcripts in gene_transcripts.items():
        # Sort by transcript length (descending) and take first
        longest = sorted(transcripts, key=lambda x: x[3], reverse=True)[0]
        gene_tss[gene_id] = (longest[0], longest[1], longest[2])  # (chrom, strand, tss)

    print(f"  Loaded {len(gene_tss):,} GENCODE genes")
    return gene_tss


def load_reftss_peaks(bed_path: Path, species_name: str) -> List[Tuple[str, str, int, str]]:
    """Load refTSS peaks from BED file.

    Returns:
        List of (chrom, strand, center, peak_id) tuples
    """
    print(f"Loading {species_name} refTSS peaks from {bed_path}...")

    peaks = []
    open_fn = gzip.open if str(bed_path).endswith('.gz') else open

    with open_fn(bed_path, 'rt') as f:
        next(f)  # Skip header

        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 6:
                continue

            chrom = fields[0]
            start = int(fields[1])  # 0-based
            end = int(fields[2])    # 0-based exclusive (BED format)
            peak_id = fields[3]
            strand = fields[5]

            # Calculate peak center
            center = (start + end) // 2

            peaks.append((chrom, strand, center, peak_id))

    print(f"  Loaded {len(peaks):,} refTSS peaks")
    return peaks


def find_reftss_in_window(
    peaks: List[Tuple],
    gene_chrom: str,
    gene_strand: str,
    gene_tss: int,
    window_bp: int
) -> List[int]:
    """Find refTSS peak centers within ±window_bp of gene TSS.

    Returns:
        List of peak centers (positions)
    """
    matching_peaks = []

    for chrom, strand, center, peak_id in peaks:
        if chrom != gene_chrom or strand != gene_strand:
            continue

        distance = abs(gene_tss - center)
        if distance <= window_bp:
            matching_peaks.append(center)

    return matching_peaks


def calculate_within_gene_distances(
    gene_tss: Dict[str, Tuple],
    peaks: List[Tuple],
    max_window: int = 3000
) -> List[Tuple[str, int]]:
    """Distance from each gene's TSS to nearest refTSS within ±max_window.

    Returns:
        List of (gene_id, distance) tuples
    """
    print(f"Calculating within-gene distances (max window: {max_window}bp)...")

    distances = []
    genes_with_reftss = 0
    genes_without_reftss = 0

    for gene_id, (chrom, strand, tss) in gene_tss.items():
        nearby_peaks = find_reftss_in_window(peaks, chrom, strand, tss, max_window)

        if nearby_peaks:
            min_distance = min(abs(tss - peak_center) for peak_center in nearby_peaks)
            distances.append((gene_id, min_distance))
            genes_with_reftss += 1
        else:
            genes_without_reftss += 1

    print(f"  Genes with refTSS within {max_window}bp: {genes_with_reftss:,}")
    print(f"  Genes without nearby refTSS: {genes_without_reftss:,}")

    return distances


def calculate_between_gene_distances(
    step2_gene_tss: Dict[str, Tuple],
    all_gencode_tss: Dict[str, Tuple]
) -> List[Tuple[str, int]]:
    """Distance from each step2 gene's TSS to nearest OTHER gene's TSS in full GENCODE.

    Args:
        step2_gene_tss: Step2 single-isoform genes {gene_id: (chrom, strand, tss)}
        all_gencode_tss: ALL GENCODE genes {gene_id: (chrom, strand, tss)}

    Returns:
        List of (gene_id, distance) tuples
    """
    print("Calculating between-gene distances (to all GENCODE genes)...")

    # Index all GENCODE genes by chromosome/strand
    chrom_genes: Dict[Tuple[str, str], List[Tuple[str, int]]] = defaultdict(list)
    for gene_id, (chrom, strand, tss) in all_gencode_tss.items():
        chrom_genes[(chrom, strand)].append((gene_id, tss))

    distances = []

    for gene_id, (chrom, strand, tss) in step2_gene_tss.items():
        # Find all OTHER genes on same chrom/strand (excluding self)
        other_genes = [
            (gid, g_tss) for gid, g_tss in chrom_genes[(chrom, strand)]
            if gid != gene_id
        ]

        if other_genes:
            min_distance = min(abs(tss - g_tss) for _, g_tss in other_genes)
            distances.append((gene_id, min_distance))

    print(f"  Calculated distances for {len(distances):,} step2 genes to {len(all_gencode_tss):,} GENCODE genes")

    return distances


def analyze_peak_count_by_window(
    gene_tss: Dict[str, Tuple],
    peaks: List[Tuple]
) -> List[Dict]:
    """Count refTSS peaks at different window sizes.

    Returns:
        List of dicts with window_size, median, q25, q75, etc.
    """
    print("Analyzing peak counts at different window sizes...")

    window_sizes = [1, 10, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 50000, 100000]
    results = []

    for window in window_sizes:
        peak_counts = []

        for gene_id, (chrom, strand, tss) in gene_tss.items():
            peaks_in_window = find_reftss_in_window(peaks, chrom, strand, tss, window)
            peak_counts.append(len(peaks_in_window))

        peak_counts = np.array(peak_counts)

        results.append({
            'window_size': window,
            'median': float(np.median(peak_counts)),
            'mean': float(np.mean(peak_counts)),
            'q25': float(np.percentile(peak_counts, 25)),
            'q75': float(np.percentile(peak_counts, 75)),
            'min': int(np.min(peak_counts)),
            'max': int(np.max(peak_counts)),
            'n_zero': int(np.sum(peak_counts == 0)),
            'n_one': int(np.sum(peak_counts == 1)),
            'n_multiple': int(np.sum(peak_counts > 1))
        })

        print(f"  Window {window:>4}bp: median={results[-1]['median']:.2f}, "
              f"0 peaks: {results[-1]['n_zero']}, "
              f"1 peak: {results[-1]['n_one']}, "
              f">1 peaks: {results[-1]['n_multiple']}")

    return results


def main():
    # Paths
    hsa_step2_gtf = REPO_ROOT / "data" / "processed" / "tusco" / "hsa" / "step2_single_isoform.gtf.gz"
    hsa_gencode_gtf = REPO_ROOT / "data" / "raw" / "reference" / "human" / "gencode.v49.annotation.gtf.gz"
    hsa_reftss_bed = REPO_ROOT / "src" / "tusco_selector" / "data" / "hsa" / "splicing_tss" / "refTSS.v4.1.human.hg38.bed.gz"

    mmu_step2_gtf = REPO_ROOT / "data" / "processed" / "tusco" / "mmu" / "step2_single_isoform.gtf.gz"
    mmu_gencode_gtf = REPO_ROOT / "data" / "raw" / "reference" / "mouse" / "gencode.vM38.annotation.gtf.gz"
    mmu_reftss_bed = REPO_ROOT / "src" / "tusco_selector" / "data" / "mmu" / "splicing_tss" / "refTSS.v4.1.mouse.mm39.bed.gz"

    output_dir = Path(__file__).parent

    # Process human
    print("\n" + "="*80)
    print("HUMAN ANALYSIS")
    print("="*80)

    hsa_gene_tss = load_step2_genes(hsa_step2_gtf, "human")
    hsa_all_gencode = load_all_gencode_genes(hsa_gencode_gtf, "human")
    hsa_peaks = load_reftss_peaks(hsa_reftss_bed, "human")

    hsa_within = calculate_within_gene_distances(hsa_gene_tss, hsa_peaks, max_window=3000)
    hsa_between = calculate_between_gene_distances(hsa_gene_tss, hsa_all_gencode)
    hsa_peak_counts = analyze_peak_count_by_window(hsa_gene_tss, hsa_peaks)

    # Process mouse
    print("\n" + "="*80)
    print("MOUSE ANALYSIS")
    print("="*80)

    mmu_gene_tss = load_step2_genes(mmu_step2_gtf, "mouse")
    mmu_all_gencode = load_all_gencode_genes(mmu_gencode_gtf, "mouse")
    mmu_peaks = load_reftss_peaks(mmu_reftss_bed, "mouse")

    mmu_within = calculate_within_gene_distances(mmu_gene_tss, mmu_peaks, max_window=3000)
    mmu_between = calculate_between_gene_distances(mmu_gene_tss, mmu_all_gencode)
    mmu_peak_counts = analyze_peak_count_by_window(mmu_gene_tss, mmu_peaks)

    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)

    # Within-gene distances
    for species, within in [("human", hsa_within), ("mouse", mmu_within)]:
        output_file = output_dir / f"tss_within_gene_distances_corrected_{species}.tsv"
        with open(output_file, 'w') as f:
            f.write("gene_id\tdistance_to_reftss\n")
            for gene_id, distance in sorted(within, key=lambda x: x[1]):
                f.write(f"{gene_id}\t{distance}\n")
        print(f"  Saved {species} within-gene: {output_file}")

    # Between-gene distances
    for species, between in [("human", hsa_between), ("mouse", mmu_between)]:
        output_file = output_dir / f"tss_between_gene_distances_corrected_{species}.tsv"
        with open(output_file, 'w') as f:
            f.write("gene_id\tdistance_to_other_gene\n")
            for gene_id, distance in sorted(between, key=lambda x: x[1]):
                f.write(f"{gene_id}\t{distance}\n")
        print(f"  Saved {species} between-gene: {output_file}")

    # Peak counts by window
    for species, peak_counts in [("human", hsa_peak_counts), ("mouse", mmu_peak_counts)]:
        output_file = output_dir / f"reftss_peak_counts_by_window_{species}.tsv"
        with open(output_file, 'w') as f:
            f.write("window_size\tmedian\tmean\tq25\tq75\tmin\tmax\tn_zero\tn_one\tn_multiple\n")
            for row in peak_counts:
                f.write(f"{row['window_size']}\t{row['median']:.4f}\t{row['mean']:.4f}\t"
                       f"{row['q25']:.4f}\t{row['q75']:.4f}\t{row['min']}\t{row['max']}\t"
                       f"{row['n_zero']}\t{row['n_one']}\t{row['n_multiple']}\n")
        print(f"  Saved {species} peak counts: {output_file}")

    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS AT 300BP THRESHOLD")
    print("="*80)

    for species, within, between in [("Human", hsa_within, hsa_between), ("Mouse", mmu_within, mmu_between)]:
        print(f"\n{species.upper()}:")

        within_distances = [d for _, d in within]
        between_distances = [d for _, d in between]

        # Within-gene stats
        n_within_300 = sum(1 for d in within_distances if d <= 300)
        pct_within_300 = 100 * n_within_300 / len(within_distances)
        print(f"  Within-gene (to refTSS):")
        print(f"    ≤ 300bp: {n_within_300:,}/{len(within_distances):,} ({pct_within_300:.1f}%)")
        print(f"    Median: {np.median(within_distances):.0f} bp")

        # Between-gene stats
        n_between_300 = sum(1 for d in between_distances if d <= 300)
        pct_between_300 = 100 * n_between_300 / len(between_distances)
        print(f"  Between-gene (to other genes):")
        print(f"    ≤ 300bp: {n_between_300:,}/{len(between_distances):,} ({pct_between_300:.1f}%)")
        print(f"    Median: {np.median(between_distances):,.0f} bp")

    print("\n" + "="*80)


if __name__ == "__main__":
    main()
