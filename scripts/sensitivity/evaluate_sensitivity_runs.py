#!/usr/bin/env python3
"""
Evaluate TUSCO sensitivity analysis runs on LRGASP classification data.

This script:
1. Loads gene sets from sensitivity run outputs
2. Evaluates them on LRGASP WTC-11 (human) and mouse ES classifications
3. Computes FP, FN, TP, PTP metrics with breakdown by category
4. Outputs evaluation results for visualization

Usage:
    python evaluate_sensitivity_runs.py [--runs-dir DIR] [--output-dir DIR]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# DATA LOADING
# =============================================================================

def read_gene_ids_from_tsv(tsv_path: Path) -> set[str]:
    """Read gene IDs from a TUSCO output TSV file."""
    gene_ids = set()
    if not tsv_path.exists():
        return gene_ids

    with tsv_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # First column is gene ID, strip version number
            gene_id = line.split("\t")[0].split(".")[0]
            gene_ids.add(gene_id)

    return gene_ids


def load_lrgasp_classification(path: Path) -> pd.DataFrame:
    """Load LRGASP classification file and clean it."""
    df = pd.read_csv(path, sep="\t", low_memory=False)

    # Clean gene and transcript IDs (remove version numbers)
    for col in ("associated_gene", "associated_transcript"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.replace(r"\.\d+$", "", regex=True)

    # Handle fusion rows (expand to multiple rows)
    if "structural_category" in df.columns:
        is_fusion = df["structural_category"].astype(str) == "fusion"
        if is_fusion.any():
            fusion = df.loc[is_fusion].copy()
            non = df.loc[~is_fusion].copy()

            if "associated_gene" in fusion.columns:
                fusion["associated_gene"] = fusion["associated_gene"].astype(str).str.split("_")
                fusion = fusion.explode("associated_gene", ignore_index=True)
            if "associated_transcript" in fusion.columns:
                fusion["associated_transcript"] = fusion["associated_transcript"].astype(str).str.split("_")
                fusion = fusion.explode("associated_transcript", ignore_index=True)

            df = pd.concat([non, fusion], ignore_index=True)

    # Ensure numeric columns
    for col in ("ref_exons", "diff_to_TSS", "diff_to_TTS", "ref_length", "min_cov", "exons"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def load_sirv_reference(gtf_path: Path) -> set[str]:
    """Load SIRV reference transcript IDs from GTF."""
    tx = set()
    with gtf_path.open() as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attrs = parts[8]
            key = 'transcript_id "'
            idx = attrs.find(key)
            if idx < 0:
                continue
            start = idx + len(key)
            end = attrs.find('"', start)
            if end >= 0:
                tx.add(attrs[start:end])
    return tx


# =============================================================================
# CLASSIFICATION RULES
# =============================================================================

def is_true_positive(row: pd.Series, monoexon_window: int = 50, long_window: int = 100) -> bool:
    """Check if a row is a True Positive based on classification rules.

    Rules:
    1. subcategory == "reference_match" -> TP
    2. FSM + single-exon (ref_exons==1) + TSS/TTS within ±monoexon_window -> TP
    3. FSM + long gene (ref_length > 3000) + TSS/TTS within ±long_window -> TP
    """
    subcat = str(row.get("subcategory", ""))
    structural = str(row.get("structural_category", ""))

    # Rule 1: Reference match
    if subcat == "reference_match":
        return True

    # Rules 2 & 3: FSM with specific criteria
    if structural != "full-splice_match":
        return False

    diff_tss = row.get("diff_to_TSS")
    diff_tts = row.get("diff_to_TTS")

    if pd.isna(diff_tss) or pd.isna(diff_tts):
        return False

    ref_exons = row.get("ref_exons")
    ref_length = row.get("ref_length")

    # Rule 2: Single-exon
    if ref_exons == 1:
        if abs(diff_tss) <= monoexon_window and abs(diff_tts) <= monoexon_window:
            return True

    # Rule 3: Long gene
    if not pd.isna(ref_length) and ref_length > 3000:
        if abs(diff_tss) <= long_window and abs(diff_tts) <= long_window:
            return True

    return False


def classify_transcript(row: pd.Series, tp_isoforms: set[str]) -> str:
    """Classify a transcript row into TP, PTP, or FP category.

    Returns one of: "TP", "PTP", "FP_NIC", "FP_NNC", "FP_genic", "FP_antisense",
                    "FP_intergenic", "FP_fusion", "FP_genic_intron", "FP_other"
    """
    isoform = str(row.get("isoform", ""))
    structural = str(row.get("structural_category", ""))
    subcat = str(row.get("subcategory", ""))

    # TP: Already classified as true positive
    if isoform in tp_isoforms:
        return "TP"

    # PTP: FSM or ISM not in TP, or mono-exon by intron retention
    if structural in ("full-splice_match", "incomplete-splice_match") or subcat == "mono-exon_by_intron_retention":
        return "PTP"

    # FP categories
    fp_mapping = {
        "novel_in_catalog": "FP_NIC",
        "novel_not_in_catalog": "FP_NNC",
        "genic": "FP_genic",
        "antisense": "FP_antisense",
        "intergenic": "FP_intergenic",
        "fusion": "FP_fusion",
        "genic_intron": "FP_genic_intron",
    }

    if subcat != "mono-exon_by_intron_retention":
        for cat, label in fp_mapping.items():
            if structural == cat:
                return label

    return "FP_other"


# =============================================================================
# EVALUATION
# =============================================================================

@dataclass
class EvaluationResult:
    """Results from evaluating a gene set on LRGASP data."""
    config_name: str
    species: str
    pipeline: str
    sweep_type: str

    # Gene set info
    n_genes: int = 0
    n_single_exon: int = 0
    n_multi_exon: int = 0

    # Classification counts
    n_observed: int = 0
    n_tp: int = 0
    n_ptp: int = 0
    n_fn: int = 0

    # FP breakdown
    n_fp_total: int = 0
    n_fp_nic: int = 0  # Novel in catalog
    n_fp_nnc: int = 0  # Novel not in catalog
    n_fp_genic: int = 0
    n_fp_antisense: int = 0
    n_fp_intergenic: int = 0
    n_fp_fusion: int = 0
    n_fp_other: int = 0

    # Metrics
    sensitivity: float = 0.0  # TP unique genes / reference genes
    precision: float = 0.0  # TP / observed
    redundant_precision: float = 0.0  # (TP + PTP) / observed
    false_discovery_rate: float = 0.0  # (observed - TP) / observed
    false_detection_rate: float = 0.0  # FP / observed

    def to_dict(self) -> dict:
        return {
            "config_name": self.config_name,
            "species": self.species,
            "pipeline": self.pipeline,
            "sweep_type": self.sweep_type,
            "n_genes": self.n_genes,
            "n_single_exon": self.n_single_exon,
            "n_multi_exon": self.n_multi_exon,
            "n_observed": self.n_observed,
            "n_tp": self.n_tp,
            "n_ptp": self.n_ptp,
            "n_fn": self.n_fn,
            "n_fp_total": self.n_fp_total,
            "n_fp_nic": self.n_fp_nic,
            "n_fp_nnc": self.n_fp_nnc,
            "n_fp_genic": self.n_fp_genic,
            "n_fp_antisense": self.n_fp_antisense,
            "n_fp_intergenic": self.n_fp_intergenic,
            "n_fp_fusion": self.n_fp_fusion,
            "n_fp_other": self.n_fp_other,
            "sensitivity": self.sensitivity,
            "precision": self.precision,
            "redundant_precision": self.redundant_precision,
            "false_discovery_rate": self.false_discovery_rate,
            "false_detection_rate": self.false_detection_rate,
        }


def evaluate_gene_set(
    gene_ids: set[str],
    classification_df: pd.DataFrame,
    config_name: str,
    species: str,
    pipeline: str,
    sweep_type: str,
    single_exon_ids: Optional[set[str]] = None,
    monoexon_window: int = 50,
    long_window: int = 100,
) -> EvaluationResult:
    """Evaluate a gene set on LRGASP classification data.

    Args:
        gene_ids: Set of gene IDs to evaluate
        classification_df: LRGASP classification dataframe
        config_name: Name of the configuration
        species: Species code ("hsa" or "mmu")
        pipeline: Pipeline name (e.g., "WTC11_drna_ont")
        sweep_type: Type of sweep ("splice", "tss", "expression")
        single_exon_ids: Optional set of single-exon gene IDs
        monoexon_window: Window size for TP classification (single-exon)
        long_window: Window size for TP classification (long genes)

    Returns:
        EvaluationResult with all metrics
    """
    result = EvaluationResult(
        config_name=config_name,
        species=species,
        pipeline=pipeline,
        sweep_type=sweep_type,
    )

    if not gene_ids:
        return result

    # Gene set info
    result.n_genes = len(gene_ids)
    if single_exon_ids:
        result.n_single_exon = len(gene_ids & single_exon_ids)
        result.n_multi_exon = result.n_genes - result.n_single_exon

    # Filter classification to our gene set
    df = classification_df[
        classification_df["associated_gene"].astype(str).isin(gene_ids)
    ].copy()

    if df.empty:
        result.n_fn = len(gene_ids)
        return result

    result.n_observed = len(df)

    # Identify TPs
    df["is_tp"] = df.apply(lambda row: is_true_positive(row, monoexon_window, long_window), axis=1)
    tp_isoforms = set(df.loc[df["is_tp"], "isoform"].astype(str))
    result.n_tp = int(df["is_tp"].sum())

    # Classify all transcripts
    df["classification"] = df.apply(lambda row: classify_transcript(row, tp_isoforms), axis=1)

    # Count categories
    class_counts = df["classification"].value_counts()
    result.n_ptp = int(class_counts.get("PTP", 0))
    result.n_fp_nic = int(class_counts.get("FP_NIC", 0))
    result.n_fp_nnc = int(class_counts.get("FP_NNC", 0))
    result.n_fp_genic = int(class_counts.get("FP_genic", 0))
    result.n_fp_antisense = int(class_counts.get("FP_antisense", 0))
    result.n_fp_intergenic = int(class_counts.get("FP_intergenic", 0))
    result.n_fp_fusion = int(class_counts.get("FP_fusion", 0))
    result.n_fp_other = int(class_counts.get("FP_other", 0) + class_counts.get("FP_genic_intron", 0))

    result.n_fp_total = (
        result.n_fp_nic + result.n_fp_nnc + result.n_fp_genic +
        result.n_fp_antisense + result.n_fp_intergenic + result.n_fp_fusion +
        result.n_fp_other
    )

    # FN: Genes in reference not detected at all
    detected_genes = set(df["associated_gene"].dropna().astype(str))
    result.n_fn = len(gene_ids - detected_genes)

    # Compute metrics
    if result.n_genes > 0:
        # Sensitivity: unique TP genes / reference genes
        tp_genes = set(df.loc[df["is_tp"], "associated_gene"].dropna().astype(str))
        result.sensitivity = len(tp_genes) / result.n_genes

    if result.n_observed > 0:
        result.precision = result.n_tp / result.n_observed
        result.redundant_precision = (result.n_tp + result.n_ptp) / result.n_observed
        result.false_discovery_rate = (result.n_observed - result.n_tp) / result.n_observed
        result.false_detection_rate = result.n_fp_total / result.n_observed

    return result


# =============================================================================
# MAIN EVALUATION
# =============================================================================

def evaluate_all_runs(
    runs_dir: Path,
    lrgasp_dir: Path,
    sirv_gtf: Path,
    output_dir: Path,
) -> list[EvaluationResult]:
    """Evaluate all sensitivity runs on LRGASP data.

    Args:
        runs_dir: Directory containing sensitivity run outputs
        lrgasp_dir: Directory containing LRGASP classification files
        sirv_gtf: Path to SIRV reference GTF
        output_dir: Directory to save results

    Returns:
        List of EvaluationResult objects
    """
    results = []

    # Define pipelines
    human_pipelines = [
        "WTC11_drna_ont",
        "WTC11_cdna_ont",
        "WTC11_cdna_pacbio",
        "WTC11_drna_ont_ls",
        "WTC11_cdna_ont_ls",
        "WTC11_cdna_pacbio_ls",
    ]
    mouse_pipelines = [
        "ES_drna_ont",
        "ES_cdna_ont",
        "ES_cdna_pacbio",
        "ES_drna_ont_ls",
        "ES_cdna_ont_ls",
        "ES_cdna_pacbio_ls",
    ]

    # Load LRGASP classifications
    print("Loading LRGASP classification files...")
    classifications = {}

    for pipeline in human_pipelines:
        path = lrgasp_dir / "human" / pipeline / f"{pipeline}_classification.txt"
        if path.exists():
            classifications[("hsa", pipeline)] = load_lrgasp_classification(path)
            print(f"  Loaded {pipeline} ({len(classifications[('hsa', pipeline)])} rows)")

    for pipeline in mouse_pipelines:
        path = lrgasp_dir / "mouse" / pipeline / f"{pipeline}_classification.txt"
        if path.exists():
            classifications[("mmu", pipeline)] = load_lrgasp_classification(path)
            print(f"  Loaded {pipeline} ({len(classifications[('mmu', pipeline)])} rows)")

    # Find all sensitivity runs
    print("\nEvaluating sensitivity runs...")
    run_dirs = sorted([d for d in runs_dir.iterdir() if d.is_dir() and not d.name.startswith(".")])

    for run_dir in run_dirs:
        config_name = run_dir.name

        # Determine species from config name
        if config_name.startswith("hsa"):
            species = "hsa"
            species_name = "human"
            pipelines = human_pipelines
        elif config_name.startswith("mmu"):
            species = "mmu"
            species_name = "mouse"
            pipelines = mouse_pipelines
        else:
            print(f"  Skipping {config_name} (unknown species)")
            continue

        # Determine sweep type
        if "_splice_" in config_name:
            sweep_type = "splice"
        elif "_tss_" in config_name:
            sweep_type = "tss"
        elif "_expr_" in config_name:
            sweep_type = "expression"
        else:
            sweep_type = "default"

        # Load gene set
        main_file = run_dir / f"tusco_{species_name}.tsv"
        single_exon_file = run_dir / f"tusco_{species_name}_single_exon.tsv"

        if not main_file.exists():
            print(f"  Skipping {config_name} (no output file)")
            continue

        gene_ids = read_gene_ids_from_tsv(main_file)
        single_exon_ids = read_gene_ids_from_tsv(single_exon_file) if single_exon_file.exists() else None

        print(f"  Evaluating {config_name} ({len(gene_ids)} genes)...")

        # Evaluate on each pipeline
        for pipeline in pipelines:
            if (species, pipeline) not in classifications:
                continue

            result = evaluate_gene_set(
                gene_ids=gene_ids,
                classification_df=classifications[(species, pipeline)],
                config_name=config_name,
                species=species,
                pipeline=pipeline,
                sweep_type=sweep_type,
                single_exon_ids=single_exon_ids,
            )
            results.append(result)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate TUSCO sensitivity analysis runs on LRGASP data"
    )
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=REPO_ROOT / "reviewer_response" / "sensitivity_runs",
        help="Directory containing sensitivity run outputs",
    )
    parser.add_argument(
        "--lrgasp-dir",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "lrgasp",
        help="Directory containing LRGASP classification files",
    )
    parser.add_argument(
        "--sirv-gtf",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "spike-ins" / "lrgasp_sirvs.gtf",
        help="Path to SIRV reference GTF",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to save results (default: runs-dir/summary)",
    )
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.runs_dir / "summary"

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Check inputs
    if not args.runs_dir.exists():
        print(f"Error: Runs directory not found: {args.runs_dir}")
        sys.exit(1)

    if not args.lrgasp_dir.exists():
        print(f"Error: LRGASP directory not found: {args.lrgasp_dir}")
        sys.exit(1)

    # Run evaluation
    results = evaluate_all_runs(
        runs_dir=args.runs_dir,
        lrgasp_dir=args.lrgasp_dir,
        sirv_gtf=args.sirv_gtf,
        output_dir=args.output_dir,
    )

    if not results:
        print("No results to save")
        sys.exit(0)

    # Save results as TSV
    results_file = args.output_dir / "evaluation_results.tsv"
    fieldnames = list(results[0].to_dict().keys())

    with results_file.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for result in results:
            writer.writerow(result.to_dict())

    print(f"\nResults saved to: {results_file}")
    print(f"Total evaluations: {len(results)}")

    # Print summary by config
    configs = sorted(set(r.config_name for r in results))
    print("\nSummary by configuration:")
    for config in configs:
        config_results = [r for r in results if r.config_name == config]
        avg_fp = sum(r.n_fp_total for r in config_results) / len(config_results)
        avg_fn = sum(r.n_fn for r in config_results) / len(config_results)
        avg_sensitivity = sum(r.sensitivity for r in config_results) / len(config_results)
        n_genes = config_results[0].n_genes

        print(f"  {config}: {n_genes} genes, avg FP={avg_fp:.1f}, avg FN={avg_fn:.1f}, "
              f"avg sensitivity={avg_sensitivity:.3f}")


if __name__ == "__main__":
    main()
