#!/usr/bin/env python3
"""
Sensitivity Analysis Configuration Profiles for TUSCO Parameter Thresholds.

This module generates configuration profiles for running the tusco_selector pipeline
with different threshold settings to address reviewer concerns about parameter sensitivity.

Three independent sweeps:
1. Splice Sweep: Vary novel_threshold in Step 3 + splice_ratio_threshold in Alphagenome
2. TSS Sweep: Vary TSS filtering (off / lenient / default / strict) in Step 3
3. Expression Sweep: Vary expression thresholds in Step 5 + Alphagenome
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]


@dataclass
class TUSCOConfig:
    """Configuration for a single TUSCO pipeline run."""

    name: str
    species: str  # "hsa" or "mmu"
    sweep_type: str  # "splice", "tss", or "expression"
    is_default: bool = False

    # Step 3: Splice junction thresholds
    novel_threshold: float = 0.01  # hsa default; mmu default is 0.05
    min_novel_length: int = 80

    # Step 3: TSS thresholds
    tss_scope: str = "single"  # "single", "both", or "none"
    tss_region_check: int = 300  # bp window for CAGE peaks

    # Step 5: Expression filtering - GTEx (human)
    gtex_prevalence_expression_cutoff: float = 0.1
    gtex_median_expression_cutoff: float = 1.0
    gtex_prevalence_threshold: float = 0.95

    # Step 5: Expression filtering - Bgee
    # Human uses: gold+silver+bronze, prevalence=0.95
    # Mouse uses: gold+silver only, prevalence=0.90
    min_genes_per_tissue: int = 25000
    bgee_quality: list[str] = field(default_factory=lambda: ["gold quality", "silver quality", "bronze quality"])
    bgee_prevalence_threshold: float = 0.95  # Human default; mouse overrides to 0.90

    # Step 5: Alphagenome Universal
    single_exon_median_rpkm_threshold: float = 1.0
    multi_exon_median_rpkm_threshold: float = 1.0
    single_exon_expression_threshold: float = 0.05
    multi_exon_expression_threshold: float = 0.05
    single_exon_prevalence_threshold: float = 0.95
    multi_exon_prevalence_threshold: float = 0.95
    splice_ratio_threshold: float = 0.001  # hsa default; mmu default is 0.01
    splice_purity_threshold: float = 0.95

    # Step 5: Alphagenome Tissue-specific
    tissue_expression_threshold: float = 2.0
    tissue_splice_ratio_threshold: float = 0.001

    def to_cli_args(self) -> list[str]:
        """Convert config to CLI arguments for tusco_selector."""
        args = [
            self.species,
            "--output-dir", f"reviewer_response/sensitivity_runs/{self.name}",
            # Step 3: Splice
            "--novel-threshold", str(self.novel_threshold),
            "--min-novel-length", str(self.min_novel_length),
            # Step 3: TSS
            "--tss-scope", self.tss_scope,
        ]

        if self.tss_region_check > 0:
            args.extend(["--tss-region-check", str(self.tss_region_check)])

        if self.species == "hsa":
            args.extend([
                "--expression-source", "auto",
                # Bgee parameters
                "--min-genes-per-tissue", str(self.min_genes_per_tissue),
                "--bgee-quality", *self.bgee_quality,
                "--bgee-prevalence-threshold", str(self.bgee_prevalence_threshold),
                # GTEx parameters
                "--gtex-prevalence-expression-cutoff", str(self.gtex_prevalence_expression_cutoff),
                "--gtex-median-expression-cutoff", str(self.gtex_median_expression_cutoff),
                "--gtex-prevalence-threshold", str(self.gtex_prevalence_threshold),
            ])
        else:  # mmu
            args.extend([
                "--expression-source", "bgee",
                "--min-genes-per-tissue", str(self.min_genes_per_tissue),
                "--bgee-quality", *self.bgee_quality,
                "--bgee-prevalence-threshold", str(self.bgee_prevalence_threshold),
            ])

        # Alphagenome Universal
        args.extend([
            "--single-exon-median-rpkm-threshold", str(self.single_exon_median_rpkm_threshold),
            "--multi-exon-median-rpkm-threshold", str(self.multi_exon_median_rpkm_threshold),
            "--single-exon-expression-threshold", str(self.single_exon_expression_threshold),
            "--multi-exon-expression-threshold", str(self.multi_exon_expression_threshold),
            "--single-exon-prevalence-threshold", str(self.single_exon_prevalence_threshold),
            "--multi-exon-prevalence-threshold", str(self.multi_exon_prevalence_threshold),
            "--splice-ratio-threshold", str(self.splice_ratio_threshold),
            "--splice-purity-threshold", str(self.splice_purity_threshold),
            # Alphagenome Tissue-specific
            "--tissue-expression-threshold", str(self.tissue_expression_threshold),
            "--tissue-splice-ratio-threshold", str(self.tissue_splice_ratio_threshold),
        ])

        return args


# =============================================================================
# DEFAULT CONFIGURATIONS (ORIGINAL THRESHOLDS)
# =============================================================================

def get_hsa_default() -> TUSCOConfig:
    """Human default configuration (ORIGINAL)."""
    return TUSCOConfig(
        name="hsa_default",
        species="hsa",
        sweep_type="default",
        is_default=True,
        # Step 3
        novel_threshold=0.01,
        tss_scope="single",
        tss_region_check=300,
        # Step 5 - Human uses defaults: bgee_quality=gold+silver+bronze, bgee_prevalence_threshold=0.95
        gtex_median_expression_cutoff=1.0,
        gtex_prevalence_expression_cutoff=0.1,
        # Alphagenome
        single_exon_median_rpkm_threshold=1.0,
        multi_exon_median_rpkm_threshold=1.0,
        single_exon_expression_threshold=0.05,
        multi_exon_expression_threshold=0.05,
        splice_ratio_threshold=0.001,
    )


def get_mmu_default() -> TUSCOConfig:
    """Mouse default configuration (ORIGINAL)."""
    return TUSCOConfig(
        name="mmu_default",
        species="mmu",
        sweep_type="default",
        is_default=True,
        # Step 3
        novel_threshold=0.05,
        tss_scope="single",
        tss_region_check=0,  # Mouse uses exon-overlap only, no CAGE window check (unlike human)
        # Step 5 - Mouse uses gold+silver only (no bronze)
        bgee_quality=["gold quality", "silver quality"],
        bgee_prevalence_threshold=0.90,
        # Alphagenome
        single_exon_median_rpkm_threshold=0.5,
        multi_exon_median_rpkm_threshold=0.25,
        single_exon_expression_threshold=0.01,
        multi_exon_expression_threshold=0.005,
        splice_ratio_threshold=0.01,
        single_exon_prevalence_threshold=0.90,
        multi_exon_prevalence_threshold=0.90,
        splice_purity_threshold=0.90,
        tissue_expression_threshold=1.5,
    )


# =============================================================================
# SWEEP A: SPLICE JUNCTION THRESHOLD
# =============================================================================

def get_splice_sweep_hsa() -> list[TUSCOConfig]:
    """Human splice threshold sweep configurations."""
    base = get_hsa_default()

    configs = [
        # Strict: 100x stricter - minimal genes
        TUSCOConfig(
            name="hsa_splice_strict",
            species="hsa",
            sweep_type="splice",
            novel_threshold=0.0001,
            splice_ratio_threshold=0.00001,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
            gtex_median_expression_cutoff=base.gtex_median_expression_cutoff,
            gtex_prevalence_expression_cutoff=base.gtex_prevalence_expression_cutoff,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
        ),
        # Default (included for reference)
        TUSCOConfig(
            name="hsa_splice_default",
            species="hsa",
            sweep_type="splice",
            is_default=True,
            novel_threshold=0.01,
            splice_ratio_threshold=0.001,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
            gtex_median_expression_cutoff=base.gtex_median_expression_cutoff,
            gtex_prevalence_expression_cutoff=base.gtex_prevalence_expression_cutoff,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
        ),
        # Lenient 1: 50x more lenient
        TUSCOConfig(
            name="hsa_splice_lenient1",
            species="hsa",
            sweep_type="splice",
            novel_threshold=0.50,
            splice_ratio_threshold=0.25,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
            gtex_median_expression_cutoff=base.gtex_median_expression_cutoff,
            gtex_prevalence_expression_cutoff=base.gtex_prevalence_expression_cutoff,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
        ),
        # Lenient 2: 75x more lenient
        TUSCOConfig(
            name="hsa_splice_lenient2",
            species="hsa",
            sweep_type="splice",
            novel_threshold=0.75,
            splice_ratio_threshold=0.50,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
            gtex_median_expression_cutoff=base.gtex_median_expression_cutoff,
            gtex_prevalence_expression_cutoff=base.gtex_prevalence_expression_cutoff,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
        ),
        # Lenient 3: Allow all splicing - maximum genes
        TUSCOConfig(
            name="hsa_splice_lenient3",
            species="hsa",
            sweep_type="splice",
            novel_threshold=1.0,
            splice_ratio_threshold=0.75,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
            gtex_median_expression_cutoff=base.gtex_median_expression_cutoff,
            gtex_prevalence_expression_cutoff=base.gtex_prevalence_expression_cutoff,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
        ),
    ]
    return configs


def get_splice_sweep_mmu() -> list[TUSCOConfig]:
    """Mouse splice threshold sweep configurations."""
    base = get_mmu_default()

    configs = [
        # Strict: 50x stricter
        TUSCOConfig(
            name="mmu_splice_strict",
            species="mmu",
            sweep_type="splice",
            novel_threshold=0.001,
            splice_ratio_threshold=0.0001,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=base.bgee_prevalence_threshold,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
            single_exon_prevalence_threshold=base.single_exon_prevalence_threshold,
            multi_exon_prevalence_threshold=base.multi_exon_prevalence_threshold,
            splice_purity_threshold=base.splice_purity_threshold,
        ),
        # Default
        TUSCOConfig(
            name="mmu_splice_default",
            species="mmu",
            sweep_type="splice",
            is_default=True,
            novel_threshold=0.05,
            splice_ratio_threshold=0.01,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=base.bgee_prevalence_threshold,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
            single_exon_prevalence_threshold=base.single_exon_prevalence_threshold,
            multi_exon_prevalence_threshold=base.multi_exon_prevalence_threshold,
            splice_purity_threshold=base.splice_purity_threshold,
        ),
        # Lenient 1: 10x more lenient
        TUSCOConfig(
            name="mmu_splice_lenient1",
            species="mmu",
            sweep_type="splice",
            novel_threshold=0.50,
            splice_ratio_threshold=0.25,
            tss_scope=base.tss_scope,
            bgee_quality=base.bgee_quality,
            tss_region_check=base.tss_region_check,
            bgee_prevalence_threshold=base.bgee_prevalence_threshold,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
            single_exon_prevalence_threshold=base.single_exon_prevalence_threshold,
            multi_exon_prevalence_threshold=base.multi_exon_prevalence_threshold,
            splice_purity_threshold=base.splice_purity_threshold,
        ),
        # Lenient 2: 15x more lenient
        TUSCOConfig(
            name="mmu_splice_lenient2",
            species="mmu",
            sweep_type="splice",
            novel_threshold=0.75,
            splice_ratio_threshold=0.50,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=base.bgee_prevalence_threshold,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
            single_exon_prevalence_threshold=base.single_exon_prevalence_threshold,
            multi_exon_prevalence_threshold=base.multi_exon_prevalence_threshold,
            splice_purity_threshold=base.splice_purity_threshold,
        ),
        # Lenient 3: Allow all splicing - maximum genes
        TUSCOConfig(
            name="mmu_splice_lenient3",
            species="mmu",
            sweep_type="splice",
            novel_threshold=1.0,
            splice_ratio_threshold=0.75,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=base.bgee_prevalence_threshold,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
            single_exon_prevalence_threshold=base.single_exon_prevalence_threshold,
            multi_exon_prevalence_threshold=base.multi_exon_prevalence_threshold,
            splice_purity_threshold=base.splice_purity_threshold,
        ),
    ]
    return configs


# =============================================================================
# SWEEP B: TSS THRESHOLD
# =============================================================================

def get_tss_sweep_hsa() -> list[TUSCOConfig]:
    """Human TSS threshold sweep configurations."""
    base = get_hsa_default()

    configs = [
        # Lenient: 10bp window (very small = most genes pass)
        TUSCOConfig(
            name="hsa_tss_lenient",
            species="hsa",
            sweep_type="tss",
            tss_scope="single",
            tss_region_check=10,
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            gtex_median_expression_cutoff=base.gtex_median_expression_cutoff,
            gtex_prevalence_expression_cutoff=base.gtex_prevalence_expression_cutoff,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
        ),
        # Default: 300bp window
        TUSCOConfig(
            name="hsa_tss_default",
            species="hsa",
            sweep_type="tss",
            is_default=True,
            tss_scope="single",
            tss_region_check=300,
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            gtex_median_expression_cutoff=base.gtex_median_expression_cutoff,
            gtex_prevalence_expression_cutoff=base.gtex_prevalence_expression_cutoff,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
        ),
        # Strict: 2000bp window (large window = fewer genes pass)
        TUSCOConfig(
            name="hsa_tss_strict",
            species="hsa",
            sweep_type="tss",
            tss_scope="single",
            tss_region_check=2000,
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            gtex_median_expression_cutoff=base.gtex_median_expression_cutoff,
            gtex_prevalence_expression_cutoff=base.gtex_prevalence_expression_cutoff,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
        ),
        # Extreme: 5000bp window (very large = minimal genes)
        TUSCOConfig(
            name="hsa_tss_extreme",
            species="hsa",
            sweep_type="tss",
            tss_scope="single",
            tss_region_check=5000,
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            gtex_median_expression_cutoff=base.gtex_median_expression_cutoff,
            gtex_prevalence_expression_cutoff=base.gtex_prevalence_expression_cutoff,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
        ),
    ]
    return configs


def get_tss_sweep_mmu() -> list[TUSCOConfig]:
    """Mouse TSS threshold sweep configurations."""
    base = get_mmu_default()

    configs = [
        # Lenient: 10bp window (very small = most genes pass)
        TUSCOConfig(
            name="mmu_tss_lenient",
            species="mmu",
            sweep_type="tss",
            tss_scope="single",
            tss_region_check=10,
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=base.bgee_prevalence_threshold,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
            single_exon_prevalence_threshold=base.single_exon_prevalence_threshold,
            multi_exon_prevalence_threshold=base.multi_exon_prevalence_threshold,
            splice_purity_threshold=base.splice_purity_threshold,
        ),
        # Default: 300bp window
        TUSCOConfig(
            name="mmu_tss_default",
            species="mmu",
            sweep_type="tss",
            is_default=True,
            tss_scope="single",
            tss_region_check=300,
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=base.bgee_prevalence_threshold,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
            single_exon_prevalence_threshold=base.single_exon_prevalence_threshold,
            multi_exon_prevalence_threshold=base.multi_exon_prevalence_threshold,
            splice_purity_threshold=base.splice_purity_threshold,
        ),
        # Strict: 2000bp window (large window = fewer genes pass)
        TUSCOConfig(
            name="mmu_tss_strict",
            species="mmu",
            sweep_type="tss",
            tss_scope="single",
            tss_region_check=2000,
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=base.bgee_prevalence_threshold,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
            single_exon_prevalence_threshold=base.single_exon_prevalence_threshold,
            multi_exon_prevalence_threshold=base.multi_exon_prevalence_threshold,
            splice_purity_threshold=base.splice_purity_threshold,
        ),
        # Extreme: 5000bp window (very large = minimal genes)
        TUSCOConfig(
            name="mmu_tss_extreme",
            species="mmu",
            sweep_type="tss",
            tss_scope="single",
            tss_region_check=5000,
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=base.bgee_prevalence_threshold,
            single_exon_median_rpkm_threshold=base.single_exon_median_rpkm_threshold,
            multi_exon_median_rpkm_threshold=base.multi_exon_median_rpkm_threshold,
            single_exon_expression_threshold=base.single_exon_expression_threshold,
            multi_exon_expression_threshold=base.multi_exon_expression_threshold,
            single_exon_prevalence_threshold=base.single_exon_prevalence_threshold,
            multi_exon_prevalence_threshold=base.multi_exon_prevalence_threshold,
            splice_purity_threshold=base.splice_purity_threshold,
        ),
    ]
    return configs


# =============================================================================
# SWEEP C: EXPRESSION THRESHOLD
# =============================================================================

def get_expression_sweep_hsa() -> list[TUSCOConfig]:
    """Human expression threshold sweep configurations."""
    base = get_hsa_default()

    configs = [
        # Lenient 2: 100x more lenient - maximum genes
        TUSCOConfig(
            name="hsa_expr_lenient2",
            species="hsa",
            sweep_type="expression",
            gtex_median_expression_cutoff=0.01,
            gtex_prevalence_expression_cutoff=0.001,
            single_exon_median_rpkm_threshold=0.01,
            multi_exon_median_rpkm_threshold=0.01,
            single_exon_expression_threshold=0.001,
            multi_exon_expression_threshold=0.001,
            # Keep other params at default
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
        ),
        # Lenient 1: 10x more lenient
        TUSCOConfig(
            name="hsa_expr_lenient1",
            species="hsa",
            sweep_type="expression",
            gtex_median_expression_cutoff=0.1,
            gtex_prevalence_expression_cutoff=0.01,
            single_exon_median_rpkm_threshold=0.1,
            multi_exon_median_rpkm_threshold=0.1,
            single_exon_expression_threshold=0.01,
            multi_exon_expression_threshold=0.01,
            # Keep other params at default
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
        ),
        # Default
        TUSCOConfig(
            name="hsa_expr_default",
            species="hsa",
            sweep_type="expression",
            is_default=True,
            gtex_median_expression_cutoff=1.0,
            gtex_prevalence_expression_cutoff=0.1,
            single_exon_median_rpkm_threshold=1.0,
            multi_exon_median_rpkm_threshold=1.0,
            single_exon_expression_threshold=0.05,
            multi_exon_expression_threshold=0.05,
            # Keep other params at default
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
        ),
        # Strict: 10x stricter - fewer genes
        TUSCOConfig(
            name="hsa_expr_strict",
            species="hsa",
            sweep_type="expression",
            gtex_median_expression_cutoff=10.0,
            gtex_prevalence_expression_cutoff=1.0,
            single_exon_median_rpkm_threshold=10.0,
            multi_exon_median_rpkm_threshold=10.0,
            single_exon_expression_threshold=0.5,
            multi_exon_expression_threshold=0.5,
            # Keep other params at default
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
        ),
    ]
    return configs


def get_expression_sweep_mmu() -> list[TUSCOConfig]:
    """Mouse expression threshold sweep configurations."""
    base = get_mmu_default()

    configs = [
        # Lenient 2: Much more lenient - maximum genes
        TUSCOConfig(
            name="mmu_expr_lenient2",
            species="mmu",
            sweep_type="expression",
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=0.50,
            single_exon_median_rpkm_threshold=0.01,
            multi_exon_median_rpkm_threshold=0.01,
            single_exon_expression_threshold=0.0001,
            multi_exon_expression_threshold=0.0001,
            single_exon_prevalence_threshold=0.50,
            multi_exon_prevalence_threshold=0.50,
            splice_purity_threshold=0.50,
            # Keep other params at default
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
        ),
        # Lenient 1: More lenient
        TUSCOConfig(
            name="mmu_expr_lenient1",
            species="mmu",
            sweep_type="expression",
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=0.75,
            single_exon_median_rpkm_threshold=0.1,
            multi_exon_median_rpkm_threshold=0.05,
            single_exon_expression_threshold=0.001,
            multi_exon_expression_threshold=0.001,
            single_exon_prevalence_threshold=0.75,
            multi_exon_prevalence_threshold=0.75,
            splice_purity_threshold=0.75,
            # Keep other params at default
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
        ),
        # Default
        TUSCOConfig(
            name="mmu_expr_default",
            species="mmu",
            sweep_type="expression",
            is_default=True,
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=0.90,
            single_exon_median_rpkm_threshold=0.5,
            multi_exon_median_rpkm_threshold=0.25,
            single_exon_expression_threshold=0.01,
            multi_exon_expression_threshold=0.005,
            single_exon_prevalence_threshold=0.90,
            multi_exon_prevalence_threshold=0.90,
            splice_purity_threshold=0.90,
            # Keep other params at default
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
        ),
        # Strict - fewer genes
        TUSCOConfig(
            name="mmu_expr_strict",
            species="mmu",
            sweep_type="expression",
            bgee_quality=base.bgee_quality,
            bgee_prevalence_threshold=0.99,
            single_exon_median_rpkm_threshold=5.0,
            multi_exon_median_rpkm_threshold=2.5,
            single_exon_expression_threshold=0.1,
            multi_exon_expression_threshold=0.05,
            single_exon_prevalence_threshold=0.99,
            multi_exon_prevalence_threshold=0.99,
            splice_purity_threshold=0.99,
            # Keep other params at default
            novel_threshold=base.novel_threshold,
            splice_ratio_threshold=base.splice_ratio_threshold,
            tss_scope=base.tss_scope,
            tss_region_check=base.tss_region_check,
        ),
    ]
    return configs


# =============================================================================
# GET ALL CONFIGURATIONS
# =============================================================================

def get_all_configs() -> list[TUSCOConfig]:
    """Get all sensitivity analysis configurations."""
    configs = []

    # Splice sweeps
    configs.extend(get_splice_sweep_hsa())
    configs.extend(get_splice_sweep_mmu())

    # TSS sweeps
    configs.extend(get_tss_sweep_hsa())
    configs.extend(get_tss_sweep_mmu())

    # Expression sweeps
    configs.extend(get_expression_sweep_hsa())
    configs.extend(get_expression_sweep_mmu())

    return configs


def get_unique_configs() -> list[TUSCOConfig]:
    """Get unique configurations (remove duplicate defaults)."""
    configs = get_all_configs()
    seen_names = set()
    unique = []

    for cfg in configs:
        if cfg.name not in seen_names:
            seen_names.add(cfg.name)
            unique.append(cfg)

    return unique


def export_to_json(output_path: Optional[Path] = None) -> Path:
    """Export all configurations to JSON."""
    if output_path is None:
        output_path = REPO_ROOT / "scripts" / "sensitivity" / "launch_configs.json"

    configs = get_unique_configs()

    data = {
        "configs": [
            {
                "name": cfg.name,
                "species": cfg.species,
                "sweep_type": cfg.sweep_type,
                "is_default": cfg.is_default,
                "cli_args": cfg.to_cli_args(),
            }
            for cfg in configs
        ],
        "summary": {
            "total_configs": len(configs),
            "splice_sweep": len([c for c in configs if c.sweep_type == "splice"]),
            "tss_sweep": len([c for c in configs if c.sweep_type == "tss"]),
            "expression_sweep": len([c for c in configs if c.sweep_type == "expression"]),
            "human_configs": len([c for c in configs if c.species == "hsa"]),
            "mouse_configs": len([c for c in configs if c.species == "mmu"]),
        }
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        json.dump(data, f, indent=2)

    print(f"Exported {len(configs)} configurations to {output_path}")
    return output_path


def main():
    """Main entry point - export configurations."""
    export_to_json()

    # Print summary
    configs = get_unique_configs()
    print(f"\nSensitivity Analysis Configuration Summary:")
    print(f"  Total configurations: {len(configs)}")
    print(f"  - Splice sweep: {len([c for c in configs if c.sweep_type == 'splice'])}")
    print(f"  - TSS sweep: {len([c for c in configs if c.sweep_type == 'tss'])}")
    print(f"  - Expression sweep: {len([c for c in configs if c.sweep_type == 'expression'])}")
    print(f"  - Human (hsa): {len([c for c in configs if c.species == 'hsa'])}")
    print(f"  - Mouse (mmu): {len([c for c in configs if c.species == 'mmu'])}")

    # Print individual configs
    print("\nConfigurations:")
    for cfg in configs:
        marker = " (DEFAULT)" if cfg.is_default else ""
        print(f"  {cfg.name}: {cfg.sweep_type}{marker}")


if __name__ == "__main__":
    main()
