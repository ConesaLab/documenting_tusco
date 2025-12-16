"""
Step 1: Select single isoform per gene.
"""

from __future__ import annotations

import logging
import os
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

# Optional dependencies
from tusco_selector.optional_deps import HAS_TQDM, tqdm

# Setup module logger
logger = logging.getLogger(__name__)

# Shared utilities ----------------------------------------------------------
# Import once here so the heavy logic lives in a single place (utils) and
# parsing is *not* duplicated across modules.

from tusco_selector.utils import (
    GeneTranscripts,
    TranscriptExons,
    TranscriptStrand,
    load_assembly_mapping,
    parse_gtf,
)

# Public API -----------------------------------------------------------------

# The public signature is consumed by ``tusco_selector.cli``.  We purposefully
# keep the interface minimal: It only needs the *list of downloaded resource
# files* (as returned by ``downloader.fetch_all``) and the *species* code so we
# can locate helper files (like the RefSeq assembly report for human).

TranscriptInfo = Dict[str, object]  # kept unchecked on purpose


# ---------------------------------------------------------------------------
# Helper – build exon structure records
# ---------------------------------------------------------------------------

ExonStructure = Tuple[Tuple[Tuple[str, int, int], ...], str]  # (exons_sorted, strand)

# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------


def _compare_single_isoform_exon_structures(
    gtf_files: List[Path], chr_mapping: Dict[str, str]
):
    """Return *gene IDs* that have identical exon structures across *all* *gtf_files*.

    For each GTF, we first identify genes that have only one transcript.
    Then we compare structures across GTFs to find genes with consistent structures.
    """

    if not gtf_files:
        # empty – nothing to do; return placeholder empty structures
        return {}, [], {}, {}, {}

    # Step 1: Parse all GTFs and identify single-isoform genes + their structures in each
    all_gtf_structures = []
    structure_to_genes_map = []

    logger.info(
        "Processing %d GTF files to identify single-isoform genes", len(gtf_files)
    )
    for idx, gtf in enumerate(
        tqdm(gtf_files, desc="Processing GTF files", disable=not HAS_TQDM)
    ):
        start_time = time.time()
        logger.info("Parsing GTF file: %s", gtf.name)
        gene_tx, tx_exons, tx_strand = parse_gtf(gtf, chr_mapping, show_progress=True)
        parse_time = time.time() - start_time
        logger.info("Parsed %s in %.2f seconds (%d genes, %d transcripts)", 
                   gtf.name, parse_time, len(gene_tx), len(tx_exons))

        # Extract structures for single-isoform genes in this GTF
        gene_structures = {}
        # Create reverse mapping: structure -> set of genes with that structure
        structure_to_genes = defaultdict(set)

        # Process genes without nested progress bar for better performance
        gene_items = list(gene_tx.items())
        analysis_start = time.time()
        logger.info("Analyzing %d genes in %s", len(gene_items), gtf.name)
        for gene_id, tx_ids in gene_items:
            if len(tx_ids) != 1:
                continue

            tx_id = next(iter(tx_ids))
            # Sort exons by start coordinate only
            exons = sorted(
                tx_exons[tx_id], key=lambda x: (x[0], x[1])
            )  # Sort by chr, start
            strand = tx_strand.get(tx_id, "+")
            structure = (tuple(exons), strand)

            gene_structures[gene_id] = structure
            structure_to_genes[structure].add(gene_id)

        analysis_time = time.time() - analysis_start
        logger.info(
            "Found %d single-transcript genes in %s (analyzed in %.2f seconds)", 
            len(gene_structures), gtf.name, analysis_time
        )
        all_gtf_structures.append(gene_structures)
        structure_to_genes_map.append(structure_to_genes)

        if idx == 0:
            gene_tx_ref = gene_tx
            tx_exons_ref = tx_exons
            tx_strand_ref = tx_strand

    if not all_gtf_structures:
        return {}, [], {}, {}, {}

    # Get reference structures and genes
    gene_structures_ref = all_gtf_structures[0]

    # Find structures that exist in all GTFs
    matched_genes = []

    logger.info("Identifying genes with matching structures across all GTFs")
    if len(all_gtf_structures) == 1:
        # If only one GTF, all single-isoform genes match by definition
        matched_genes = list(gene_structures_ref.keys())
        logger.info(
            "Using all %d single-isoform genes from single GTF", len(matched_genes)
        )
    else:
        # Find structures that exist in all GTFs
        gene_items = list(gene_structures_ref.items())
        logger.info("Comparing %d gene structures across all GTFs", len(gene_items))
        for gene_id, ref_structure in gene_items:
            # Check if this structure exists in all other GTFs
            if all(
                ref_structure in structure_to_genes
                for structure_to_genes in structure_to_genes_map[1:]
            ):
                matched_genes.append(gene_id)

    logger.info(
        "Found %d genes with matching structures across all %d GTF files",
        len(matched_genes),
        len(gtf_files),
    )

    return (
        gene_structures_ref,
        matched_genes,
        gene_tx_ref,
        tx_exons_ref,
        tx_strand_ref,
    )


# ---------------------------------------------------------------------------
# Helper – filter mRNA transcripts
# ---------------------------------------------------------------------------


def filter_genes_present_in_mapping(
    transcripts: Dict[str, TranscriptInfo], mapping_file: str
) -> Dict[str, TranscriptInfo]:
    """Filter transcripts to keep only genes present in the mapping file.
    
    Note: This function checks if gene IDs exist in the first column of the mapping file.
    It does NOT specifically verify mRNA evidence (refseq_mrna column).

    Args:
        transcripts: Dictionary mapping gene IDs to transcript information
        mapping_file: Path to the mapping file containing gene IDs

    Returns:
        Dictionary containing only genes that exist in the mapping file
    """
    if not os.path.exists(mapping_file):
        logger.warning(
            "Mapping file not found at %s - skipping gene filtering",
            mapping_file,
        )
        return transcripts

    # Filter out genes that are not in the mapping file
    transcript_ids = set(transcripts.keys())
    mrna_transcript_ids = set()

    # Read the mapping file and check first column
    with open(mapping_file, "r") as f:
        for line in f:
            # Skip header if present
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if cols and cols[0] in transcript_ids:
                mrna_transcript_ids.add(cols[0])

    # Keep only genes that exist in the mapping file
    filtered_transcripts = {
        gene_id: info
        for gene_id, info in transcripts.items()
        if gene_id in mrna_transcript_ids
    }
    logger.info(
        f"Filtered out {len(transcript_ids) - len(mrna_transcript_ids)} genes that are not mRNA transcripts"
    )
    logger.info(f"Retained {len(mrna_transcript_ids)} mRNA transcripts")

    return filtered_transcripts


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def select_matched_single_isoforms(
    resource_files: List[str], *, species: str
) -> Dict[str, TranscriptInfo]:  # noqa: D401
    """Return a mapping *gene_id → transcript information*.

    The function automatically identifies GTF files among *resource_files*,
    performs the single-isoform cross-reference procedure and returns basic
    metadata for the matched genes.
    """
    start_time = time.time()
    logger.info("Starting single isoform selection process")

    gtf_files: List[Path] = [
        Path(f) for f in resource_files if re.search(r"\.gtf(?:\.gz)?$", str(f))
    ]

    if len(gtf_files) < 1:
        logger.warning("No GTF files provided – skipping single isoform selection")
        return {}

    logger.info("Found %d GTF files for single isoform selection", len(gtf_files))
    for gtf in gtf_files:
        logger.info("Using GTF file: %s", gtf)

    # Load assembly mapping for human (hsa). For other species this is either
    # not required or not provided.
    assembly_mapping: Dict[str, str] = {}

    # Heuristic: look for first '*.assembly_report.txt' among resources
    report_candidates = [
        Path(f) for f in resource_files if f.endswith("assembly_report.txt")
    ]
    assembly_report = report_candidates[0] if report_candidates else None
    if assembly_report:
        logger.info("Loading chromosome mapping from %s", assembly_report)
        assembly_mapping = load_assembly_mapping(assembly_report)
        logger.info(
            "Loaded %d chromosome mappings for species %s",
            len(assembly_mapping),
            species,
        )
    else:
        logger.warning("No assembly report found for human (hsa)")

    logger.info("Comparing single isoform exon structures across GTF files")
    (
        gene_structures_ref,
        matched_genes,
        gene_tx_ref,
        tx_exons_ref,
        tx_strand_ref,
    ) = _compare_single_isoform_exon_structures(gtf_files, assembly_mapping)

    # Build the return mapping using transcript/exon information from the
    # reference GTF (index 0).
    logger.info("Building transcript information for matched genes")
    transcripts: Dict[str, TranscriptInfo] = {}

    logger.info("Building transcript information for %d matched genes", len(matched_genes))
    for gene_id in matched_genes:
        tx_id = next(iter(gene_tx_ref[gene_id]))
        transcripts[gene_id] = {
            "transcript_id": tx_id,
            # Sort exons by start coordinate only when building final output
            "exons": sorted(
                tx_exons_ref[tx_id], key=lambda x: (x[0], x[1])
            ),  # Sort by chr, start
            "strand": tx_strand_ref[tx_id],
            # Store chromosome (assumes all exons on same chr)
            "chrom": tx_exons_ref[tx_id][0][0],
        }

    # Replace versioned gene IDs with cleaned IDs (without version numbers)
    cleaned_transcripts: Dict[str, TranscriptInfo] = {}
    replaced_count = 0
    for gene_id, info in transcripts.items():
        gene_id_clean = gene_id.split(".")[0]
        cleaned_transcripts[gene_id_clean] = info
        if gene_id_clean != gene_id:
            replaced_count += 1

    logger.info(
        "Replaced %d gene IDs with their version-stripped equivalents", replaced_count
    )
    
    total_time = time.time() - start_time
    logger.info(
        "Single isoform selection complete: selected %d genes in %.2f seconds", 
        len(cleaned_transcripts), total_time
    )
    return cleaned_transcripts
