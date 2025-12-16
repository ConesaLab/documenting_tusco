"""Pipeline runner module for orchestrating TUSCO selector steps."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Set, Tuple, TYPE_CHECKING

from tusco_selector.logging_utils import get_logger
from tusco_selector.utils import (
    subset_gtf_by_gene_ids,
    subset_mapping_file,
)

if TYPE_CHECKING:
    import pandas as pd

logger = get_logger(__name__)

# Total number of effective pipeline steps
TOTAL_STEPS = 6  # Step 6 integrated into Step 5


def get_step_paths(output_dir: str, step_num: int, slug: str) -> Tuple[str, str]:
    """Generate standard GTF and mapping file paths for a pipeline step.
    
    Args:
        output_dir: Output directory path
        step_num: Step number (1-based)
        slug: Step identifier slug
        
    Returns:
        Tuple of (gtf_path, mapping_path)
    """
    gtf_filename = f"step{step_num}_{slug}.gtf.gz"
    mapping_filename = f"step{step_num}_{slug}_mapping.tsv"
    return (
        os.path.join(output_dir, gtf_filename),
        os.path.join(output_dir, mapping_filename)
    )


def save_step_outputs(
    source_gtf_path: str,
    gene_ids_to_keep: Set[str],
    output_gtf_path: str,
    output_mapping_path: str,
    main_mapping_file: str | None,
    step_description: str,
) -> None:
    """Save GTF subset and corresponding mapping subset for a pipeline step.
    
    Args:
        source_gtf_path: Path to source GTF file
        gene_ids_to_keep: Set of gene IDs to retain
        output_gtf_path: Path for output GTF file
        output_mapping_path: Path for output mapping file
        main_mapping_file: Path to main mapping file (optional)
        step_description: Description of the step for logging
    """
    # Save subset GTF
    subset_gtf_by_gene_ids(source_gtf_path, gene_ids_to_keep, output_gtf_path)
    
    # Save subset mapping file if the main mapping exists
    if main_mapping_file and os.path.exists(main_mapping_file):
        subset_mapping_file(
            list(gene_ids_to_keep),
            main_mapping_file,
            output_mapping_path,
            step_description,
        )
    elif not main_mapping_file:
        logger.warning(
            f"Main mapping file path not provided for step: {step_description}. "
            f"Cannot create subset mapping."
        )


def print_step(step: int, total: int, description: str) -> None:
    """Print step header with progress indication.
    
    Args:
        step: Current step number
        total: Total number of steps
        description: Step description
    """
    logger.info(f"\n{'='*60}")
    logger.info(f"Step {step}/{total}: {description}")
    logger.info(f"{'='*60}")


def print_summary(gtf_path: str) -> None:
    """Print summary statistics for a GTF file.
    
    Args:
        gtf_path: Path to GTF file
    """
    if not os.path.exists(gtf_path):
        logger.warning(f"GTF file not found: {gtf_path}")
        return
    
    from tusco_selector.utils import count_exons_in_gtf
    
    exon_counts = count_exons_in_gtf(gtf_path)
    total_genes = len(exon_counts)
    single_exon = sum(1 for count in exon_counts.values() if count == 1)
    multi_exon = total_genes - single_exon
    
    logger.info(f"Genes: {total_genes:,} | Single-exon: {single_exon:,} | Multi-exon: {multi_exon:,}")


def run_pipeline(
    species: str,
    output_dir: str,
    cache_dir: str,
    thresholds: dict,
    args: object,
) -> Tuple[str, str, pd.DataFrame]:
    """Run the complete TUSCO selector pipeline.
    
    This is the main orchestrator that runs all pipeline steps in sequence.
    
    Args:
        species: Species code (hsa, mmu, dre)
        output_dir: Output directory path
        cache_dir: Cache directory for downloaded resources
        thresholds: Dictionary of filtering thresholds
        args: Command-line arguments object with additional parameters
        
    Returns:
        Tuple of (final_gtf_path, final_mapping_path, filter_log_df)
    """
    # Import pipeline steps (lazy imports to avoid circular dependencies)
    from tusco_selector.pipeline.resource_downloader import download_resources
    from tusco_selector.pipeline.single_isoform_selector import select_single_isoform_genes
    from tusco_selector.pipeline.splice_tss_check import check_splice_and_tss
    from tusco_selector.pipeline.manual_filter import apply_manual_filters
    from tusco_selector.pipeline.expression_filter import filter_by_expression
    from tusco_selector.pipeline.introverse_filter import filter_by_introverse
    from tusco_selector.cli import (
        initialize_filter_log,
        save_filter_log,
        summarize_results,
    )
    
    logger.info(f"Starting TUSCO selector pipeline for {species}")
    logger.info(f"Output directory: {output_dir}")
    
    # Step 1: Download resources
    print_step(1, TOTAL_STEPS, "Downloading and preparing resources")
    files, template_gtf = download_resources(species, cache_dir, args)
    print_summary(template_gtf)
    
    # Step 2: Select single-isoform genes
    print_step(2, TOTAL_STEPS, "Selecting matched single-isoform genes across annotations")
    transcripts, step2_gtf = select_single_isoform_genes(
        files, species, output_dir, template_gtf, cache_dir
    )
    print_summary(step2_gtf)
    
    # Initialize filter log
    filter_log_df = initialize_filter_log(output_dir, transcripts)
    
    # Step 3: Check splice junctions and TSS
    print_step(3, TOTAL_STEPS, "Checking splice junctions and transcription start sites")
    transcripts, step3_gtf, filter_log_df = check_splice_and_tss(
        transcripts,
        cache_dir,
        species,
        output_dir,
        filter_log_df,
        thresholds,
        args,
    )
    print_summary(step3_gtf)
    
    # Step 4: Apply manual filters
    if hasattr(args, 'skip_manual_filter') and not args.skip_manual_filter:
        print_step(4, TOTAL_STEPS, "Applying manual gene filters")
        transcripts, step4_gtf, filter_log_df = apply_manual_filters(
            transcripts,
            species,
            output_dir,
            step3_gtf,
            filter_log_df,
        )
        print_summary(step4_gtf)
        current_gtf = step4_gtf
    else:
        logger.info("Skipping manual filter step")
        current_gtf = step3_gtf
    
    # Step 5: Expression filtering (includes AlphaGenome)
    if hasattr(args, 'skip_expression_filter') and not args.skip_expression_filter:
        print_step(5, TOTAL_STEPS, "Filtering by expression data (Bgee/GTEx and AlphaGenome)")
        transcripts, step5_gtf, filter_log_df = filter_by_expression(
            transcripts,
            species,
            cache_dir,
            output_dir,
            current_gtf,
            filter_log_df,
            thresholds,
            args,
        )
        print_summary(step5_gtf)
        current_gtf = step5_gtf
    else:
        logger.info("Skipping expression filter step")
    
    # Step 6: IntroVerse filtering
    if hasattr(args, 'skip_introverse_filter') and not args.skip_introverse_filter:
        print_step(6, TOTAL_STEPS, "Filtering by IntroVerse mis-splicing data")
        transcripts, final_gtf, filter_log_df = filter_by_introverse(
            transcripts,
            species,
            cache_dir,
            output_dir,
            current_gtf,
            filter_log_df,
            thresholds,
        )
        print_summary(final_gtf)
    else:
        logger.info("Skipping IntroVerse filter step")
        final_gtf = current_gtf
    
    # Save filter log
    save_filter_log(filter_log_df, output_dir)
    
    # Generate final summary
    final_mapping = os.path.join(output_dir, "step7_final_filtered_mapping.tsv")
    summarize_results(species, output_dir, final_gtf, filter_log_df)
    
    return final_gtf, final_mapping, filter_log_df


__all__ = [
    'get_step_paths',
    'save_step_outputs',
    'print_step',
    'print_summary',
    'run_pipeline',
    'TOTAL_STEPS',
]