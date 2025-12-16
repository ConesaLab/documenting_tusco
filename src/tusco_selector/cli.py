from __future__ import annotations

# NOTE: Import *minimal* standard-library modules at import-time to keep
# dependency footprint light and avoid import errors when optional runtime
# dependencies (like pandas) are still missing during test discovery.

import argparse
import importlib
import logging
import os
from pathlib import Path
from typing import Dict, List, Set, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

# Shared utilities
from tusco_selector.utils import (
    build_transcripts_dict_from_gtf,
    count_exons_in_gtf,
    load_junctions_bed,
    subset_gtf_by_gene_ids,
    subset_mapping_file,
    generate_tusco_tables,
)

# Define custom logging levels
TRACE = 5
VERBOSE = 15

# Add custom levels to logging module
logging.addLevelName(TRACE, "TRACE")
logging.addLevelName(VERBOSE, "VERBOSE")

# Add convenience methods for custom levels
def _log_trace(self, message, *args, **kwargs):
    if self.isEnabledFor(TRACE):
        self._log(TRACE, message, args, **kwargs)

def _log_verbose(self, message, *args, **kwargs):
    if self.isEnabledFor(VERBOSE):
        self._log(VERBOSE, message, args, **kwargs)

logging.Logger.trace = _log_trace
logging.Logger.verbose = _log_verbose

# Configure basic logging (will be overridden by setup_logging)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Setup module logger
logger = logging.getLogger(__name__)

# Optional dependencies
# Removed colorized console dependencies; rely on logger for console output

# ---------------------------------------------------------------------------
# Global constants and small shared helpers
# ---------------------------------------------------------------------------

# Total number of pipeline steps (used for consistent progress output)
# Note: Step 6 (AlphaGenome) is now integrated into Step 5 (Expression filtering)
TOTAL_STEPS: int = 6

# (Removed verbose step-specific debug configuration to simplify CLI)

def setup_logging(log_level: str = "INFO", file_log_level: str = "VERBOSE", log_file: str | None = None) -> None:
    """Setup enhanced logging configuration with custom levels.
    
    Parameters
    ----------
    log_level : str
        Console logging level (TRACE, DEBUG, VERBOSE, INFO, WARNING, ERROR)
    file_log_level : str
        File logging level (TRACE, DEBUG, VERBOSE, INFO, WARNING, ERROR)
    log_file : str, optional
        Path to log file. If None, no file logging.
    """
    # Clear existing handlers
    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    
    # Console handler with cleaner format
    console_handler = logging.StreamHandler()
    console_level = getattr(logging, log_level.upper(), None)
    if console_level is None:
        console_level = {'TRACE': TRACE, 'VERBOSE': VERBOSE}.get(log_level.upper(), logging.INFO)
    console_handler.setLevel(console_level)
    
    # Simple format for console
    console_format = logging.Formatter(
        "%(message)s" if console_level >= logging.INFO else "%(levelname)s - %(message)s"
    )
    console_handler.setFormatter(console_format)
    root_logger.addHandler(console_handler)
    
    # File handler with detailed format
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w')
        file_level = getattr(logging, file_log_level.upper(), None)
        if file_level is None:
            file_level = {'TRACE': TRACE, 'VERBOSE': VERBOSE}.get(file_log_level.upper(), VERBOSE)
        file_handler.setLevel(file_level)
        
        file_format = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        file_handler.setFormatter(file_format)
        root_logger.addHandler(file_handler)
    
    # Set root logger to lowest level to allow handlers to filter
    root_logger.setLevel(TRACE)


# Built-in mapping from common tissue names to UBERON IDs per species.
# Extend this as additional tissues are validated.
TISSUE_TO_UBERON: Dict[str, Dict[str, str]] = {
    "mmu": {
        "kidney": "UBERON:0002113",
        # Add more validated mappings here as needed
    },
    "hsa": {
        # Example: "kidney": "UBERON:0002113",
    },
}

def _resolve_uberon_id(species: str, tissue_name: str) -> str | None:
    """Map a user-specified tissue name to an UBERON ID for the given species.

    Tissue names are compared case-insensitively. Returns None if not found.
    """
    sp = (species or "").lower()
    t = (tissue_name or "").strip().lower()
    return TISSUE_TO_UBERON.get(sp, {}).get(t)


def _project_root() -> str:
    """Return the absolute project root directory.

    This resolves to the repository root by walking up two levels from this
    module's file (src/tusco_selector/...). Centralizing this avoids
    duplicating long ``os.path`` expressions across helpers.
    """
    # src/tusco_selector/cli.py -> src -> project root
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _mapping_file_path(cache_dir: str) -> str:
    """Return the canonical path to the Ensembl↔RefSeq mapping TSV.

    Centralizing this avoids repeating string joins across steps, keeping the
    code concise and consistent.
    """
    return os.path.join(cache_dir, "annotation", "mapping", "ensembl_refseq.tsv")


def _subset_mapping_if_missing(transcripts, mapping_file: str, dest_mapping: str, description: str) -> None:
    """Create a subset mapping if not present, based on current transcripts.

    This helper reduces boilerplate across step functions.
    """
    from tusco_selector.utils import subset_mapping_file

    if not os.path.exists(dest_mapping) and os.path.exists(mapping_file):
        subset_mapping_file(
            list(transcripts.keys()), mapping_file, dest_mapping, description
        )


def _copy_file_if_missing(src_path: str, dest_path: str, reason: str) -> None:
    """Copy a file from src to dest if dest doesn't exist and src does.

    Parameters
    ----------
    src_path
        Source file path.
    dest_path
        Destination file path.
    reason
        Short log message to clarify why the copy is performed.
    """
    if os.path.exists(src_path) and not os.path.exists(dest_path):
        import shutil

        logger.info(reason)
        shutil.copy(src_path, dest_path)

"""Minimal logger-driven console helpers."""


def print_header(title: str) -> None:
    logger.info(title)


def print_step(step_num: int, total_steps: int, description: str) -> None:
    logger.info(f"[{step_num}/{total_steps}] {description}")


def print_success(message: str) -> None:
    logger.info(message)


def print_warning(message: str) -> None:
    logger.warning(message)


def print_info(message: str) -> None:
    logger.info(message)


def print_summary(gtf_path: str | os.PathLike[str]) -> None:
    """Count and log a brief gene/exon summary for a GTF file.

    Logs a warning if the file does not exist or an error occurs while
    computing the summary. This function never raises so that it is safe for
    pipeline progress logging.
    """
    if not os.path.exists(gtf_path):
        logger.warning(f"GTF file not found, cannot generate summary: {gtf_path}")
        return
    try:
        total_genes, single_exon, multi_exon = count_exons_in_gtf(gtf_path)
        if total_genes > 0:
            logger.info(
                f"  Summary: {total_genes} genes total ({single_exon} single-exon, {multi_exon} multi-exon)"
            )
        else:
            logger.info(f"  Summary: 0 genes found in GTF.")
    except Exception as e:
        logger.warning(f"Could not generate summary for {gtf_path}: {e}")
        logger.error(f"Failed to count exons in {gtf_path}: {e}")


def log_step_summary(step_name: str, genes_in: int, genes_out: int, reasons: dict = None, time_elapsed: float = None) -> None:
    """Log a clean summary for a pipeline step.
    
    Parameters
    ----------
    step_name : str
        Name of the step
    genes_in : int
        Number of genes entering the step
    genes_out : int
        Number of genes after the step
    reasons : dict, optional
        Dictionary of removal reasons and counts
    time_elapsed : float, optional
        Time taken for the step in seconds
    """
    logger.info("=" * 60)
    logger.info(f"Step Summary: {step_name}")
    logger.info(f"  Genes in:  {genes_in:,}")
    logger.info(f"  Genes out: {genes_out:,}")
    logger.info(f"  Removed:   {genes_in - genes_out:,}")
    
    if reasons:
        logger.info("  Removal reasons:")
        for reason, count in sorted(reasons.items(), key=lambda x: x[1], reverse=True):
            logger.info(f"    - {reason}: {count}")
    
    if time_elapsed:
        if time_elapsed < 60:
            logger.info(f"  Time:      {time_elapsed:.1f}s")
        else:
            logger.info(f"  Time:      {time_elapsed/60:.1f}m")
    
    logger.info("=" * 60)


# (BioMart helpers removed as part of deprecating BioMart-based validation.)


# ---------------------------------------------------------------------------
# Public helper functions required by the test-suite
# ---------------------------------------------------------------------------


def list_species() -> List[str]:
    """Return the supported three-letter species codes."""
    return [
        "hsa",  # Human (Homo sapiens)
        "mmu",  # Mouse (Mus musculus)
        "dre",  # Zebrafish (Danio rerio)
    ]


def run(species: str, output_dir: str | None = None, **threshold_kwargs) -> None:
    """Programmatic entry-point mirroring the CLI *main* behaviour.

    Parameters
    ----------
    species : str
        Three-letter species code (hsa, mmu, dre)
    output_dir : str, optional
        Directory to store output files
    **threshold_kwargs
        Threshold parameters for filtering steps. Available options:

        Splice/TSS filtering:
        - novel_threshold : float (default: 0.01)
        - min_novel_length : int (default: 80)

        GTEx expression filtering (human tissue-specific):
        - gtex_prevalence_expression_cutoff : float (default: 0.1)
        - gtex_median_expression_cutoff : float (default: 1.0)
        - gtex_prevalence_threshold : float (default: 0.95)

        Bgee expression filtering (mouse and other species):
        - min_genes_per_tissue : int (default: 25000)
        - bgee_quality : list (default: ["gold quality", "silver quality"])
        - bgee_prevalence_threshold : float (default: 0.95)

        Alphagenome filtering - Universal genes:
        - single_exon_median_rpkm_threshold : float (default: 1.0)
        - multi_exon_median_rpkm_threshold : float (default: 1.0)
        - single_exon_expression_threshold : float (default: 0.1)  # RPKM
        - multi_exon_expression_threshold : float (default: 0.1)    # RPKM
        - single_exon_prevalence_threshold : float (default: 0.95)
        - multi_exon_prevalence_threshold : float (default: 0.95)
        - splice_ratio_threshold : float (default: 0.05)
        - splice_purity_threshold : float (default: 0.95)

        Alphagenome filtering - Tissue-specific genes (stricter):
        - tissue_expression_threshold : float (default: 1.0)
        - tissue_splice_ratio_threshold : float (default: 0.01)
    """
    argv = [species]
    if output_dir is not None:
        argv.extend(["--output-dir", output_dir])

    # Add threshold arguments
    for key, value in threshold_kwargs.items():
        arg_name = f"--{key.replace('_', '-')}"
        if isinstance(value, list):
            argv.extend([arg_name] + value)
        else:
            argv.extend([arg_name, str(value)])

    main(argv)


def _load_external_resources(species: str):
    """Dynamically import the *resources* module for *species*."""
    return importlib.import_module(
        f"tusco_selector.species.{species}.resources", package=None
    )


# ---------------------------------------------------------------------------
# Helper Functions within CLI
# ---------------------------------------------------------------------------
def _update_filter_log(
    filter_log_df: "pd.DataFrame" | None,
    removed_items: Dict[str, str] | Set[str] | List[str],
    reason_prefix: str,
) -> "pd.DataFrame" | None:
    """
    Updates the filter log DataFrame with reasons for removed genes.

    Args:
        filter_log_df: The DataFrame tracking gene filtering reasons.
        removed_items: Either a dictionary {gene_id: specific_reason} or a set/list of gene_ids.
        reason_prefix: The general reason for filtering in this step.

    Returns:
        The updated DataFrame or None if the input df was None.
    """
    if filter_log_df is None:
        return None

    if isinstance(removed_items, dict):
        # Handle dictionary input {gene_id: specific_reason}
        for gene_id, specific_reason in removed_items.items():
            if gene_id in filter_log_df.index:
                if filter_log_df.loc[gene_id, "FilterReason"] == "Passed":
                    filter_log_df.loc[gene_id, "FilterReason"] = (
                        f"{reason_prefix} ({specific_reason})"
                    )
            else:
                logger.warning(
                    f"Gene {gene_id} from {reason_prefix} filter not found in initial filter log."
                )
    elif isinstance(removed_items, (set, list)):
        # Handle set or list input (gene_ids)
        for gene_id in removed_items:
            if gene_id in filter_log_df.index:
                if filter_log_df.loc[gene_id, "FilterReason"] == "Passed":
                    filter_log_df.loc[gene_id, "FilterReason"] = reason_prefix
            else:
                logger.warning(
                    f"Gene {gene_id} from {reason_prefix} filter not found in initial filter log."
                )
    else:
        logger.error(
            f"Unsupported type for removed_items in _update_filter_log: {type(removed_items)}"
        )

    return filter_log_df


def _get_step_paths(output_dir: str, step_num: int, slug: str) -> Tuple[str, str]:
    """Generates the standard GTF and mapping file paths for a pipeline step."""
    gtf_filename = f"step{step_num}_{slug}.gtf.gz"
    mapping_filename = f"step{step_num}_{slug}_mapping.tsv"
    return os.path.join(output_dir, gtf_filename), os.path.join(
        output_dir, mapping_filename
    )


def _save_step_outputs(
    source_gtf_path: str,
    gene_ids_to_keep: Set[str],
    output_gtf_path: str,
    output_mapping_path: str,
    main_mapping_file: str | None,
    step_description: str,
) -> None:
    """Handles saving the GTF subset and the corresponding mapping subset for a step."""
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
            f"Main mapping file path not provided for step: {step_description}. Cannot create subset mapping."
        )
    # No warning if main_mapping_file was provided but doesn't exist, subset_mapping_file handles that internally (logs info).


# ---------------------------------------------------------------------------
# Helper functions for expression filtering
# ---------------------------------------------------------------------------


def _process_bgee_expression_data_unified(cache_dir, species, bgee_thresholds):
    """Process Bgee expression data for both human and mouse."""
    from tusco_selector.pipeline.expression_filter import (
        get_universal_high_expression_genes_bgee_with_ids,
    )

    universal_genes_bgee, high_exp_per_tissue_bgee = set(), {}  # Initialize
    bgee_expr_data_file_path = None
    bgee_mapping_tsv_path = None

    # Find expression data directory
    expr_dir = os.path.join(cache_dir, "expression")
    if not os.path.isdir(expr_dir):
        expr_dir = cache_dir  # Fallback to cache_dir

    try:
        resources_module = _load_external_resources(species)
        bgee_resources = (
            getattr(resources_module, "EXTERNAL_RESOURCES", {})
            .get("expression", {})
            .get("bgee", {})
        )

        if bgee_resources:
            bgee_expr_info = bgee_resources.get("bgee_expression", {})
            bgee_mapping_info = bgee_resources.get("bgee_mapping", {})

            # Define base path for Bgee files, assuming they are in a 'bgee' subdirectory within expr_dir
            bgee_subdir_path = os.path.join(expr_dir, "bgee")

            if not os.path.isdir(bgee_subdir_path):
                logger.warning(
                    f"Bgee subdirectory not found at {bgee_subdir_path}. Bgee files might not be located correctly."
                )

            # Locate Bgee expression data file (.tsv.gz)
            if bgee_expr_info and "output" in bgee_expr_info:
                expr_file_name = bgee_expr_info["output"]
                potential_expr_path = os.path.join(bgee_subdir_path, expr_file_name)
                if os.path.exists(potential_expr_path) and os.path.isfile(
                    potential_expr_path
                ):
                    bgee_expr_data_file_path = potential_expr_path
                    print_info(
                        f"Found Bgee expression data file: {os.path.basename(bgee_expr_data_file_path)}"
                    )
                else:
                    logger.warning(
                        f"Bgee expression data file not found at expected path: {potential_expr_path}"
                    )
            else:
                logger.warning("Bgee expression resource 'output' info missing.")

            # Locate Bgee mapping file - different structure for human vs mouse
            if bgee_mapping_info and "output" in bgee_mapping_info:
                mapping_archive_name = bgee_mapping_info[
                    "output"
                ]  # e.g., "bgee.rna_seq_experiments_libraries.tar.gz"
                extracted_folder_name = mapping_archive_name.replace(".tar.gz", "")

                # Determine the correct species-specific mapping file name
                if species == "hsa":
                    mapping_file_name = "Homo_sapiens_RNA-Seq_libraries.tsv"
                elif species == "mmu":
                    mapping_file_name = "Mus_musculus_RNA-Seq_libraries.tsv"
                else:
                    mapping_file_name = None
                    logger.warning(f"Unsupported species for Bgee mapping: {species}")

                if mapping_file_name:
                    # First check if the mapping file exists directly in the bgee directory
                    direct_mapping_path = os.path.join(
                        bgee_subdir_path, mapping_file_name
                    )

                    if os.path.exists(direct_mapping_path) and os.path.isfile(
                        direct_mapping_path
                    ):
                        bgee_mapping_tsv_path = direct_mapping_path
                        print_info(
                            f"Found Bgee mapping file: {os.path.basename(bgee_mapping_tsv_path)}"
                        )
                    else:
                        # If not found directly, try to extract from archive
                        archive_path = os.path.join(
                            bgee_subdir_path, mapping_archive_name
                        )
                        extracted_folder_path = os.path.join(
                            bgee_subdir_path, extracted_folder_name
                        )
                        potential_mapping_path = os.path.join(
                            extracted_folder_path, mapping_file_name
                        )

                        if os.path.exists(archive_path) and not os.path.exists(
                            extracted_folder_path
                        ):
                            print_info(
                                f"Extracting Bgee mapping archive: {os.path.basename(archive_path)}"
                            )
                            try:
                                import tarfile

                                with tarfile.open(archive_path, "r:gz") as tar:
                                    tar.extractall(path=bgee_subdir_path)
                                print_success(
                                    f"Extracted Bgee mapping archive to: {extracted_folder_path}"
                                )
                            except Exception as e:
                                logger.error(
                                    f"Failed to extract Bgee mapping archive: {e}"
                                )
                                return set(), {}

                        if os.path.exists(potential_mapping_path) and os.path.isfile(
                            potential_mapping_path
                        ):
                            bgee_mapping_tsv_path = potential_mapping_path
                            print_info(
                                f"Found Bgee mapping file: {os.path.basename(bgee_mapping_tsv_path)}"
                            )
                        else:
                            logger.warning(
                                f"Bgee mapping file not found at expected paths:"
                            )
                            logger.warning(f"  Direct path: {direct_mapping_path}")
                            logger.warning(
                                f"  Extracted path: {potential_mapping_path}"
                            )
                            if os.path.exists(extracted_folder_path):
                                logger.warning(
                                    f"Extracted folder exists but mapping file not found. Contents: {os.listdir(extracted_folder_path)}"
                                )
                            if os.path.exists(bgee_subdir_path):
                                logger.warning(
                                    f"Bgee directory contents: {os.listdir(bgee_subdir_path)}"
                                )
            else:
                logger.warning("Bgee mapping resource 'output' info missing.")
        else:
            logger.warning(
                "Bgee resource configuration not found in species resources."
            )

    except Exception as e:
        logger.error(f"Error locating Bgee files: {e}")

    if not bgee_expr_data_file_path:
        print_warning(
            "Bgee expression data file (*_expr_simple.tsv.gz) could not be located. Skipping expression filtering."
        )
        return set(), {}
    elif not bgee_mapping_tsv_path:
        print_warning(
            f"Bgee mapping file ({species}_RNA-Seq_libraries.tsv) could not be located. Skipping expression filtering."
        )
        return set(), {}

    # Process Bgee data if both files are found
    print_info("Processing Bgee expression data...")

    # Check for GTEx mapping file (for human only)
    gtex_mapping_file_path = None
    if species == "hsa":
        gtex_mapping_file = os.path.join(cache_dir, "expression", "gtex", "uberon-map-gtex-subgroup.tsv")
        if os.path.exists(gtex_mapping_file):
            gtex_mapping_file_path = gtex_mapping_file
            logger.info(f"Found GTEx tissue mapping file: {gtex_mapping_file_path}")
        else:
            logger.warning(f"GTEx tissue mapping file not found: {gtex_mapping_file}")

    bgee_processing_kwargs = {
        "bgee_mapping_file_path": bgee_mapping_tsv_path,
        "min_genes_per_tissue": bgee_thresholds["min_genes_per_tissue"],
        "qual_ok": bgee_thresholds["qual_ok"],
        "prevalence_threshold": bgee_thresholds["prevalence_threshold"],
        "species": species,
        "gtex_mapping_file_path": gtex_mapping_file_path,
    }

    try:
        universal_genes_bgee, high_exp_per_tissue_bgee = (
            get_universal_high_expression_genes_bgee_with_ids(
                bgee_expr_data_file_path, **bgee_processing_kwargs
            )
        )
        print_success(
            f"Bgee processing complete: {len(universal_genes_bgee)} universal genes found"
        )
    except Exception as e:
        logger.error(f"Failed processing Bgee file: {e}")
        return set(), {}

    return universal_genes_bgee, high_exp_per_tissue_bgee


def _process_hybrid_expression_data_human(cache_dir, species, expression_thresholds):
    """Process both Bgee (for universal genes) and GTEx (for tissue-specific genes) for human.
    
    Returns:
        Tuple of (universal_genes_from_bgee, tissue_specific_genes_from_gtex)
    """
    print_info("Processing hybrid expression data for human...")
    print_info("  • Using Bgee for universal gene identification")
    print_info("  • Using GTEx for tissue-specific gene identification")
    
    # Get universal genes from Bgee
    universal_genes_bgee, _ = _process_bgee_expression_data_unified(
        cache_dir, species, expression_thresholds["bgee"]
    )
    
    # Get tissue-specific genes from GTEx (ignore universal genes from GTEx)
    _, high_exp_per_tissue_gtex = _process_gtex_expression_data_unified(
        cache_dir, species, expression_thresholds["gtex"]
    )
    
    print_info(f"Hybrid processing complete:")
    print_info(f"  • Bgee universal genes: {len(universal_genes_bgee)}")
    print_info(f"  • GTEx tissue-specific sets: {len(high_exp_per_tissue_gtex)} tissues")
    
    return universal_genes_bgee, high_exp_per_tissue_gtex


def _process_gtex_expression_data_unified(cache_dir, species, gtex_thresholds):
    """Process GTEx (.gct.gz) expression data for human.

    Returns a tuple of (universal_gene_ids, high_exp_per_tissue) where the
    second element maps (UBERON/EFO anatomical ID, tissue name) to the set of
    gene IDs passing per-tissue thresholds. Behaviour is unchanged; this only
    centralizes a few path utilities and adds clearer logging.
    """
    from tusco_selector.pipeline.expression_filter import (
        get_universal_high_expression_genes,
    )
    import json
    from pathlib import Path
    
    universal_genes_gtex, high_exp_per_tissue_gtex = set(), {}
    
    # Find expression data directory
    expr_dir = os.path.join(cache_dir, "expression", "gtex")
    if not os.path.isdir(expr_dir):
        print_warning(f"GTEx expression directory not found: {expr_dir}")
        return set(), {}
    
    # Load ontology mapping to convert tissue names to anatomical IDs
    ontology_mapping_path = os.path.join(
        _project_root(), "data", "common", "ontology_term_to_tissue_mapping.json"
    )
    
    if not os.path.exists(ontology_mapping_path):
        print_warning(f"Ontology mapping file not found: {ontology_mapping_path}")
        return set(), {}
    
    try:
        with open(ontology_mapping_path, "r") as f:
            ontology_mapping = json.load(f)
        
        # We'll use the ontology mapping directly (ontology_id -> tissue_name)
        # since GTEx filenames are already UBERON/EFO IDs
    except Exception as e:
        print_warning(f"Failed to load ontology mapping: {e}")
        return set(), {}
    
    # Find all GTEx .gct.gz files
    gct_files = []
    try:
        for file_path in Path(expr_dir).glob("*.gct.gz"):
            gct_files.append(file_path)
        
        if not gct_files:
            print_warning(f"No GTEx .gct.gz files found in {expr_dir}")
            return set(), {}
        
        print_info(f"Found {len(gct_files)} GTEx expression files")
        
        # Process GTEx files using the provided function
        universal_genes_gtex, high_exp_per_tissue_raw = get_universal_high_expression_genes(
            gct_files,
            prevalence_expression_cutoff=gtex_thresholds.get("prevalence_expression_cutoff", 0.1),
            median_expression_cutoff=gtex_thresholds.get("median_expression_cutoff", 1.0),
            prevalence_threshold=gtex_thresholds.get("prevalence_threshold", 0.95),
            universal_tissue_fraction=1.0,  # Set to 1.0 for semantic correctness (universal genes discarded in hybrid mode)
        )
        
        print_success(f"GTEx processing complete: {len(universal_genes_gtex)} universal genes found")
        
        # Convert tissue names to anatomical IDs for compatibility with downstream processing
        high_exp_per_tissue_gtex = {}
        for tissue_identifier, gene_set in high_exp_per_tissue_raw.items():
            # The tissue_identifier is already the UBERON/EFO ID from the filename
            anatomical_id = tissue_identifier
            
            # Look up the human-readable tissue name from the ontology mapping
            tissue_name = None
            if "human" in ontology_mapping and anatomical_id in ontology_mapping["human"]:
                tissue_name = ontology_mapping["human"][anatomical_id]
                print_info(f"Mapped {anatomical_id} -> {tissue_name} ({len(gene_set)} genes)")
            else:
                # Fallback: use the anatomical ID as the tissue name
                tissue_name = anatomical_id
                print_info(f"Using anatomical ID as name: {anatomical_id} ({len(gene_set)} genes)")
            
            # Use tuple format (anatomical_id, tissue_name) as expected downstream
            tissue_info = (anatomical_id, tissue_name)
            high_exp_per_tissue_gtex[tissue_info] = gene_set
        
    except Exception as e:
        logger.error(f"Failed processing GTEx files: {e}")
        return set(), {}
    
    return universal_genes_gtex, high_exp_per_tissue_gtex


def _load_alphagenome_data(cache_dir, species):
    """Load and process Alphagenome data for gene filtering."""
    import json

    # Load ontology mapping
    ontology_mapping_path = os.path.join(
        _project_root(), "data", "common", "ontology_term_to_tissue_mapping.json"
    )

    if not os.path.exists(ontology_mapping_path):
        logger.warning(f"Ontology mapping file not found: {ontology_mapping_path}")
        return None, None

    try:
        with open(ontology_mapping_path, "r") as f:
            ontology_mapping = json.load(f)
    except Exception as e:
        logger.error(f"Failed to load ontology mapping: {e}")
        return None, None

    # Get species-specific mapping
    species_name = (
        "human" if species == "hsa" else "mouse" if species == "mmu" else None
    )
    if species_name not in ontology_mapping:
        logger.warning(f"Species {species_name} not found in ontology mapping")
        return None, None

    tissue_mapping = ontology_mapping[species_name]

    # Load Alphagenome data
    alphagenome_file = f"alphagenome_{species_name}.json"
    alphagenome_path = os.path.join(cache_dir, alphagenome_file)

    if not os.path.exists(alphagenome_path):
        logger.warning(f"Alphagenome file not found: {alphagenome_path}")
        return None, None

    try:
        with open(alphagenome_path, "r") as f:
            alphagenome_data = json.load(f)
    except Exception as e:
        logger.error(f"Failed to load Alphagenome data: {e}")
        return None, None

    print_success(f"Loaded Alphagenome data for {len(alphagenome_data)} genes")
    return alphagenome_data, tissue_mapping


def _process_alphagenome_gene_data(gene_data):
    """Process a single gene's Alphagenome data."""

    def _clean_tissue_key(key):
        """Remove trailing labels after first space."""
        return key.split()[0] if " " in key else key

    # Process tissue expression data
    tissue_expression = {}
    for key, value in gene_data.get("tissue_expression", {}).items():
        clean_key = _clean_tissue_key(key)
        # If we already have this tissue, take the maximum value (any >= 0.1 TPM counts as expressed)
        if clean_key in tissue_expression:
            tissue_expression[clean_key] = max(tissue_expression[clean_key], value)
        else:
            tissue_expression[clean_key] = value

    # Process tissue splice ratios
    tissue_splice_ratios = {}
    for key, value in gene_data.get("tissue_splice_ratios", {}).items():
        # Remove "junction_" prefix and trailing labels
        if key.startswith("junction_"):
            clean_key = _clean_tissue_key(key[9:])  # Remove "junction_" prefix
            # If we already have this tissue, take the minimum ratio (most conservative)
            if clean_key in tissue_splice_ratios:
                tissue_splice_ratios[clean_key] = min(
                    tissue_splice_ratios[clean_key], value
                )
            else:
                tissue_splice_ratios[clean_key] = value

    return tissue_expression, tissue_splice_ratios


def _analyze_independent_alphagenome_filters(
    genes_to_filter,
    alphagenome_data,
    tissue_mapping,
    alphagenome_thresholds,
    single_exon_genes: Set[str] | None = None,
    multi_exon_genes: Set[str] | None = None,
):
    """Analyze how many genes each Alphagenome filter would remove if applied independently.

    This provides insight into the true restrictiveness of each filter criterion,
    as opposed to the sequential filtering used in _apply_alphagenome_filter.

    Returns:
        dict: Independent filter results with keys 'median', 'prevalence', 'splice'
              Each value is the number of genes that would be filtered out by that criterion alone.
    """

    if not alphagenome_data or not tissue_mapping:
        return {"median": 0, "prevalence": 0, "splice": 0}

    # Create a lookup dict for faster access
    alphagenome_lookup = {gene["gene_id"]: gene for gene in alphagenome_data}

    # Get all available tissue ontology terms from the mapping
    available_tissues = set(tissue_mapping.keys())

    # Track genes that would fail each filter independently
    median_failures = 0
    prevalence_failures = 0
    splice_failures = 0

    genes_with_data = 0  # Count genes that have alphagenome data

    for gene_id in genes_to_filter:
        # Strip version from gene ID if present for lookup
        lookup_gene_id = gene_id.split(".")[0]

        if lookup_gene_id not in alphagenome_lookup:
            continue  # Skip genes not in alphagenome data for this analysis

        gene_data = alphagenome_lookup[lookup_gene_id]
        tissue_expression, tissue_splice_ratios = _process_alphagenome_gene_data(
            gene_data
        )

        # Filter to only tissues that we have mapping for
        filtered_expression = {
            k: v for k, v in tissue_expression.items() if k in available_tissues
        }
        filtered_splice_ratios = {
            k: v for k, v in tissue_splice_ratios.items() if k in available_tissues
        }

        if not filtered_expression:
            continue  # Skip genes without expression data for mapped tissues

        genes_with_data += 1

        # Test each filter independently

        # 1. Median RPKM threshold (class-specific)
        expression_values = list(filtered_expression.values())
        median_rpkm = (
            sorted(expression_values)[len(expression_values) // 2]
            if expression_values
            else 0
        )
        # Determine class-specific thresholds with conservative fallback (min of both) if unknown
        base_id = gene_id.split(".")[0]
        is_single = single_exon_genes is not None and base_id in single_exon_genes
        is_multi = multi_exon_genes is not None and base_id in multi_exon_genes
        med_thr_single = float(alphagenome_thresholds.get("single_exon_median_rpkm_threshold", 1.0))
        med_thr_multi = float(alphagenome_thresholds.get("multi_exon_median_rpkm_threshold", 1.0))
        expr_thr_single = float(alphagenome_thresholds.get("single_exon_expression_threshold", 0.1))
        expr_thr_multi = float(alphagenome_thresholds.get("multi_exon_expression_threshold", 0.1))
        prev_thr_single = float(alphagenome_thresholds.get("single_exon_prevalence_threshold", 0.95))
        prev_thr_multi = float(alphagenome_thresholds.get("multi_exon_prevalence_threshold", 0.95))

        median_threshold = med_thr_single if is_single else med_thr_multi if is_multi else min(med_thr_single, med_thr_multi)
        expression_threshold = expr_thr_single if is_single else expr_thr_multi if is_multi else min(expr_thr_single, expr_thr_multi)
        prevalence_threshold = prev_thr_single if is_single else prev_thr_multi if is_multi else min(prev_thr_single, prev_thr_multi)

        if median_rpkm < median_threshold:
            median_failures += 1

        # 2. Prevalence: fraction of tissues with RPKM > expression_threshold
        expressed_tissues = sum(1 for rpkm in expression_values if rpkm > expression_threshold)
        prevalence = (
            expressed_tissues / len(expression_values) if expression_values else 0
        )

        if prevalence < prevalence_threshold:
            prevalence_failures += 1

        # 3. Splice purity: fraction of tissues with splice ratio < splice_ratio_threshold
        splice_ratios = []
        for tissue_id in filtered_expression.keys():
            ratio = filtered_splice_ratios.get(
                tissue_id, 0.0
            )  # Default to 0 if missing
            splice_ratios.append(ratio)

        pure_splice_tissues = sum(
            1
            for ratio in splice_ratios
            if ratio < alphagenome_thresholds["splice_ratio_threshold"]
        )
        splice_purity = pure_splice_tissues / len(splice_ratios) if splice_ratios else 0

        if splice_purity < alphagenome_thresholds["splice_purity_threshold"]:
            splice_failures += 1

    return {
        "median": median_failures,
        "prevalence": prevalence_failures,
        "splice": splice_failures,
        "total_analyzed": genes_with_data,
    }


def _apply_alphagenome_filter(
    genes_to_filter,
    alphagenome_data,
    tissue_mapping,
    alphagenome_thresholds,
    single_exon_genes: Set[str] | None = None,
    multi_exon_genes: Set[str] | None = None,
):
    """Apply Alphagenome filtering criteria to a set of genes."""

    if not alphagenome_data or not tissue_mapping:
        logger.warning(
            "Alphagenome data or tissue mapping not available, skipping filter"
        )
        return genes_to_filter, {}

    # Create a lookup dict for faster access
    alphagenome_lookup = {gene["gene_id"]: gene for gene in alphagenome_data}

    passing_genes = set()
    removed_genes = {}

    # Get all available tissue ontology terms from the mapping
    available_tissues = set(tissue_mapping.keys())

    for gene_id in genes_to_filter:
        # Strip version from gene ID if present for lookup
        lookup_gene_id = gene_id.split(".")[0]

        if lookup_gene_id not in alphagenome_lookup:
            removed_genes[gene_id] = "not in Alphagenome data"
            continue

        gene_data = alphagenome_lookup[lookup_gene_id]
        tissue_expression, tissue_splice_ratios = _process_alphagenome_gene_data(
            gene_data
        )

        # Filter to only tissues that we have mapping for
        filtered_expression = {
            k: v for k, v in tissue_expression.items() if k in available_tissues
        }
        filtered_splice_ratios = {
            k: v for k, v in tissue_splice_ratios.items() if k in available_tissues
        }

        if not filtered_expression:
            removed_genes[gene_id] = "no expression data for mapped tissues"
            continue

        # Apply filtering criteria

        # 1. Median RPKM threshold (class-specific)
        expression_values = list(filtered_expression.values())
        median_rpkm = (
            sorted(expression_values)[len(expression_values) // 2]
            if expression_values
            else 0
        )
        # Determine class-specific thresholds with conservative fallback (min of both) if unknown
        base_id = gene_id.split(".")[0]
        is_single = single_exon_genes is not None and base_id in single_exon_genes
        is_multi = multi_exon_genes is not None and base_id in multi_exon_genes
        med_thr_single = float(alphagenome_thresholds.get("single_exon_median_rpkm_threshold", 1.0))
        med_thr_multi = float(alphagenome_thresholds.get("multi_exon_median_rpkm_threshold", 1.0))
        expr_thr_single = float(alphagenome_thresholds.get("single_exon_expression_threshold", 0.1))
        expr_thr_multi = float(alphagenome_thresholds.get("multi_exon_expression_threshold", 0.1))
        prev_thr_single = float(alphagenome_thresholds.get("single_exon_prevalence_threshold", 0.95))
        prev_thr_multi = float(alphagenome_thresholds.get("multi_exon_prevalence_threshold", 0.95))
        median_threshold = med_thr_single if is_single else med_thr_multi if is_multi else min(med_thr_single, med_thr_multi)
        expression_threshold = expr_thr_single if is_single else expr_thr_multi if is_multi else min(expr_thr_single, expr_thr_multi)
        prevalence_threshold = prev_thr_single if is_single else prev_thr_multi if is_multi else min(prev_thr_single, prev_thr_multi)

        if median_rpkm < median_threshold:
            removed_genes[gene_id] = (
                f"median RPKM {median_rpkm:.3f} < {median_threshold}"
            )
            continue

        # 2. Prevalence: fraction of tissues with RPKM > expression_threshold
        expressed_tissues = sum(1 for rpkm in expression_values if rpkm > expression_threshold)
        prevalence = (
            expressed_tissues / len(expression_values) if expression_values else 0
        )

        if prevalence < prevalence_threshold:
            removed_genes[gene_id] = (
                f"prevalence {prevalence:.3f} < {prevalence_threshold}"
            )
            continue

        # 3. Splice purity: fraction of tissues with splice ratio < splice_ratio_threshold
        # For tissues without splice data, assume ratio = 0 (fully annotated)
        splice_ratios = []
        for tissue_id in filtered_expression.keys():
            ratio = filtered_splice_ratios.get(
                tissue_id, 0.0
            )  # Default to 0 if missing
            splice_ratios.append(ratio)

        pure_splice_tissues = sum(
            1
            for ratio in splice_ratios
            if ratio < alphagenome_thresholds["splice_ratio_threshold"]
        )
        splice_purity = pure_splice_tissues / len(splice_ratios) if splice_ratios else 0

        if splice_purity < alphagenome_thresholds["splice_purity_threshold"]:
            removed_genes[gene_id] = (
                f"splice purity {splice_purity:.3f} < {alphagenome_thresholds['splice_purity_threshold']}"
            )
            continue

        # Gene passes all filters
        passing_genes.add(gene_id)

    return passing_genes, removed_genes


def _determine_tissue_data_sources(tissues, alphagenome_data, tissue_mapping, expression_source="bgee"):
    """Determine which data sources are available for each tissue.
    
    Args:
        tissues: List of (tissue_id, tissue_name) tuples
        alphagenome_data: Alphagenome data if available
        tissue_mapping: Tissue mapping data if available  
        expression_source: "bgee" or "gtex" indicating the source of tissue data
    
    Returns:
        dict: {tissue_id: {'sources': ['bgee'|'gtex'|'alphagenome'], 'name': tissue_name}}
    """
    tissue_data_sources = {}
    
    # Determine primary expression source based on parameter
    primary_source = 'gtex' if expression_source == 'gtex' else 'bgee'
    
    # Get tissues from the expression source (tissue_id, tissue_name tuples)
    for tissue_id, tissue_name in tissues:
        tissue_data_sources[tissue_id] = {
            'sources': [primary_source],
            'name': tissue_name
        }
    
    # Check which expression source tissues also have Alphagenome data
    if alphagenome_data and tissue_mapping:
        # Get all tissues that have data in Alphagenome
        alphagenome_tissue_ids = set()
        for gene_data in alphagenome_data:
            tissue_expression, _ = _process_alphagenome_gene_data(gene_data)
            alphagenome_tissue_ids.update(tissue_expression.keys())
        
        # Add Alphagenome as a source for tissues that have it
        for tissue_id in tissue_data_sources:
            if tissue_id in alphagenome_tissue_ids:
                tissue_data_sources[tissue_id]['sources'].append('alphagenome')
        
        # Add tissues that only exist in Alphagenome
        for tissue_id in alphagenome_tissue_ids:
            if tissue_id not in tissue_data_sources:
                # Get tissue name from mapping if available, ensure it's a string
                tissue_name = tissue_mapping.get(tissue_id, tissue_id)
                if not isinstance(tissue_name, str):
                    tissue_name = str(tissue_name) if tissue_name is not None else str(tissue_id)
                tissue_data_sources[tissue_id] = {
                    'sources': ['alphagenome'],
                    'name': tissue_name
                }
    
    return tissue_data_sources


def _apply_tissue_specific_filtering(
    tissue_genes,
    tissue_id,
    tissue_name,
    data_sources,
    expression_genes,
    alphagenome_data,
    tissue_mapping,
    alphagenome_thresholds,
):
    """Apply appropriate filtering based on available data sources for this tissue.
    
    Logic:
    - Bgee/GTEx only: Use expression source genes (already filtered)
    - Alphagenome only: Apply Alphagenome filtering
    - Both expression source and Alphagenome: Gene must pass BOTH sources (intersection)
    """
    
    # Log initial status
    logger.info(f"\n=== Tissue: {tissue_name} ({tissue_id}) ===")
    logger.info(f"  Input genes (from previous filters): {len(tissue_genes)}")
    logger.info(f"  Data sources available: {', '.join(data_sources)}")
    
    if ('bgee' in data_sources or 'gtex' in data_sources) and 'alphagenome' not in data_sources:
        # Expression source only (Bgee or GTEx) - use expression filtered genes
        expr_present = tissue_genes & expression_genes
        result = expr_present
        
        source_name = 'GTEx' if 'gtex' in data_sources else 'Bgee'
        logger.info(f"  {source_name} filtering:")
        logger.info(f"    - Genes present in {source_name}: {len(expression_genes)}")
        logger.info(f"    - Overlap with input genes: {len(expr_present)}")
        logger.info(f"  Final genes for this tissue: {len(result)}")
        
        print_info(f"  • {tissue_name}: {len(result)}/{len(tissue_genes)} genes ({source_name} only)")
        return result
        
    elif 'alphagenome' in data_sources and 'bgee' not in data_sources and 'gtex' not in data_sources:
        # Alphagenome only - apply Alphagenome filtering
        alphagenome_passing = _apply_alphagenome_filter_only(
            tissue_genes, tissue_id, alphagenome_data, tissue_mapping, alphagenome_thresholds
        )
        result = alphagenome_passing
        
        logger.info(f"  Alphagenome filtering:")
        logger.info(f"    - Genes passing Alphagenome thresholds: {len(alphagenome_passing)}")
        logger.info(f"  Final genes for this tissue: {len(result)}")
        
        print_info(f"  • {tissue_name}: {len(result)}/{len(tissue_genes)} genes (Alphagenome only)")
        return result
        
    elif ('bgee' in data_sources or 'gtex' in data_sources) and 'alphagenome' in data_sources:
        # Both expression source and Alphagenome available - gene must pass BOTH sources (intersection)
        expr_present = tissue_genes & expression_genes
        alphagenome_passing = _apply_alphagenome_filter_only(
            tissue_genes, tissue_id, alphagenome_data, tissue_mapping, alphagenome_thresholds
        )
        result = expr_present & alphagenome_passing
        
        source_name = 'GTEx' if 'gtex' in data_sources else 'Bgee'
        logger.info(f"  {source_name} filtering:")
        logger.info(f"    - Genes present in {source_name}: {len(expression_genes)}")
        logger.info(f"    - Overlap with input genes: {len(expr_present)}")
        logger.info(f"  Alphagenome filtering:")
        logger.info(f"    - Genes passing Alphagenome thresholds: {len(alphagenome_passing)}")
        logger.info(f"  Combined filtering (intersection):")
        logger.info(f"    - Genes passing BOTH {source_name} AND Alphagenome: {len(result)}")
        logger.info(f"  Final genes for this tissue: {len(result)}")
        
        print_info(f"  • {tissue_name}: {len(result)}/{len(tissue_genes)} genes ({source_name}: {len(expr_present)}, Alpha: {len(alphagenome_passing)}, Both: {len(result)})")
        return result
    
    else:
        # No data sources available - shouldn't happen
        logger.warning(f"  No data sources available for tissue {tissue_name}")
        print_warning(f"  • {tissue_name}: No data sources available")
        return set()


def _apply_alphagenome_filter_only(tissue_genes, tissue_id, alphagenome_data, tissue_mapping, alphagenome_thresholds):
    """Apply only Alphagenome filtering to genes."""
    if not alphagenome_data:
        return tissue_genes
        
    alphagenome_lookup = {gene["gene_id"]: gene for gene in alphagenome_data}
    passing_genes = set()
    
    # Track statistics for logging
    genes_not_in_alphagenome = 0
    genes_without_tissue_data = 0
    genes_failed_tpm = 0
    genes_failed_splice = 0
    genes_passed = 0
    
    for gene_id in tissue_genes:
        lookup_gene_id = gene_id.split(".")[0]
        
        if lookup_gene_id not in alphagenome_lookup:
            # Gene not in Alphagenome - pass through
            passing_genes.add(gene_id)
            genes_not_in_alphagenome += 1
            continue
            
        gene_data = alphagenome_lookup[lookup_gene_id]
        tissue_expression, tissue_splice_ratios = _process_alphagenome_gene_data(gene_data)
        
        if tissue_id not in tissue_expression:
            # No tissue data - pass through  
            passing_genes.add(gene_id)
            genes_without_tissue_data += 1
            continue
            
        tpm = tissue_expression[tissue_id]
        splice_ratio = tissue_splice_ratios.get(tissue_id, 0.0)
        
        if (
            tpm > alphagenome_thresholds["expression_threshold"]
            and splice_ratio < alphagenome_thresholds["splice_ratio_threshold"]
        ):
            passing_genes.add(gene_id)
            genes_passed += 1
        else:
            # Track why genes failed
            if tpm <= alphagenome_thresholds["expression_threshold"]:
                genes_failed_tpm += 1
            if splice_ratio >= alphagenome_thresholds["splice_ratio_threshold"]:
                genes_failed_splice += 1
    
    # Log detailed statistics
    logger.debug(f"    Alphagenome filter details:")
    logger.debug(f"      - Genes not in Alphagenome DB (passed through): {genes_not_in_alphagenome}")
    logger.debug(f"      - Genes without tissue data (passed through): {genes_without_tissue_data}")
    logger.debug(f"      - Genes passing TPM > {alphagenome_thresholds['expression_threshold']} and splice < {alphagenome_thresholds['splice_ratio_threshold']}: {genes_passed}")
    logger.debug(f"      - Genes failed TPM threshold: {genes_failed_tpm}")
    logger.debug(f"      - Genes failed splice ratio threshold: {genes_failed_splice}")
    
    return passing_genes


def _apply_alphagenome_tissue_filter(
    tissue_genes,
    tissue_id,
    tissue_name,
    alphagenome_data,
    tissue_mapping,
    alphagenome_thresholds,
):
    """Apply Alphagenome filtering to genes in a specific tissue."""

    if not alphagenome_data or not tissue_mapping:
        logger.info(f"No Alphagenome data available for tissue {tissue_name} ({tissue_id}) - passing through all {len(tissue_genes)} genes")
        return tissue_genes

    # Create a lookup dict for faster access
    alphagenome_lookup = {gene["gene_id"]: gene for gene in alphagenome_data}

    passing_genes = set()
    genes_without_alphagenome_data = 0
    genes_without_tissue_expression = 0
    genes_filtered_out = 0
    genes_failed_tpm = 0
    genes_failed_splice = 0

    for gene_id in tissue_genes:
        # Strip version from gene ID if present for lookup
        lookup_gene_id = gene_id.split(".")[0]

        if lookup_gene_id not in alphagenome_lookup:
            # Pass through genes not in Alphagenome data without filtering
            genes_without_alphagenome_data += 1
            passing_genes.add(gene_id)
            continue

        gene_data = alphagenome_lookup[lookup_gene_id]
        tissue_expression, tissue_splice_ratios = _process_alphagenome_gene_data(
            gene_data
        )

        # Check if this specific tissue has data
        if tissue_id not in tissue_expression:
            # Pass through genes without tissue expression data - no filtering applied
            genes_without_tissue_expression += 1
            passing_genes.add(gene_id)
            continue

        # For tissue-specific filtering, only check expression and splicing in this specific tissue
        tpm = tissue_expression[tissue_id]
        splice_ratio = tissue_splice_ratios.get(
            tissue_id, 0.0
        )  # Default to 0 if missing

        # Apply tissue-specific criteria and track reasons
        if tpm <= alphagenome_thresholds["expression_threshold"]:
            genes_filtered_out += 1
            genes_failed_tpm += 1
        elif splice_ratio >= alphagenome_thresholds["splice_ratio_threshold"]:
            genes_filtered_out += 1
            genes_failed_splice += 1
        else:
            passing_genes.add(gene_id)

    # Log the filtering results with detailed breakdown
    total_genes = len(tissue_genes)
    genes_evaluated = total_genes - genes_without_alphagenome_data - genes_without_tissue_expression
    
    logger.verbose(f"Alphagenome tissue filter for {tissue_name} ({tissue_id}):")
    logger.verbose(f"  Total genes: {total_genes}")
    logger.verbose(f"  Genes passing: {len(passing_genes)} ({100*len(passing_genes)/total_genes:.1f}%)")
    logger.verbose(f"  Genes without Alphagenome data (passed): {genes_without_alphagenome_data}")
    logger.verbose(f"  Genes without tissue expression (passed): {genes_without_tissue_expression}")
    if genes_filtered_out > 0:
        logger.verbose(f"  Genes filtered out: {genes_filtered_out}")
        if genes_failed_tpm > 0:
            logger.verbose(f"    - Below TPM threshold (≤{alphagenome_thresholds['expression_threshold']}): {genes_failed_tpm}")
        if genes_failed_splice > 0:
            logger.verbose(f"    - Above splice ratio threshold (≥{alphagenome_thresholds['splice_ratio_threshold']}): {genes_failed_splice}")

    return passing_genes


def _load_ontology_tissue_mapping(species: str) -> Dict[str, str]:
    """Load tissue name mapping from ontology_term_to_tissue_mapping.json and GTEx mapping for human."""
    project_root = _project_root()
    ontology_file = os.path.join(project_root, "data", "common", "ontology_term_to_tissue_mapping.json")
    
    try:
        import json
        with open(ontology_file, 'r') as f:
            ontology_data = json.load(f)
        
        if species == "hsa":
            species_key = "human"
        elif species == "mmu":
            species_key = "mouse"
        else:
            species_key = species
            
        mapping = ontology_data.get(species_key, {})
        
        # For human, also load GTEx mapping file to handle UBERON/EFO codes
        if species == "hsa":
            gtex_mapping_file = os.path.join(project_root, "data", "hsa", "expression", "gtex", "uberon-map-gtex-subgroup.tsv")
            if os.path.exists(gtex_mapping_file):
                import pandas as pd
                df = pd.read_csv(gtex_mapping_file, sep='\t')
                for _, row in df.iterrows():
                    tissue_name = row['uberon_description_gtex_subgroup']
                    # Add UBERON mapping with both formats
                    if pd.notna(row['uberon_code']):
                        uberon_code = row['uberon_code']
                        mapping[uberon_code] = tissue_name  # UBERON_xxxx
                        mapping[uberon_code.replace('_', ':')] = tissue_name  # UBERON:xxxx
                    # Add EFO mapping with both formats
                    if pd.notna(row['efo_code']):
                        efo_code = row['efo_code']
                        mapping[efo_code] = tissue_name  # EFO_xxxx
                        mapping[efo_code.replace('_', ':')] = tissue_name  # EFO:xxxx
        
        logger.verbose(f"Loaded {len(mapping)} tissue mappings for species {species}")
        return mapping
    except Exception as e:
        logger.error(f"Failed to load ontology tissue mapping: {e}")
        return {}


def _create_tissue_tsv(output_file: str, tissue_name: str, gene_ids: Set[str], mapping_file: str = None) -> None:
    """Create a TSV file with tissue name and corresponding TUSCO genes in the 6-column format."""
    import csv
    
    try:
        # Load mapping data if available
        mapping_data = {}
        if mapping_file and os.path.exists(mapping_file):
            with open(mapping_file, 'r') as mf:
                reader = csv.reader(mf, delimiter='\t')
                for row in reader:
                    if row and not row[0].startswith('#'):
                        # Mapping format: Gene_ID, Transcript_ID, Gene_Symbol, EntrezGene_ID, RefSeq_mRNA_ID, RefSeq_Protein_ID
                        if len(row) >= 6:
                            gene_id = row[0].split('.')[0]  # Remove version if present
                            mapping_data[gene_id] = row
        
        with open(output_file, 'w') as f:
            # Write header
            f.write(f"# Tissue: {tissue_name}\n")
            f.write(f"# Gene Count: {len(gene_ids)}\n")
            f.write("Gene_ID\tTranscript_ID\tGene_Symbol\tEntrezGene_ID\tRefSeq_mRNA_ID\tRefSeq_Protein_ID\n")
            
            # Write gene data with all 6 columns
            for gene_id in sorted(gene_ids):
                lookup_id = gene_id.split('.')[0]  # Remove version for lookup
                if lookup_id in mapping_data:
                    # Use the full row from mapping
                    row = mapping_data[lookup_id]
                    f.write('\t'.join(row[:6]) + '\n')  # Ensure we only write 6 columns
                else:
                    # If not in mapping, write gene_id with empty fields
                    f.write(f"{gene_id}\t\t\t\t\t\n")
                
        logger.verbose(f"Created tissue TSV file: {output_file} ({len(gene_ids)} genes)")
    except Exception as e:
        logger.error(f"Failed to create tissue TSV file {output_file}: {e}")


def _create_tissue_statistics_tsv(output_file: str, tissue_stats: List[Tuple[str, str, int]]) -> None:
    """Create a statistical TSV file with tissue names and corresponding TUSCO gene counts.
    
    Parameters
    ----------
    output_file : str
        Path to the output statistics TSV file
    tissue_stats : List[Tuple[str, str, int]]
        List of tuples containing (tissue_name, anatomical_id, gene_count)
    """
    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write("Anatomical_Entity_Name\tNumber_of_TUSCO_Genes\n")
            
            # Write tissue statistics (sorted by gene count, descending)
            for tissue_name, anatomical_id, gene_count in tissue_stats:
                f.write(f"{tissue_name}\t{gene_count}\n")
                
        logger.info(f"Created tissue statistics TSV file: {output_file} ({len(tissue_stats)} tissues)")
    except Exception as e:
        logger.error(f"Failed to create tissue statistics TSV file {output_file}: {e}")


def _create_comprehensive_statistics_tsv(output_dir: str, species: str) -> None:
    """Create comprehensive statistics TSV file including universal and tissue-specific gene counts.
    
    Parameters
    ----------
    output_dir : str
        Output directory containing the TUSCO files
    species : str  
        Species code (hsa, mmu, dre)
    """
    import glob
    import os
    
    try:
        # Get universal gene count from the main TUSCO file
        universal_file = None
        if species == "hsa":
            universal_file = os.path.join(output_dir, "tusco_human.tsv")
        elif species == "mmu": 
            universal_file = os.path.join(output_dir, "tusco_mouse.tsv")
        
        universal_count = 0
        universal_genes = {}
        
        if universal_file and os.path.exists(universal_file):
            with open(universal_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        universal_count += 1
                        gene_id = parts[0]
                        universal_genes[gene_id] = parts
        
        # Collect tissue-specific data
        data = [['Universal', universal_count]]
        
        # Map species codes to full names (same as in tissue file creation)
        species_name_map = {"hsa": "human", "mmu": "mouse", "dre": "zebrafish"}
        species_name = species_name_map.get(species, species)
        
        # Find all tissue-specific TSV files using the full species name (exclude statistics and detailed files)
        tissue_pattern = os.path.join(output_dir, f"tusco_{species_name}_*.tsv")
        tissue_files = glob.glob(tissue_pattern)
        tissue_files = [f for f in tissue_files if not f.endswith('tissue_statistics.tsv') and not f.endswith('comprehensive_statistics.tsv')]
        
        # Load ontology mapping for proper tissue name resolution
        ontology_mapping = _load_ontology_tissue_mapping(species)
        
        for tissue_file in sorted(tissue_files):
            # Extract tissue name from filename using the full species name
            raw_name = os.path.basename(tissue_file).replace(f'tusco_{species_name}_', '').replace('.tsv', '')
            
            # Check if it's a UBERON/EFO code and map it, otherwise clean the name
            if raw_name.startswith('UBERON:') or raw_name.startswith('EFO:'):
                # The mapping now handles both UBERON: and UBERON_ formats
                tissue_name = ontology_mapping.get(raw_name, raw_name)
            else:
                tissue_name = raw_name.replace('_', ' ')
            
            # Count genes in file (excluding header lines starting with #)
            count = 0
            with open(tissue_file, 'r') as f:
                for line in f:
                    if not line.strip().startswith('#') and line.strip() and not line.strip().startswith('Gene_ID'):
                        count += 1
            
            data.append([tissue_name, count])
        
        # Write comprehensive statistics file
        stats_file = os.path.join(output_dir, f"tusco_{species}_comprehensive_statistics.tsv")
        with open(stats_file, 'w') as f:
            f.write('Anatomical_Entity_Name\tNumber_of_TUSCO_Genes\n')
            for row in data:
                f.write(f'{row[0]}\t{row[1]}\n')
        
        logger.info(f"Created comprehensive statistics TSV file: {stats_file} ({len(data)} entities)")
        print_success(f"Created comprehensive statistics file: {stats_file}")
        
    except Exception as e:
        logger.error(f"Failed to create comprehensive statistics TSV file: {e}")
        print_warning(f"Failed to create comprehensive statistics: {e}")


def _save_expression_outputs_with_ids(
    passing_universal_genes,
    high_exp_per_tissue,
    initial_genes,
    step4_gtf,
    output_dir,
    mapping_file,
    species,
):
    """Save universal and tissue-specific expression outputs using anatomical entity IDs."""
    from tusco_selector.utils import subset_gtf_by_gene_ids, subset_mapping_file, generate_tusco_tables

    if passing_universal_genes:
        # Write universal gene output
        universal_gtf = os.path.join(
            output_dir, "step5_universal_high_expression.gtf.gz"
        )
        subset_gtf_by_gene_ids(step4_gtf, passing_universal_genes, universal_gtf)

        # Subset mapping file for universal genes
        if os.path.exists(mapping_file):
            universal_mapping = os.path.join(
                output_dir, "step5_universal_high_expression_mapping.tsv"
            )
            subset_mapping_file(
                list(passing_universal_genes),
                mapping_file,
                universal_mapping,
                "universal expression (present in all Bgee tissues)",
            )

        print_success(
            f"Universal genes: {len(passing_universal_genes)} genes expressed in all tissues"
        )

        # Load tissue name mapping for proper naming
        ontology_mapping = _load_ontology_tissue_mapping(species)
        
        # Save tissue-specific high-expression genes with proper naming
        from tusco_selector.optional_deps import tqdm, HAS_TQDM
        
        # Map species codes to full names for file naming
        species_name_map = {"hsa": "human", "mmu": "mouse", "dre": "zebrafish"}
        species_name = species_name_map.get(species, species)

        tissue_stats = []
        logger.info(f"Generating tissue-specific outputs for {len(high_exp_per_tissue)} tissues...")
        
        # Use progress bar if available and progress flag is set
        iterator = tqdm(high_exp_per_tissue.items(), desc="Tissues", disable=not HAS_TQDM) if HAS_TQDM else high_exp_per_tissue.items()
        
        for tissue_info, genes in iterator:
            # tissue_info is a tuple of (anatomical_entity_id, tissue_name)
            anatomical_id, tissue_name = tissue_info
            
            # Ensure both anatomical_id and tissue_name are strings
            anatomical_id = str(anatomical_id) if anatomical_id is not None else ""
            tissue_name = str(tissue_name) if tissue_name is not None else ""

            tissue_passing_genes = (
                genes & initial_genes
            )  # Filter against genes entering this step
            if not tissue_passing_genes:
                continue

            # Get proper tissue name from ontology mapping
            clean_tissue_name = ontology_mapping.get(anatomical_id, tissue_name)
            # Ensure tissue name is a string
            if not isinstance(clean_tissue_name, str):
                clean_tissue_name = str(clean_tissue_name) if clean_tissue_name is not None else str(tissue_name)
            # Sanitize filename - replace spaces and special characters
            safe_tissue_name = clean_tissue_name.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "").replace("-", "_")
            
            # Track stats for reporting
            tissue_stats.append((clean_tissue_name, anatomical_id, len(tissue_passing_genes)))

            # Use new naming scheme: tusco_{species_name}_{tissue_name}.gtf/.tsv
            tissue_gtf = os.path.join(output_dir, f"tusco_{species_name}_{safe_tissue_name}.gtf.gz")
            subset_gtf_by_gene_ids(step4_gtf, tissue_passing_genes, tissue_gtf)

            # Create TSV file with tissue name and genes
            tissue_tsv = os.path.join(output_dir, f"tusco_{species_name}_{safe_tissue_name}.tsv")
            _create_tissue_tsv(tissue_tsv, clean_tissue_name, tissue_passing_genes, mapping_file)

        # Create statistical summary TSV file
        if tissue_stats:
            tissue_stats.sort(key=lambda x: x[2], reverse=True)
            
            # Create statistical summary file
            stats_file = os.path.join(output_dir, f"tusco_{species_name}_tissue_statistics.tsv")
            _create_tissue_statistics_tsv(stats_file, tissue_stats)
            
            # Calculate statistics for summary
            gene_counts = [count for _, _, count in tissue_stats]
            total_tissues = len(tissue_stats)
            min_genes = min(gene_counts) if gene_counts else 0
            max_genes = max(gene_counts) if gene_counts else 0
            avg_genes = sum(gene_counts) / total_tissues if total_tissues > 0 else 0
            
            # Print compact summary
            logger.info(f"Generated {total_tissues} tissue-specific files:")
            logger.info(f"  • Gene range: {min_genes}-{max_genes} genes/tissue")
            logger.info(f"  • Average: {avg_genes:.1f} genes/tissue")
            
            # Show top tissues at VERBOSE level
            for tissue_name, anatomical_id, count in tissue_stats[:5]:
                logger.verbose(f"  • {tissue_name} ({anatomical_id}): {count} genes")
            if len(tissue_stats) > 5:
                logger.verbose(f"  ... and {len(tissue_stats) - 5} more tissues")
            
            logger.info(f"Created tissue statistics TSV file: {stats_file} ({total_tissues} tissues)")
    else:
        print_warning("No universal genes detected in Bgee data")


# ---------------------------------------------------------------------------
# ENCODE augmentation for tissue-specific sets (mouse kidney)
# ---------------------------------------------------------------------------

def _augment_tissue_specific_with_encode(
    cache_dir: str,
    species: str,
    high_exp_per_tissue: Dict[tuple, Set[str]],
    expression_thresholds: Dict,
    encode_tissues: List[str] | None,
) -> Dict[tuple, Set[str]]:
    """Augment tissue-specific gene sets with ENCODE TSVs for selected tissues.

    Parameters
    - cache_dir: project cache directory for the species (data/<species>)
    - species: three-letter code
    - high_exp_per_tissue: existing mapping {(UBERON, name) -> set(genes)} from Bgee/GTEx
    - expression_thresholds: thresholds dict, uses alphagenome.tissue_specific.expression_threshold
    - encode_tissues: list of tissue names provided by user; each will be mapped to UBERON via TISSUE_TO_UBERON
    """
    if not encode_tissues:
        return high_exp_per_tissue

    # Thresholds: use ENCODE-specific if provided, else fall back to Alphagenome tissue-specific
    tpm_prevalence_cutoff = expression_thresholds.get("encode", {}).get(
        "prevalence_expression_cutoff",
        expression_thresholds.get("alphagenome", {}).get("tissue_specific", {}).get("expression_threshold", 0.1),
    )
    tpm_median_cutoff = expression_thresholds.get("encode", {}).get(
        "median_expression_cutoff",
        expression_thresholds.get("alphagenome", {}).get("tissue_specific", {}).get("expression_threshold", 1.0),
    )
    prevalence_threshold = expression_thresholds.get("encode", {}).get(
        "prevalence_threshold", 0.95
    )

    # Load species resources to locate ENCODE outputs
    try:
        resources_module = _load_external_resources(species)
        encode_resources = (
            getattr(resources_module, "EXTERNAL_RESOURCES", {})
            .get("expression", {})
            .get("encode", {})
        )
    except Exception:
        encode_resources = {}

    for tissue in encode_tissues:
        uberon_id = _resolve_uberon_id(species, tissue)
        if not uberon_id:
            print_warning(f"Unknown tissue mapping for '{tissue}' in species '{species}'. Skipping.")
            continue

        # Determine encoded filename key in resources (by convention it is the tissue key)
        filename = None
        if isinstance(encode_resources, dict):
            info = encode_resources.get(tissue, {})
            filename = info.get("output") if isinstance(info, dict) else None
        if not filename:
            # Fallback to standard naming
            filename = f"encode.{tissue}.expression.tsv"

        encode_path = os.path.join(cache_dir, "expression", "encode", filename)
        if not os.path.exists(encode_path):
            print_warning(f"ENCODE file for '{tissue}' not found: {encode_path}")
            continue

        try:
            # Lazy import to keep pandas optional at import-time
            import pandas as pd

            print_info(
                (
                    f"Applying ENCODE augmentation for '{tissue}' ({uberon_id}) from {encode_path} "
                    f"with thresholds: median >= {tpm_median_cutoff}, prevalence >= {prevalence_threshold} @ TPM > {tpm_prevalence_cutoff}"
                )
            )

            df = pd.read_csv(encode_path, sep="\t", skiprows=1)
            if "Feature ID" not in df.columns or "TPM" not in df.columns:
                print_warning(
                    f"ENCODE file has unexpected columns, skipping: {encode_path}"
                )
                continue

            # Determine Ensembl gene prefix by species
            ens_prefix = "ENSG" if species == "hsa" else "ENSMUSG" if species == "mmu" else "ENS"

            # Compute per-gene prevalence and median using all ENCODE entries in the file
            df = df[df["Feature ID"].astype(str).str.startswith(ens_prefix)].copy()
            df["GeneID"] = df["Feature ID"].astype(str).str.split(".").str[0]
            df["TPM"] = df["TPM"].astype(float)

            grouped = df.groupby("GeneID")["TPM"]
            n_entries = grouped.count()
            prevalence = (df["TPM"] > float(tpm_prevalence_cutoff)).groupby(df["GeneID"]).mean()
            median_expr = grouped.median()

            passing = (prevalence >= float(prevalence_threshold)) & (median_expr >= float(tpm_median_cutoff))
            gene_ids = set(passing[passing].index)
            if not gene_ids:
                print_warning(
                    f"ENCODE {tissue} ({uberon_id}): 0/{df['GeneID'].nunique()} genes passed thresholds"
                )
                continue

            tissue_key = (uberon_id, tissue.lower())
            prev = high_exp_per_tissue.get(tissue_key, set())
            if not prev:
                # If there is no Bgee/GTEx set for this tissue, intersection yields empty
                intersect_genes = set()
            else:
                intersect_genes = prev & gene_ids
            high_exp_per_tissue[tissue_key] = intersect_genes

            considered = int(df["GeneID"].nunique())
            sample_e = ", ".join(sorted(list(gene_ids))[:5])
            if prev:
                sample_prev = ", ".join(sorted(list(prev))[:5])
                sample_int = ", ".join(sorted(list(intersect_genes))[:5])
                print_info(
                    f"ENCODE {tissue} ({uberon_id}): ENCODE-pass {len(gene_ids)}/{considered}; "
                    f"Bgee/GTEx tissue set {len(prev)}; intersection {len(intersect_genes)}; sample INT: {sample_int}"
                )
            else:
                print_info(
                    f"ENCODE {tissue} ({uberon_id}): ENCODE-pass {len(gene_ids)}/{considered}; "
                    f"no existing Bgee/GTEx set found for intersection"
                )
        except Exception as e:
            logger.warning(f"Failed to augment {tissue} with ENCODE data: {e}")

    return high_exp_per_tissue

# ---------------------------------------------------------------------------
# Pipeline execution
# ---------------------------------------------------------------------------


def run_pipeline_steps(args, cache_dir, files, template_gtf, thresholds):
    """Execute all pipeline steps and return final results."""
    from tusco_selector.utils import (
        build_transcripts_dict_from_gtf,
        subset_gtf_by_gene_ids,
        subset_mapping_file,
    )

    # Use module-level TOTAL_STEPS for consistent progress output

    # Step 2: Select single-isoform genes
    print_step(
        2, TOTAL_STEPS, "Selecting matched single-isoform genes across annotations"
    )
    transcripts, step2_gtf = select_single_isoform_genes(
        files, args.species, args.output_dir, template_gtf, cache_dir
    )
    print_summary(step2_gtf)

    # Initialize filter log DataFrame
    filter_log_df = initialize_filter_log(args.output_dir, transcripts)

    # Step 3: Check splice junctions and TSS
    print_step(
        3, TOTAL_STEPS, "Checking splice junctions and transcription start sites"
    )
    transcripts, step3_gtf, filter_log_df = check_splice_and_tss(
        transcripts,
        cache_dir,
        args.species,
        args.output_dir,
        step2_gtf,
        filter_log_df,
        thresholds["splice_tss"],
    )
    print_summary(step3_gtf)

    # Step 4: Apply IntroVerse filter (human only)
    if args.species == "hsa":
        print_step(4, TOTAL_STEPS, "Filtering novel splice junctions (IntroVerse)")
        transcripts, step4_gtf, filter_log_df = apply_introverse_filter(
            transcripts,
            args.species,
            cache_dir,
            args.output_dir,
            step3_gtf,
            filter_log_df,
        )
        print_summary(step4_gtf)
    else:
        print_warning("IntroVerse filter not available for this species")
        step4_gtf = step3_gtf

    # Step 5: Filter by expression level
    print_step(5, TOTAL_STEPS, "Filtering by expression level")
    transcripts, step5_gtf, filter_log_df = filter_by_expression(
        transcripts,
        args.species,
        cache_dir,
        args.output_dir,
        step4_gtf,
        filter_log_df,
        thresholds["expression"],
        expression_source=args.expression_source,
        encode_tissues=args.encode_tissues,
    )
    print_summary(step5_gtf)

    # Step 6: Alphagenome filtering is now integrated into Step 5
    # Universal genes require both Bgee AND Alphagenome criteria
    # Tissue-specific genes use per-tissue logic based on data source availability
    print_step(6, TOTAL_STEPS, "Alphagenome filtering (integrated into Step 5)")
    print_info("Alphagenome filtering has been integrated into expression filtering (Step 5)")
    print_info("- Universal genes: Must pass BOTH Bgee AND Alphagenome criteria")  
    print_info("- Tissue-specific genes: Per-tissue logic (Bgee only, Alphagenome only, or BOTH when available)")
    step6_gtf = step5_gtf  # Use the same GTF as Step 5 since filtering is integrated
    print_summary(step6_gtf)

    # Step 7: Apply manual filters
    print_step(
        7, TOTAL_STEPS, "Applying manual filters and generating annotation files"
    )
    transcripts, step7_gtf, filter_log_df = apply_manual_filters(
        transcripts, cache_dir, args.output_dir, step6_gtf, args.species, filter_log_df
    )
    print_summary(step7_gtf)

    # Write filter log
    if filter_log_df is not None:
        write_filter_log(filter_log_df, args.output_dir)

    return transcripts


def initialize_filter_log(output_dir, transcripts):
    """Initialize the filter log DataFrame from the step 2 results."""
    # Lazy import to keep pandas optional at import-time
    import pandas as pd
    step2_mapping_path = os.path.join(output_dir, "step2_single_isoform_mapping.tsv")
    filter_log_df = None
    initial_gene_ids = set(transcripts.keys())

    if os.path.exists(step2_mapping_path):
        try:
            # Read the mapping file, skipping comments and using the first column (ENSG IDs)
            mapping_df = pd.read_csv(
                step2_mapping_path,
                sep="\t",
                comment="#",
                header=None,
                usecols=[0],
                engine="python",
            )
            mapping_df.columns = ["GeneID_Version"]
            # Strip version from Gene IDs if present
            mapping_df["GeneID"] = mapping_df["GeneID_Version"].apply(
                lambda x: x.split(".")[0]
            )

            # Use the gene IDs from the initial transcripts dict as the base
            filter_log_df = pd.DataFrame({"GeneID": list(initial_gene_ids)})
            filter_log_df["FilterReason"] = "Passed"
            filter_log_df.set_index("GeneID", inplace=True)
            logger.info(
                f"Initialized filter log with {len(filter_log_df)} genes from step 2."
            )

        except Exception as e:
            logger.error(f"Error reading step 2 mapping file for filter log: {e}")
            # Fallback: Use gene IDs directly from transcripts dictionary if mapping file fails
            if initial_gene_ids:
                filter_log_df = pd.DataFrame({"GeneID": list(initial_gene_ids)})
                filter_log_df["FilterReason"] = "Passed"
                filter_log_df.set_index("GeneID", inplace=True)
                logger.warning(
                    f"Fallback: Initialized filter log with {len(filter_log_df)} genes directly from transcript keys."
                )
            else:
                logger.error(
                    "Could not initialize filter log: No initial transcripts found and mapping file failed."
                )

    elif initial_gene_ids:
        # Fallback if mapping file doesn't exist but transcripts do
        filter_log_df = pd.DataFrame({"GeneID": list(initial_gene_ids)})
        filter_log_df["FilterReason"] = "Passed"
        filter_log_df.set_index("GeneID", inplace=True)
        logger.warning(
            f"Step 2 mapping file not found. Initialized filter log with {len(filter_log_df)} genes directly from transcript keys."
        )
    else:
        logger.error(
            "Could not initialize filter log: No initial transcripts or mapping file found."
        )

    return filter_log_df


def setup_file_and_console_logging(args, output_dir):
    """Set up file and console logging with custom levels."""
    log_file_path = os.path.join(output_dir, f"tusco_{args.species}.log")
    
    # Handle backward compatibility with --verbose flag
    if args.verbose and args.log_level == "info":
        console_level = "debug"
    else:
        console_level = args.log_level
    
    # Setup the enhanced logging
    setup_logging(
        log_level=console_level.upper(),
        file_log_level=args.file_log_level.upper(),
        log_file=log_file_path
    )
    
    logger.verbose(f"Log file created at: {os.path.abspath(log_file_path)}")
    return log_file_path


def display_welcome_banner(species: str) -> None:
    """Log a concise start message for the selected species."""
    logger.info("TUSCO Selector – species: %s", species.upper())


def download_resources(species, cache_dir):
    """Step 1: Download required resources."""
    # Import here to avoid requiring these dependencies for tests
    from tusco_selector.downloader import fetch_all

    resources_module = _load_external_resources(species)
    files = fetch_all(resources_module.EXTERNAL_RESOURCES, cache_dir)

    print_success(f"Downloaded {len(files)} resources to {cache_dir}")

    # Identify first GTF as template
    gtf_files = [f for f in files if f.lower().endswith((".gtf", ".gtf.gz"))]
    if not gtf_files:
        print_warning(
            "No GTF files found in the downloaded resources – aborting after step 1"
        )
        return None, None

    template_gtf = gtf_files[0]
    logger.info(f"Using {template_gtf} as template GTF for sub-setting")

    return files, template_gtf


def select_single_isoform_genes(files, species, output_dir, template_gtf, cache_dir):
    """Step 2: Select matched single-isoform genes."""
    from tusco_selector.pipeline.single_isoform import (
        filter_genes_present_in_mapping,
        select_matched_single_isoforms,
    )
    from tusco_selector.utils import (
        build_transcripts_dict_from_gtf,
        subset_gtf_by_gene_ids,
        subset_mapping_file,
    )

    step2_gtf, step2_mapping = _get_step_paths(output_dir, 2, "single_isoform")

    # Get path to the mapping file
    mapping_file = _mapping_file_path(cache_dir)

    if os.path.exists(step2_gtf):
        print_info("Found existing Step-2 GTF. Loading transcripts.")
        transcripts = build_transcripts_dict_from_gtf(step2_gtf)
    else:
        transcripts = select_matched_single_isoforms(files, species=species)

        # Filter out non-mRNA transcripts using the mapping file
        transcripts = filter_genes_present_in_mapping(transcripts, mapping_file)

    # Save step-2 outputs reflecting the selected single-isoform set
    _save_step_outputs(
        source_gtf_path=template_gtf,
        gene_ids_to_keep=set(transcripts.keys()),
        output_gtf_path=step2_gtf,
        output_mapping_path=step2_mapping,
        main_mapping_file=mapping_file,
        step_description="single isoform selection",
    )

    print_success(f"Selected {len(transcripts)} single-isoform genes")
    return transcripts, step2_gtf


def check_splice_and_tss(
    transcripts,
    cache_dir,
    species,
    output_dir,
    step2_gtf,
    filter_log_df,
    splice_tss_thresholds,
):
    """Step 3: Check splice junctions and TSS."""
    from tusco_selector.pipeline.splice_tss_check import (
        _load_tss_bed,
        check_splice_junctions_tss,
    )
    from tusco_selector.utils import (
        build_transcripts_dict_from_gtf,
        subset_gtf_by_gene_ids,
        subset_mapping_file,
    )

    step3_gtf, step3_mapping = _get_step_paths(output_dir, 3, "splice_tss_filtered")

    # Simplified: omit verbose per-gene debug logging in Step 3

    # Get path to the mapping file
    mapping_file = _mapping_file_path(cache_dir)

    force_step3 = str(os.environ.get("TUSCO_FORCE_STEP3", "")).lower() in {"1", "true", "yes"}
    if os.path.exists(step3_gtf) and not force_step3:
        print_info(f"Found existing Step-3 GTF. Loading transcripts.")
        transcripts = build_transcripts_dict_from_gtf(step3_gtf)

        _subset_mapping_if_missing(
            transcripts, mapping_file, step3_mapping, "splice junction and TSS checks"
        )
        # No additional debug output for cached Step-3
    else:
        # Load resources module for this species
        resources_module = _load_external_resources(species)

        # Load splice junction and TSS data
        junctions = None
        tss_df = None

        try:
            splicing_tss_resources = getattr(
                resources_module, "EXTERNAL_RESOURCES", {}
            ).get("splicing_tss", {})

            # Load recount3 data if available
            recount3_info = splicing_tss_resources.get("recount3", {})
            if recount3_info:
                recount3_file = os.path.join(
                    cache_dir, "splicing_tss", recount3_info.get("output", "")
                )
                if os.path.exists(recount3_file):
                    junctions = load_junctions_bed(recount3_file, species)

            # Load refTSS data if available
            reftss_info = splicing_tss_resources.get("refTSS", {})
            if reftss_info:
                reftss_file = os.path.join(
                    cache_dir, "splicing_tss", reftss_info.get("output", "")
                )
                if os.path.exists(reftss_file):
                    tss_df = _load_tss_bed(reftss_file)

        except Exception as e:
            print_warning(f"Error loading splice/TSS data: {e}")
            logger.exception("Failed to load splice/TSS data")

        # Run checks
        transcripts, removed_info = check_splice_junctions_tss(
            transcripts,
            junctions=junctions,
            tss_df=tss_df,
            novel_threshold=splice_tss_thresholds["novel_threshold"],
            min_novel_length=splice_tss_thresholds["min_novel_length"],
            tss_scope=splice_tss_thresholds["tss_scope"],
            tss_region_bp=splice_tss_thresholds.get("tss_region_bp"),
        )

        # Update filter log using the helper function
        filter_log_df = _update_filter_log(
            filter_log_df, removed_info, "Splice/TSS Filter"
        )

        # Save using the helper function
        _save_step_outputs(
            source_gtf_path=step2_gtf,
            gene_ids_to_keep=set(transcripts.keys()),
            output_gtf_path=step3_gtf,
            output_mapping_path=step3_mapping,
            main_mapping_file=mapping_file,
            step_description="splice junction and TSS checks",
        )

    print_success("Completed splice junction and TSS checks")
    # Return filter_log_df
    return transcripts, step3_gtf, filter_log_df


def apply_introverse_filter(
    transcripts, species, cache_dir, output_dir, step3_gtf, filter_log_df
):
    """Step 4: Apply IntroVerse filtering (human-only)."""
    from tusco_selector.utils import (
        build_transcripts_dict_from_gtf,
        subset_gtf_by_gene_ids,
        subset_mapping_file,
    )

    step4_gtf, step4_mapping = _get_step_paths(output_dir, 4, "introverse_filtered")

    # Get path to the mapping file
    mapping_file = _mapping_file_path(cache_dir)

    if species == "hsa":
        try:
            from tusco_selector.pipeline.introverse_filter import (
                apply_introverse_filter as apply_introverse_filter_func,
            )
            from tusco_selector.pipeline.introverse_filter import (
                read_csv_and_get_ensembl_ids,
            )

            resources_module = _load_external_resources(species)
            splicing_tss_resources = getattr(
                resources_module, "EXTERNAL_RESOURCES", {}
            ).get("splicing_tss", {})
            introverse_info = splicing_tss_resources.get("introverse", {})
            introverse_file = os.path.join(
                cache_dir, "splicing_tss", introverse_info.get("output", "")
            )

            if os.path.exists(introverse_file):
                if os.path.exists(step4_gtf):
                    print_info(f"Found existing Step-4 GTF. Loading transcripts.")
                    transcripts = build_transcripts_dict_from_gtf(step4_gtf)

                    _subset_mapping_if_missing(
                        transcripts, mapping_file, step4_mapping, "IntroVerse filtering"
                    )
                else:
                    introverse_genes = read_csv_and_get_ensembl_ids(introverse_file)
                    # Modify apply_introverse_filter_func to return removed IDs
                    transcripts, stats, removed_ids = apply_introverse_filter_func(
                        transcripts, introverse_genes
                    )

                    # Update filter log using the helper function
                    filter_log_df = _update_filter_log(
                        filter_log_df, removed_ids, "IntroVerse Filter"
                    )

                    # Save using the helper function
                    _save_step_outputs(
                        source_gtf_path=step3_gtf,
                        gene_ids_to_keep=set(transcripts.keys()),
                        output_gtf_path=step4_gtf,
                        output_mapping_path=step4_mapping,
                        main_mapping_file=mapping_file,
                        step_description="IntroVerse filtering",
                    )
            else:
                print_warning(f"IntroVerse CSV not found. Skipping filtering.")
                step4_gtf = step3_gtf

                # Copy mapping file from step 3 if it exists
                step3_mapping = os.path.join(
                    output_dir, "step3_splice_tss_filtered_mapping.tsv"
                )
                _copy_file_if_missing(
                    step3_mapping,
                    step4_mapping,
                    "Copying step 3 mapping to step 4 as IntroVerse was skipped.",
                )
                # Ensure the mapping file reflects previous step
                _, prev_step_mapping = _get_step_paths(
                    output_dir, 3, "splice_tss_filtered"
                )
                _copy_file_if_missing(
                    prev_step_mapping,
                    step4_mapping,
                    "Copying step 3 mapping to step 4 as IntroVerse was skipped or failed.",
                )
        except Exception as e:
            print_warning(f"Error during IntroVerse filtering: {e}")
            logger.exception("IntroVerse filtering failed")
            step4_gtf = step3_gtf

            # Copy mapping file from step 3 if it exists
            step3_mapping = os.path.join(
                output_dir, "step3_splice_tss_filtered_mapping.tsv"
            )
            _copy_file_if_missing(
                step3_mapping,
                step4_mapping,
                "Copying step 3 mapping to step 4 due to IntroVerse error.",
            )
    else:
        # Non-human species - skip this step
        print_info("Skipping IntroVerse filtering (not applicable for this species)")
        step4_gtf = step3_gtf

        # Copy mapping file from step 3 if it exists
        step3_mapping = os.path.join(
            output_dir, "step3_splice_tss_filtered_mapping.tsv"
        )
        _copy_file_if_missing(
            step3_mapping,
            step4_mapping,
            "Copying step 3 mapping to step 4 for non-human species.",
        )
        # If step is skipped, copy previous step's mapping file if it doesn't exist
        _, prev_step_mapping = _get_step_paths(output_dir, 3, "splice_tss_filtered")
        _copy_file_if_missing(
            prev_step_mapping,
            step4_mapping,
            "Copying step 3 mapping to step 4 as IntroVerse was skipped for non-human species.",
        )

    # Return filter_log_df
    return transcripts, step4_gtf, filter_log_df


def filter_by_expression(
    transcripts,
    species,
    cache_dir,
    output_dir,
    step4_gtf,
    filter_log_df,
    expression_thresholds,
    expression_source: str = "auto",
    encode_tissues: List[str] | None = None,
):
    """Step 5: Filter by expression level."""
    from tusco_selector.utils import (
        build_transcripts_dict_from_gtf,
        subset_gtf_by_gene_ids,
        subset_mapping_file,
    )

    step5_gtf, step5_mapping = _get_step_paths(output_dir, 5, "expression_filtered")

    # Get path to the mapping file
    mapping_file = _mapping_file_path(cache_dir)
    initial_genes = set(transcripts.keys())  # Genes entering this step

    if os.path.exists(step5_gtf):
        print_info(f"Found existing Step-5 GTF. Loading transcripts.")
        transcripts = build_transcripts_dict_from_gtf(step5_gtf)

        _subset_mapping_if_missing(
            transcripts, mapping_file, step5_mapping, "expression filtering (universal)"
        )
    else:
        from tusco_selector.pipeline.expression_filter import (
            get_universal_high_expression_genes_bgee_with_ids,
        )

        # Choose expression source based on user selection
        selected_source = (expression_source or "auto").lower()

        if selected_source == "gtex":
            if species != "hsa":
                print_warning(
                    f"GTEx selected explicitly but species is {species}. Proceeding with GTEx; results may be empty if data is unavailable."
                )
            print_info("Processing GTEx expression data...")
            universal_genes, high_exp_per_tissue = _process_gtex_expression_data_unified(
                cache_dir,
                species,
                expression_thresholds.get(
                    "gtex",
                    {
                        "prevalence_expression_cutoff": 0.1,
                        "median_expression_cutoff": 1.0,
                        "prevalence_threshold": 0.95,
                    },
                ),
            )
        elif selected_source == "bgee":
            print_info(f"Processing Bgee expression data for {species}...")
            universal_genes, high_exp_per_tissue = _process_bgee_expression_data_unified(
                cache_dir, species, expression_thresholds["bgee"]
            )
        else:  # auto behavior (new hybrid approach for human)
            if species == "hsa":
                # Use hybrid approach: Bgee for universal genes, GTEx for tissue-specific genes
                universal_genes, high_exp_per_tissue = _process_hybrid_expression_data_human(
                    cache_dir, species, expression_thresholds
                )
            else:
                print_info(f"Processing Bgee expression data for {species}...")
                universal_genes, high_exp_per_tissue = _process_bgee_expression_data_unified(
                    cache_dir, species, expression_thresholds["bgee"]
                )

        if not universal_genes:
            print_warning(
                "No genes met the universal expression criteria – retaining all genes for downstream steps"
            )
            universal_genes = initial_genes
            high_exp_per_tissue = {}

        # Try to load housekeeping genes if available
        expr_dir = os.path.join(cache_dir, "expression")
        housekeeping_genes_file = os.path.join(expr_dir, "hrt_atlas.housekeeping_genes.csv")
        if os.path.exists(housekeeping_genes_file):
            try:
                from tusco_selector.pipeline.expression_filter import load_housekeeping_genes
                housekeeping_genes = load_housekeeping_genes(housekeeping_genes_file)
                if housekeeping_genes:
                    original_count = len(universal_genes)
                    universal_genes.update(housekeeping_genes)
                    added_count = len(universal_genes) - original_count
                    print_success(f"Added {added_count} housekeeping genes to universal genes set")
            except Exception as e:
                print_warning(f"Error processing housekeeping genes: {e}")

        # For universal genes, require passing BOTH expression source and Alphagenome filters
        # First get expression-based universal genes
        expression_universal_genes = universal_genes & initial_genes
        
        # Determine which expression source was used for logging
        if selected_source == "gtex":
            expression_source_name = "GTEx"
        elif selected_source == "bgee":
            expression_source_name = "Bgee"
        else:  # auto mode
            if species == "hsa":
                expression_source_name = "Bgee (universal) + GTEx (tissue-specific)"
            else:
                expression_source_name = "Bgee"
        print_info(f"{expression_source_name} universal genes (expression-based): {len(expression_universal_genes)}")
        
        # Apply Alphagenome filter to get Alphagenome universal genes  
        alphagenome_data, tissue_mapping = _load_alphagenome_data(cache_dir, species)
        if alphagenome_data and tissue_mapping:
            # Classify genes by exon counts for Alphagenome filtering
            try:
                from tusco_selector.utils import classify_genes_by_exons
                single_exon_genes, multi_exon_genes = classify_genes_by_exons(step4_gtf)
            except Exception:
                single_exon_genes, multi_exon_genes = set(), set()
            
            # Apply Alphagenome universal filter to all initial genes
            alphagenome_universal_genes, alphagenome_removed = _apply_alphagenome_filter(
                initial_genes,
                alphagenome_data,
                tissue_mapping,
                expression_thresholds["alphagenome"]["universal"],
                single_exon_genes,
                multi_exon_genes,
            )
            
            # Analyze removal reasons for Alphagenome
            thr = expression_thresholds["alphagenome"]["universal"]
            alpha_reason_counts = {
                "median": 0,
                "prevalence": 0,
                "splice_purity": 0,
                "not_in_db": 0,
                "no_tissue": 0,
            }
            for reason in alphagenome_removed.values():
                r = str(reason).lower()
                if "median" in r and "rpkm" in r:
                    alpha_reason_counts["median"] += 1
                elif r.startswith("prevalence"):
                    alpha_reason_counts["prevalence"] += 1
                elif r.startswith("splice purity"):
                    alpha_reason_counts["splice_purity"] += 1
                elif "not in alphagenome" in r:
                    alpha_reason_counts["not_in_db"] += 1
                elif "no expression data" in r:
                    alpha_reason_counts["no_tissue"] += 1
            
            logger.info("=" * 60)
            logger.info("Expression Filter Summary (Step 5):")
            logger.info(f"  {expression_source_name} universal genes: {len(expression_universal_genes):,}")
            logger.info(f"  Alphagenome universal genes: {len(alphagenome_universal_genes):,}")
            logger.info(f"  Alphagenome filter breakdown:")
            if alpha_reason_counts["median"] > 0:
                logger.info(f"    - Below median RPKM: {alpha_reason_counts['median']:,}")
            if alpha_reason_counts["prevalence"] > 0:
                logger.info(f"    - Below prevalence: {alpha_reason_counts['prevalence']:,}")
            if alpha_reason_counts["splice_purity"] > 0:
                logger.info(f"    - Below splice purity: {alpha_reason_counts['splice_purity']:,}")
            if alpha_reason_counts["not_in_db"] > 0:
                logger.info(f"    - Not in Alphagenome: {alpha_reason_counts['not_in_db']:,}")
            if alpha_reason_counts["no_tissue"] > 0:
                logger.info(f"    - No tissue data: {alpha_reason_counts['no_tissue']:,}")
            
            # Universal genes must pass BOTH expression source AND Alphagenome criteria
            passing_universal_genes = expression_universal_genes & alphagenome_universal_genes
            logger.info(f"  Final universal genes (intersection): {len(passing_universal_genes):,}")
            logger.info(f"  Total genes removed in Step 5: {len(initial_genes) - len(passing_universal_genes):,}")
            logger.info("=" * 60)
            
        else:
            # Fallback: only expression source filtering available
            passing_universal_genes = expression_universal_genes
            print_warning(f"No Alphagenome data available - universal genes based only on {expression_source_name}")

        # Determine removed IDs from combined expression and Alphagenome filtering  
        removed_ids = initial_genes - passing_universal_genes

        # Update filter log using the helper function
        filter_log_df = _update_filter_log(
            filter_log_df, removed_ids, "Expression+Alphagenome Filter"
        )
        
        # Optionally augment tissue-specific sets with ENCODE data before Alphagenome filtering
        try:
            high_exp_per_tissue = _augment_tissue_specific_with_encode(
                cache_dir,
                species,
                high_exp_per_tissue,
                expression_thresholds,
                encode_tissues=encode_tissues,
            )
        except Exception as e:
            print_warning(f"Skipping ENCODE augmentation due to error: {e}")

        # Apply tissue-specific filtering based on data source availability
        filtered_high_exp_per_tissue = {}
        try:
            alphagenome_data, tissue_mapping = _load_alphagenome_data(cache_dir, species)
            
            if high_exp_per_tissue:
                # Determine which data sources are available for each tissue
                tissue_list = list(high_exp_per_tissue.keys())  # (tissue_id, tissue_name) tuples
                
                # Determine the expression source used for tissue-specific genes
                tissue_expression_source = "bgee"  # default
                if selected_source == "gtex":
                    tissue_expression_source = "gtex"
                elif selected_source == "auto" and species == "hsa":
                    tissue_expression_source = "gtex"  # hybrid mode uses GTEx for tissue-specific
                
                tissue_data_sources = _determine_tissue_data_sources(
                    tissue_list, alphagenome_data, tissue_mapping, tissue_expression_source
                )
                
                tissue_thresholds = expression_thresholds["alphagenome"]["tissue_specific"]
                kept_counts = []
                
                print_info(f"\n=== Processing {len(high_exp_per_tissue)} tissues with tissue-specific filtering ===")
                print_info(f"Genes entering tissue-specific filtering: {len(initial_genes)}")
                print_info(f"Expression thresholds: TPM > {tissue_thresholds.get('expression_threshold', 'N/A')}, Splice ratio < {tissue_thresholds.get('splice_ratio_threshold', 'N/A')}")
                
                for tissue_info, expression_genes in high_exp_per_tissue.items():
                    anatomical_id, tissue_name = tissue_info
                    candidate_genes = expression_genes & initial_genes
                    
                    if anatomical_id not in tissue_data_sources:
                        logger.warning(f"Tissue {tissue_name} ({anatomical_id}) not found in data source mapping")
                        continue
                    
                    data_sources = tissue_data_sources[anatomical_id]['sources']
                    
                    # Apply appropriate filtering based on available data sources
                    passing_genes_in_tissue = _apply_tissue_specific_filtering(
                        candidate_genes,
                        anatomical_id, 
                        tissue_name,
                        data_sources,
                        expression_genes,
                        alphagenome_data,
                        tissue_mapping,
                        tissue_thresholds,
                    )
                    
                    # Always keep tissue (even if empty) to ensure all tissues get output files
                    filtered_high_exp_per_tissue[tissue_info] = passing_genes_in_tissue
                    if passing_genes_in_tissue:
                        kept_counts.append((tissue_name, anatomical_id, len(passing_genes_in_tissue), '+'.join(data_sources)))
                
                if kept_counts:
                    kept_counts.sort(key=lambda x: x[2], reverse=True)
                    print_info(f"\n=== Tissue-specific filtering summary ===")
                    print_info(f"Total tissues processed: {len(high_exp_per_tissue)}")
                    print_info(f"Tissues with genes: {len(kept_counts)}")
                    print_info(f"Top tissues by gene count:")
                    for tissue_name, tid, count, sources in kept_counts[:10]:  # Show top 10
                        print_info(f"  • {tissue_name}: {count} genes [{sources}]")
                    if len(kept_counts) > 10:
                        print_info(f"  ... and {len(kept_counts) - 10} more tissues with genes")
                    
                    empty_tissues = len(filtered_high_exp_per_tissue) - len(kept_counts)
                    if empty_tissues > 0:
                        print_info(f"\nTissues with no genes after filtering: {empty_tissues}")
                        print_info(f"  (Empty output files will be created for these tissues)")
                else:
                    print_warning("\nNo tissues have genes after tissue-specific filtering")
            else:
                print_info("No tissue-specific data from expression filtering")
                
        except Exception as e:
            # Fail open: keep original sets if filtering errors out
            print_warning(f"Tissue-specific filtering skipped due to error: {e}")
            filtered_high_exp_per_tissue = high_exp_per_tissue

        # Save universal high-expression genes and tissue-specific outputs (after Alphagenome tissue-specific filter)
        _save_expression_outputs_with_ids(
            passing_universal_genes,
            filtered_high_exp_per_tissue,
            initial_genes,
            step4_gtf,
            output_dir,
            mapping_file,
            species,
        )

        # Update transcripts dict to only contain the universally passing genes for the next step
        if passing_universal_genes:
            transcripts = {
                g: info
                for g, info in transcripts.items()
                if g in passing_universal_genes
            }
        else:
            transcripts = {}  # No genes pass

        # Save final outputs for this step (genes passing the universal filter)
        subset_gene_ids = set(transcripts.keys())
        subset_gtf_by_gene_ids(step4_gtf, subset_gene_ids, step5_gtf)
        # Subset mapping file if it exists
        if os.path.exists(mapping_file):
            subset_mapping_file(
                list(transcripts.keys()),
                mapping_file,
                step5_mapping,
                "expression filtering (universal)",
            )

    final_gene_count = len(transcripts)
    print_success(f"Expression filtering complete: {final_gene_count} genes remain")
    # Return filter_log_df
    return transcripts, step5_gtf, filter_log_df


def apply_alphagenome_filter(
    transcripts,
    species,
    cache_dir,
    output_dir,
    step5_gtf,
    filter_log_df,
    alphagenome_thresholds,
):
    """Step 6: Apply Alphagenome filtering criteria.
    
    NOTE: Universal gene Alphagenome filtering is now applied in Step 5 (expression filtering)
    to ensure universal genes must pass BOTH Bgee AND Alphagenome criteria.
    This step now only handles any remaining Alphagenome filtering for non-universal genes.
    """
    from tusco_selector.utils import (
        build_transcripts_dict_from_gtf,
        subset_gtf_by_gene_ids,
        subset_mapping_file,
    )

    step6_gtf, step6_mapping = _get_step_paths(output_dir, 6, "alphagenome_filtered")

    # Get path to the mapping file
    mapping_file = _mapping_file_path(cache_dir)
    initial_genes = set(transcripts.keys())  # Genes entering this step

    if os.path.exists(step6_gtf):
        print_info(f"Found existing Step-6 GTF. Loading transcripts.")
        transcripts = build_transcripts_dict_from_gtf(step6_gtf)

        _subset_mapping_if_missing(
            transcripts, mapping_file, step6_mapping, "Alphagenome filtering"
        )
    else:
        # Apply Alphagenome filter
        print_info("Applying Alphagenome filter (RPKM, prevalence, splice purity)...")
        alphagenome_data, tissue_mapping = _load_alphagenome_data(cache_dir, species)

        if alphagenome_data and tissue_mapping:
            # Classify genes by exon counts using step-5 GTF
            try:
                from tusco_selector.utils import classify_genes_by_exons
                single_exon_genes, multi_exon_genes = classify_genes_by_exons(step5_gtf)
            except Exception:
                single_exon_genes, multi_exon_genes = set(), set()

            if logger.isEnabledFor(logging.DEBUG):
                # Optional debug analytics
                all_alphagenome_genes = {gene["gene_id"] for gene in alphagenome_data}
                debug_passing_genes, debug_removed_genes = _apply_alphagenome_filter(
                    all_alphagenome_genes,
                    alphagenome_data,
                    tissue_mapping,
                    alphagenome_thresholds["universal"],
                    single_exon_genes,
                    multi_exon_genes,
                )

                independent_results = _analyze_independent_alphagenome_filters(
                    all_alphagenome_genes,
                    alphagenome_data,
                    tissue_mapping,
                    alphagenome_thresholds["universal"],
                    single_exon_genes,
                    multi_exon_genes,
                )

                logger.debug(
                    "Alphagenome debug: total=%d, pass=%d, removed=%d, pass_rate=%.1f%%",
                    len(all_alphagenome_genes),
                    len(debug_passing_genes),
                    len(debug_removed_genes),
                    (len(debug_passing_genes) / max(1, len(all_alphagenome_genes))) * 100.0,
                )
                if debug_removed_genes:
                    reason_counts = {}
                    for reason in debug_removed_genes.values():
                        key = reason.split(" ")[0] if " " in reason else reason
                        reason_counts[key] = reason_counts.get(key, 0) + 1
                    logger.debug("Sequential removal breakdown: %s", reason_counts)

                if independent_results.get("total_analyzed", 0) > 0:
                    logger.debug("Independent breakdown: %s", independent_results)

            # Apply Alphagenome filter to the genes from previous step
            alpha_passing_genes, alpha_removed_genes = _apply_alphagenome_filter(
                initial_genes,
                alphagenome_data,
                tissue_mapping,
                alphagenome_thresholds["universal"],
                single_exon_genes,
                multi_exon_genes,
            )

            # Update filter log for Alphagenome filtering
            filter_log_df = _update_filter_log(
                filter_log_df, alpha_removed_genes, "Alphagenome Filter"
            )

            # Summarize detailed removal reasons
            try:
                thr = alphagenome_thresholds["universal"]
                reason_counts = {
                    "median": 0,
                    "prevalence": 0,
                    "splice_purity": 0,
                    "not_in_alphagenome": 0,
                    "no_tissue_data": 0,
                    "other": 0,
                }
                for reason in alpha_removed_genes.values():
                    r = str(reason).lower()
                    if "median" in r and "rpkm" in r:
                        reason_counts["median"] += 1
                    elif r.startswith("prevalence"):
                        reason_counts["prevalence"] += 1
                    elif r.startswith("splice purity"):
                        reason_counts["splice_purity"] += 1
                    elif "not in alphagenome data" in r:
                        reason_counts["not_in_alphagenome"] += 1
                    elif "no expression data for mapped tissues" in r:
                        reason_counts["no_tissue_data"] += 1
                    else:
                        reason_counts["other"] += 1

                # Log detailed breakdown
                logger.info("=" * 60)
                logger.info("Alphagenome Filter Summary (Step 6):")
                logger.info(f"  Genes entering filter: {len(initial_genes):,}")
                logger.info(f"  Genes passing filter: {len(alpha_passing_genes):,}")
                logger.info(f"  Genes removed: {len(alpha_removed_genes):,}")
                logger.info("  Removal breakdown:")
                logger.info(f"    - Below median RPKM threshold: {reason_counts['median']:,}")
                logger.info(f"      (single-exon >= {thr['single_exon_median_rpkm_threshold']}, multi-exon >= {thr['multi_exon_median_rpkm_threshold']})")
                logger.info(f"    - Below prevalence threshold: {reason_counts['prevalence']:,}")
                logger.info(f"      (single-exon >= {thr['single_exon_prevalence_threshold']} @ RPKM > {thr['single_exon_expression_threshold']})")
                logger.info(f"      (multi-exon >= {thr['multi_exon_prevalence_threshold']} @ RPKM > {thr['multi_exon_expression_threshold']})")
                logger.info(f"    - Below splice purity threshold: {reason_counts['splice_purity']:,}")
                logger.info(f"      (>= {thr['splice_purity_threshold']} tissues with splice-ratio < {thr['splice_ratio_threshold']})")
                logger.info(f"    - Not in Alphagenome database: {reason_counts['not_in_alphagenome']:,}")
                logger.info(f"    - No expression data for mapped tissues: {reason_counts['no_tissue_data']:,}")
                if reason_counts['other'] > 0:
                    logger.info(f"    - Other reasons: {reason_counts['other']:,}")
                logger.info("=" * 60)
            except Exception as e:
                logger.debug(f"Failed to summarize Alphagenome removal reasons: {e}")

            print_success(
                f"Alphagenome filter: {len(alpha_passing_genes)} genes pass all criteria"
            )

            # Update transcripts dict to only contain passing genes
            transcripts = {
                g: info
                for g, info in transcripts.items()
                if g in alpha_passing_genes
            }
        else:
            print_warning("Alphagenome filtering skipped - data not available")

        # Save using the helper function
        _save_step_outputs(
            source_gtf_path=step5_gtf,
            gene_ids_to_keep=set(transcripts.keys()),
            output_gtf_path=step6_gtf,
            output_mapping_path=step6_mapping,
            main_mapping_file=mapping_file,
            step_description="Alphagenome filtering",
        )

    final_gene_count = len(transcripts)
    print_success(f"Alphagenome filtering complete: {final_gene_count} genes remain")
    return transcripts, step6_gtf, filter_log_df


def apply_manual_filters(
    transcripts, cache_dir, output_dir, step6_gtf, species, filter_log_df
):
    """Step 7: Apply manual filters and generate final annotation."""
    from tusco_selector.pipeline.manual_filter import (
        apply_manual_filter,
        load_manual_filter,
    )
    from tusco_selector.utils import (
        build_transcripts_dict_from_gtf,
        subset_gtf_by_gene_ids,
        subset_mapping_file,
    )

    step7_gtf, step7_mapping = _get_step_paths(output_dir, 7, "final_filtered")

    # Get path to the mapping file
    mapping_file = _mapping_file_path(cache_dir)

    if os.path.exists(step7_gtf):
        print_info(f"Found existing Step-7 GTF. Loading final transcripts.")
        transcripts = build_transcripts_dict_from_gtf(step7_gtf)

        _subset_mapping_if_missing(
            transcripts, mapping_file, step7_mapping, "manual filtering"
        )
    else:
        # Apply manual filters if present
        filter_file = os.path.join(cache_dir, "manual_filtered.txt")
        removed_ids = set()
        if os.path.isfile(filter_file):
            filter_genes_set = load_manual_filter(Path(filter_file))
            # Modify apply_manual_filter to return removed IDs
            transcripts, removed_ids = apply_manual_filter(
                transcripts, filter_genes_set
            )
            print_success(f"Applied manual filters from {filter_file}")

            # Update filter log
            if filter_log_df is not None:
                for gene_id in removed_ids:
                    if gene_id in filter_log_df.index:
                        # Only update if the gene hasn't been filtered already
                        if filter_log_df.loc[gene_id, "FilterReason"] == "Passed":
                            filter_log_df.loc[gene_id, "FilterReason"] = "Manual Filter"
                    else:
                        logger.warning(
                            f"Gene {gene_id} from Manual filter not found in initial filter log."
                        )
            # Update filter log using the helper function
            filter_log_df = _update_filter_log(
                filter_log_df, removed_ids, "Manual Filter"
            )
        else:
            print_warning("No manual filter file found, skipping manual filtering")

        # Save using the helper function
        _save_step_outputs(
            source_gtf_path=step6_gtf,
            gene_ids_to_keep=set(transcripts.keys()),
            output_gtf_path=step7_gtf,
            output_mapping_path=step7_mapping,
            main_mapping_file=mapping_file,
            step_description="manual filtering",
        )

    final_gene_count = len(transcripts)
    print_success(
        f"Manual filtering complete: {final_gene_count} genes in final output"
    )
    # Return filter_log_df
    return transcripts, step7_gtf, filter_log_df


def summarize_results(output_dir):
    """Summarize available output files (logged)."""
    logger.info("Available annotation files:")

    output_files = [
        ("Step-2: single-isoform", "step2_single_isoform.gtf.gz"),
        ("Step-3: splice/TSS filtered", "step3_splice_tss_filtered.gtf.gz"),
        ("Step-4: manual filtered", "step4_introverse_filtered.gtf.gz"),
        ("Step-5: expression filtered (includes AlphaGenome)", "step5_expression_filtered.gtf.gz"),
        ("Step-5: universal high expression", "step5_universal_high_expression.gtf.gz"),
        ("Step-6: IntroVerse filtered (final)", "step6_introverse_filtered.gtf.gz"),
        ("Final output", "step7_final_filtered.gtf.gz"),
    ]

    mapping_files = [
        ("Step-2 mapping", "step2_single_isoform_mapping.tsv"),
        ("Step-3 mapping", "step3_splice_tss_filtered_mapping.tsv"),
        ("Step-4 mapping", "step4_introverse_filtered_mapping.tsv"),
        ("Step-5 mapping", "step5_expression_filtered_mapping.tsv"),
        ("Step-5 universal mapping", "step5_universal_high_expression_mapping.tsv"),
        ("Step-6 mapping", "step6_introverse_filtered_mapping.tsv"),
        ("Final mapping", "step7_final_filtered_mapping.tsv"),
    ]

    for label, filename in output_files:
        path = os.path.join(output_dir, filename)
        if os.path.exists(path):
            logger.info("- %s → %s", label, os.path.abspath(path))

    logger.info("Available mapping files (Ensembl ↔ RefSeq ↔ NCBI):")
    for label, filename in mapping_files:
        path = os.path.join(output_dir, filename)
        if os.path.exists(path):
            logger.info("- %s → %s", label, os.path.abspath(path))

    logger.info("TUSCO selection complete! Results in: %s", os.path.abspath(output_dir))


# ---------------------------------------------------------------------------
# TUSCO tables generation (final deliverables)
# ---------------------------------------------------------------------------


def generate_tusco_deliverables(output_dir: str, *, human: str = "hsa", mouse: str = "mmu", species: str | None = None) -> None:
    """Generate TUSCO TSV deliverables from step-7 outputs in ``output_dir``.

    Looks for ``step7_final_filtered_mapping.tsv`` and ``step7_final_filtered.gtf.gz``
    and writes three tables: ``tusco_{human}{mouse}.tsv``, ``..._multi_exon.tsv``,
    and ``..._single_exon.tsv`` to the same directory.
    """

    mapping_path = os.path.join(output_dir, "step7_final_filtered_mapping.tsv")
    gtf_path = os.path.join(output_dir, "step7_final_filtered.gtf.gz")

    if not os.path.exists(mapping_path) or not os.path.exists(gtf_path):
        print_warning(
            "Final step-7 files not found; skipping TUSCO table generation"
        )
        return

    try:
        # Decide output naming scheme
        species_label = None
        if species:
            sp = species.lower()
            if sp == "hsa":
                species_label = "human"
            elif sp == "mmu":
                species_label = "mouse"

        all_p, multi_p, single_p = generate_tusco_tables(
            mapping_path,
            gtf_path,
            out_dir=output_dir,
            human=human,
            mouse=mouse,
            species_label=species_label,
        )
        print_success(
            f"Generated TUSCO tables: {os.path.basename(all_p)}, "
            f"{os.path.basename(multi_p)}, {os.path.basename(single_p)}"
        )
    except Exception as exc:
        logger.exception("Failed generating TUSCO tables: %s", exc)
        print_warning(f"Failed generating TUSCO tables: {exc}")


def export_species_named_gtf(output_dir: str, species: str) -> None:
    """Export final step-7 GTF to species-named file ``tusco_human.gtf`` or ``tusco_mouse.gtf``.

    Only applies to human (hsa) and mouse (mmu). If the step-7 compressed GTF is
    missing, this is a no-op with a warning.
    """
    src_gtf_gz = os.path.join(output_dir, "step7_final_filtered.gtf.gz")
    if not os.path.exists(src_gtf_gz):
        print_warning(
            f"Final GTF not found, cannot write species-named file: {src_gtf_gz}"
        )
        return

    label = "human" if species == "hsa" else "mouse" if species == "mmu" else None
    if label is None:
        # Only requested for human/mouse
        return

    dest_gtf = os.path.join(output_dir, f"tusco_{label}.gtf")

    try:
        import gzip
        import shutil

        with gzip.open(src_gtf_gz, mode="rt") as fin, open(dest_gtf, mode="wt") as fout:
            shutil.copyfileobj(fin, fout)
        print_success(f"Wrote species-named GTF: {os.path.abspath(dest_gtf)}")
    except Exception as exc:
        logger.exception("Failed writing species-named GTF: %s", exc)
        print_warning(f"Failed writing species-named GTF: {exc}")

# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main(cli_args: list[str] | None = None):
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Download resources and output GTFs for TUSCO gene selection"
    )
    parser.add_argument(
        "species",
        choices=["hsa", "mmu", "dre"],
        help="Three-letter species code: hsa (human), mmu (mouse), dre (zebrafish)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="output_gtfs",
        help="Directory to list or copy resulting GTF files",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show detailed progress information (deprecated, use --log-level)",
    )
    parser.add_argument(
        "--log-level",
        choices=["trace", "debug", "verbose", "info", "warning", "error"],
        default="info",
        help="Console logging verbosity (default: info)",
    )
    parser.add_argument(
        "--file-log-level",
        choices=["trace", "debug", "verbose", "info", "warning", "error"],
        default="verbose",
        help="File logging verbosity (default: verbose)",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Show progress bars for long operations",
    )

    # Optional: ENCODE augmentation selection
    parser.add_argument(
        "--encode-tissues",
        nargs="+",
        default=[],
        help=(
            "List of tissue names to augment tissue-specific sets with ENCODE (e.g., 'kidney'). "
            "Tissues are resolved via a built-in tissue→UBERON mapping."
        ),
    )

    # Step 3: Splice junction and TSS filtering thresholds
    splice_group = parser.add_argument_group(
        "Step 3: Splice junction and TSS filtering"
    )
    splice_group.add_argument(
        "--novel-threshold",
        type=float,
        default=0.01,
        help="Novel junction coverage threshold relative to average annotated coverage (default: 0.01)",
    )
    splice_group.add_argument(
        "--min-novel-length",
        type=int,
        default=80,
        help="Minimum length (bp) for non-annotated junctions to be considered (default: 80)",
    )
    splice_group.add_argument(
        "--tss-scope",
        choices=["single", "multi", "both"],
        default="both",
        help=(
            "Apply TSS filtering to: 'single' (single-exon only), 'multi' (multi-exon only), "
            "or 'both' (default)."
        ),
    )
    splice_group.add_argument(
        "--tss-region-check",
        type=int,
        default=0,
        help=(
            "Window size (bp) around the TSS to check for CAGE peaks. "
            "If more than one peak exists within ±N bp, the gene is filtered. "
            "Set 0 to disable and use exon-overlap mode (default: 0)."
        ),
    )

    # Step 5: Expression filtering thresholds
    expr_group = parser.add_argument_group("Step 5: Expression filtering")

    # Allow user to choose expression data source explicitly
    expr_group.add_argument(
        "--expression-source",
        choices=["auto", "bgee", "gtex"],
        default="auto",
        help=(
            "Choose expression data source for filtering: 'gtex' or 'bgee'. "
            "Default 'auto' uses Bgee (universal) + GTEx (tissue-specific) for human and Bgee for others."
        ),
    )

    # GTEx processing thresholds (for human)
    expr_group.add_argument(
        "--gtex-prevalence-expression-cutoff",
        type=float,
        default=0.1,
        help="GTEx TPM threshold for prevalence test (default: 0.1)",
    )
    expr_group.add_argument(
        "--gtex-median-expression-cutoff",
        type=float,
        default=1.0,
        help="GTEx minimum median TPM across tissue samples (default: 1.0)",
    )
    expr_group.add_argument(
        "--gtex-prevalence-threshold",
        type=float,
        default=0.95,
        help="GTEx minimum fraction of samples that must be expressed (default: 0.95)",
    )

    # ENCODE processing thresholds (for optional augmentation)
    expr_group.add_argument(
        "--encode-prevalence-expression-cutoff",
        type=float,
        default=0.1,
        help="ENCODE TPM threshold used for prevalence test (default: 0.1)",
    )
    expr_group.add_argument(
        "--encode-median-expression-cutoff",
        type=float,
        default=1.0,
        help="ENCODE median TPM cutoff per gene across entries (default: 1.0)",
    )
    expr_group.add_argument(
        "--encode-prevalence-threshold",
        type=float,
        default=0.95,
        help="ENCODE minimum fraction of entries with TPM > prevalence cutoff (default: 0.95)",
    )

    # Bgee processing thresholds (for mouse and other species)
    expr_group.add_argument(
        "--min-genes-per-tissue",
        type=int,
        default=25000,
        help="Minimum number of genes a tissue must have in Bgee data to be considered (default: 25000)",
    )
    expr_group.add_argument(
        "--bgee-quality",
        nargs="+",
        default=["gold quality", "silver quality"],
        help="Bgee call quality values to accept (default: 'gold quality' 'silver quality')",
    )
    expr_group.add_argument(
        "--bgee-prevalence-threshold",
        type=float,
        default=0.95,
        help="Minimum fraction of tissues where a gene must be expressed to be considered universal in Bgee (default: 0.95)",
    )

    # Alphagenome filtering thresholds - Universal genes
    universal_group = parser.add_argument_group(
        "Step 5: Alphagenome filtering - Universal genes"
    )
    # RPKM-based median thresholds (with TPM aliases for backward compatibility)
    universal_group.add_argument(
        "--single-exon-median-rpkm-threshold",
        type=float,
        default=1.0,
        help="Minimum median RPKM across tissues for single-exon genes (default: 1.0)",
    )
    universal_group.add_argument(
        "--multi-exon-median-rpkm-threshold",
        type=float,
        default=1.0,
        help="Minimum median RPKM across tissues for multi-exon genes (default: 1.0)",
    )
    universal_group.add_argument(
        "--single-exon-median-tpm-threshold",
        dest="single_exon_median_rpkm_threshold",
        type=float,
        help=argparse.SUPPRESS,
    )
    universal_group.add_argument(
        "--multi-exon-median-tpm-threshold",
        dest="multi_exon_median_rpkm_threshold",
        type=float,
        help=argparse.SUPPRESS,
    )
    universal_group.add_argument(
        "--single-exon-expression-threshold",
        type=float,
        default=0.1,
        help="Minimum RPKM for a gene to be considered expressed in a tissue (single-exon) (default: 0.1)",
    )
    universal_group.add_argument(
        "--multi-exon-expression-threshold",
        type=float,
        default=0.1,
        help="Minimum RPKM for a gene to be considered expressed in a tissue (multi-exon) (default: 0.1)",
    )
    universal_group.add_argument(
        "--single-exon-prevalence-threshold",
        type=float,
        default=0.95,
        help="Minimum fraction of tissues where the single-exon gene must be expressed (default: 0.95)",
    )
    universal_group.add_argument(
        "--multi-exon-prevalence-threshold",
        type=float,
        default=0.95,
        help="Minimum fraction of tissues where the multi-exon gene must be expressed (default: 0.95)",
    )
    universal_group.add_argument(
        "--splice-ratio-threshold",
        type=float,
        default=0.05,
        help="Maximum splice ratio for unannoated junction vs annotated splicing (default: 0.05) - universal genes",
    )
    universal_group.add_argument(
        "--splice-purity-threshold",
        type=float,
        default=0.95,
        help="Minimum fraction of tissues with good splice ratios for universal genes (default: 0.95)",
    )

    # Alphagenome filtering thresholds - Tissue-specific genes
    tissue_group = parser.add_argument_group(
        "Step 5: Alphagenome filtering - Tissue-specific genes"
    )
    tissue_group.add_argument(
        "--tissue-expression-threshold",
        type=float,
        default=1.0,
        help="Minimum RPKM for a gene to be considered expressed in the specific tissue (default: 1.0)",
    )
    tissue_group.add_argument(
        "--tissue-splice-ratio-threshold",
        type=float,
        default=0.01,
        help="Maximum splice ratio for unannotated junction vs annotated splicing (default: 0.01)",
    )

    args = parser.parse_args(cli_args)

    # Create threshold configuration from parsed arguments
    thresholds = {
        "splice_tss": {
            "novel_threshold": args.novel_threshold,
            "min_novel_length": args.min_novel_length,
            "tss_scope": args.tss_scope,
            "tss_region_bp": args.tss_region_check,
        },
        "expression": {
            "gtex": {
                "prevalence_expression_cutoff": args.gtex_prevalence_expression_cutoff,
                "median_expression_cutoff": args.gtex_median_expression_cutoff,
                "prevalence_threshold": args.gtex_prevalence_threshold,
            },
            "bgee": {
                "min_genes_per_tissue": args.min_genes_per_tissue,
                "qual_ok": tuple(args.bgee_quality),
                "prevalence_threshold": args.bgee_prevalence_threshold,
            },
            "encode": {
                "prevalence_expression_cutoff": args.encode_prevalence_expression_cutoff,
                "median_expression_cutoff": args.encode_median_expression_cutoff,
                "prevalence_threshold": args.encode_prevalence_threshold,
            },
            "alphagenome": {
                "universal": {
                    "single_exon_median_rpkm_threshold": args.single_exon_median_rpkm_threshold,
                    "multi_exon_median_rpkm_threshold": args.multi_exon_median_rpkm_threshold,
                    "single_exon_expression_threshold": args.single_exon_expression_threshold,
                    "multi_exon_expression_threshold": args.multi_exon_expression_threshold,
                    "single_exon_prevalence_threshold": args.single_exon_prevalence_threshold,
                    "multi_exon_prevalence_threshold": args.multi_exon_prevalence_threshold,
                    "splice_ratio_threshold": args.splice_ratio_threshold,
                    "splice_purity_threshold": args.splice_purity_threshold,
                },
                "tissue_specific": {
                    "expression_threshold": args.tissue_expression_threshold,
                    "splice_ratio_threshold": args.tissue_splice_ratio_threshold,
                },
            },
        },
    }

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Setup logging
    setup_file_and_console_logging(args, args.output_dir)

    # Welcome banner
    display_welcome_banner(args.species)

    # Determine data directory path
    src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    project_root = os.path.dirname(src_dir)
    cache_dir = os.path.join(project_root, "data", args.species)

    # Step 1: Download resources
    print_step(1, 7, "Downloading required resources")
    files, template_gtf = download_resources(args.species, cache_dir)
    if not files or not template_gtf:
        return

    # Execute pipeline steps 2-7
    transcripts = run_pipeline_steps(args, cache_dir, files, template_gtf, thresholds)

    # Summarize results
    summarize_results(args.output_dir)
    # Export species-named final GTF
    export_species_named_gtf(args.output_dir, args.species)
    # Generate final TUSCO TSV deliverables (human+mouse naming)
    generate_tusco_deliverables(args.output_dir, human="hsa", mouse="mmu", species=args.species)
    # Generate comprehensive statistics file including universal and tissue-specific counts
    _create_comprehensive_statistics_tsv(args.output_dir, args.species)
    logger.info(
        f"TUSCO selection completed successfully – {len(transcripts)} final transcripts"
    )


# Add function to write the filter log
def write_filter_log(filter_log_df, output_dir):
    """Writes the filter log DataFrame to a TSV file."""
    log_file_path = os.path.join(output_dir, "gene_filter_log.tsv")
    try:
        filter_log_df.reset_index().to_csv(log_file_path, sep="	", index=False)
        logger.info(f"Gene filter log saved to: {os.path.abspath(log_file_path)}")
        print_info(f"Gene filter log saved to: {os.path.abspath(log_file_path)}")
    except Exception as e:
        logger.error(f"Failed to write gene filter log: {e}")
        print_warning(f"Failed to write gene filter log: {e}")


if __name__ == "__main__":
    main()
