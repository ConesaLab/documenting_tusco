"""Bgee expression data processing."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Dict, Set, Tuple

import pandas as pd

from tusco_selector.logging_utils import get_logger

logger = get_logger(__name__)


def load_anatomical_mapping(mapping_file_path: str | Path) -> pd.DataFrame:
    """Load anatomical entity ID to name mapping from Bgee TSV file.
    
    Expected columns in the mapping file:
    - Anatomical entity ID (e.g., "UBERON:0000001")
    - Anatomical entity name (e.g., "anatomical entity")
    
    Args:
        mapping_file_path: Path to Bgee anatomical mapping TSV file
        
    Returns:
        DataFrame with anatomical entity mappings
    """
    mapping_path = Path(mapping_file_path)
    
    if not mapping_path.exists():
        logger.warning(f"Anatomical mapping file not found: {mapping_path}")
        return pd.DataFrame()
    
    logger.info(f"Loading anatomical mapping from: {mapping_path}")
    
    try:
        # Read the mapping file
        df = pd.read_csv(
            mapping_path,
            sep='\t',
            dtype=str,
            usecols=['Anatomical entity ID', 'Anatomical entity name'],
        )
        
        # Remove duplicates if any
        df = df.drop_duplicates(subset=['Anatomical entity ID'])
        
        logger.info(f"Loaded {len(df)} anatomical entity mappings")
        return df
        
    except Exception as e:
        logger.error(f"Error loading anatomical mapping: {e}")
        return pd.DataFrame()


def read_expr_simple(
    file_path: str | Path,
    anatomical_mapping_df: pd.DataFrame,
    min_genes_per_tissue: int = 10000,
    qual_ok: Tuple[str, ...] = ("gold quality", "silver quality"),
    species: str = None,
    gtex_mapping_df: pd.DataFrame = None,
    include_all_tissues_for_tissue_specific: bool = False,
) -> pd.DataFrame:
    """Read and pre-filter a Bgee *_expr_simple.tsv.gz file.
    
    Simple filtering logic:
    - Call quality in qual_ok (gold/silver quality)
    - Expression == "present"
    - For universal genes: Keep tissues with at least min_genes_per_tissue distinct genes
    - For tissue-specific genes: Species-specific logic
    
    Args:
        file_path: Path to the Bgee expression file
        anatomical_mapping_df: DataFrame with anatomical entity ID to name mappings
        min_genes_per_tissue: Minimum genes per tissue for universal gene selection
        qual_ok: Acceptable quality levels
        species: Species code (e.g., 'hsa', 'mmu')
        gtex_mapping_df: GTEx tissue mapping for human tissues
        include_all_tissues_for_tissue_specific: Include all tissues for tissue-specific selection
        
    Returns:
        Filtered DataFrame with expression data
    """
    file_path = Path(file_path)
    logger.info(f"Reading Bgee expression file: {file_path}")
    
    try:
        with gzip.open(file_path, "rt") as handle:
            df = pd.read_csv(
                handle,
                sep="\t",
                usecols=[
                    "Gene ID",
                    "Anatomical entity ID",
                    "Expression",
                    "Call quality",
                ],
                dtype={
                    "Gene ID": str,
                    "Anatomical entity ID": str,
                    "Expression": str,
                    "Call quality": str,
                },
                engine="python",
            )
    except Exception as e:
        logger.error(f"Error reading Bgee file {file_path}: {e}")
        return pd.DataFrame()
    
    logger.info(f"Initial Bgee records: {len(df)}")
    
    # Simple filtering: good quality and expression == "present"
    df_filtered = df[
        df["Call quality"].isin(qual_ok) & (df["Expression"] == "present")
    ]
    logger.info(f"Records after quality and expression filters: {len(df_filtered)}")
    
    if df_filtered.empty:
        logger.warning("No data left after Bgee quality filters.")
        return pd.DataFrame()
    
    # Merge with anatomical names
    df_with_names = pd.merge(
        df_filtered, anatomical_mapping_df, on="Anatomical entity ID", how="left"
    )
    logger.info(f"Records after anatomical name mapping: {len(df_with_names)}")
    
    if df_with_names.empty:
        logger.warning("No data left after anatomical name mapping.")
        return pd.DataFrame()
    
    # Count genes per tissue
    tissue_sizes = (
        df_with_names.groupby("Anatomical entity name")["Gene ID"]
        .nunique()
        .reset_index(name="n_genes")
    )
    
    # Apply species-specific tissue selection logic
    if species == "hsa" and gtex_mapping_df is not None and not gtex_mapping_df.empty:
        # For humans: Include all GTEx tissues regardless of gene count
        gtex_uberon_codes = set(
            gtex_mapping_df["uberon_code"].dropna().str.replace("_", ":")
        )
        
        logger.info(
            f"GTEx mapping: {len(gtex_mapping_df)} tissues mapped to "
            f"{len(gtex_uberon_codes)} unique UBERON codes"
        )
        
        # Get tissues matching GTEx UBERON codes
        gtex_tissues = df_with_names[
            df_with_names["Anatomical entity ID"].isin(gtex_uberon_codes)
        ]["Anatomical entity name"].unique()
        
        # Combine criteria: GTEx tissues + high-gene-count tissues
        high_gene_tissues = tissue_sizes[tissue_sizes["n_genes"] >= min_genes_per_tissue][
            "Anatomical entity name"
        ].tolist()
        
        valid_tissues = set(gtex_tissues) | set(high_gene_tissues)
        logger.info(
            f"Including {len(valid_tissues)} tissues "
            f"({len(gtex_tissues)} GTEx + {len(high_gene_tissues)} high-gene)"
        )
        
    elif include_all_tissues_for_tissue_specific:
        # Include all tissues for tissue-specific genes
        valid_tissues = set(tissue_sizes["Anatomical entity name"])
        logger.info(f"Including all {len(valid_tissues)} tissues for tissue-specific analysis")
        
    else:
        # Default: only tissues with sufficient genes
        valid_tissues = set(
            tissue_sizes[tissue_sizes["n_genes"] >= min_genes_per_tissue][
                "Anatomical entity name"
            ]
        )
        logger.info(
            f"Keeping {len(valid_tissues)} tissues with >= {min_genes_per_tissue} genes"
        )
    
    # Filter to valid tissues
    df_final = df_with_names[df_with_names["Anatomical entity name"].isin(valid_tissues)]
    logger.info(f"Final records after tissue filtering: {len(df_final)}")
    
    return df_final


def get_universal_high_expression_genes_bgee_with_ids(
    file_path: str | Path,
    bgee_mapping_file_path: str | Path,
    *,
    prevalence_ubiquitous_cutoff: float = 0.9,
    min_tissues_for_high_expression: int = 10,
    min_genes_per_tissue: int = 10000,
    qual_ok: Tuple[str, ...] = ("gold quality", "silver quality"),
    species: str = None,
    gtex_mapping_df: pd.DataFrame = None,
    filter_tissues_for_tissue_specific: bool = True,
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """Identify universally high-expressed genes from Bgee data.
    
    A gene is considered universally high-expressed if it is expressed
    in >= prevalence_ubiquitous_cutoff fraction of tissues.
    
    Args:
        file_path: Path to Bgee expression file
        bgee_mapping_file_path: Path to anatomical mapping file
        prevalence_ubiquitous_cutoff: Minimum fraction of tissues for universal expression
        min_tissues_for_high_expression: Minimum absolute number of tissues
        min_genes_per_tissue: Minimum genes per tissue to include it
        qual_ok: Acceptable quality levels
        species: Species code
        gtex_mapping_df: GTEx tissue mapping for humans
        filter_tissues_for_tissue_specific: Whether to filter tissues for tissue-specific analysis
        
    Returns:
        Tuple of (universal_genes, high_exp_per_tissue_with_ids)
    """
    # Load anatomical mapping
    anatomical_mapping_df = load_anatomical_mapping(bgee_mapping_file_path)
    
    if anatomical_mapping_df.empty:
        logger.error("Failed to load anatomical mapping")
        return set(), {}
    
    # Read and filter expression data
    df = read_expr_simple(
        file_path=file_path,
        anatomical_mapping_df=anatomical_mapping_df,
        min_genes_per_tissue=min_genes_per_tissue,
        qual_ok=qual_ok,
        species=species,
        gtex_mapping_df=gtex_mapping_df,
        include_all_tissues_for_tissue_specific=not filter_tissues_for_tissue_specific,
    )
    
    if df.empty:
        logger.warning("No valid Bgee data after filtering")
        return set(), {}
    
    # Get unique tissues
    all_tissues = df["Anatomical entity name"].unique()
    n_tissues = len(all_tissues)
    
    # Calculate minimum tissues required
    min_tissues_required = max(
        min_tissues_for_high_expression,
        int(n_tissues * prevalence_ubiquitous_cutoff)
    )
    
    logger.info(
        f"Requiring presence in {min_tissues_required}/{n_tissues} tissues "
        f"for universal expression"
    )
    
    # Count tissue occurrences per gene
    gene_tissue_counts = (
        df.groupby("Gene ID")["Anatomical entity name"]
        .nunique()
        .reset_index(name="n_tissues")
    )
    
    # Select universally expressed genes
    universal_genes = set(
        gene_tissue_counts[gene_tissue_counts["n_tissues"] >= min_tissues_required]["Gene ID"]
    )
    
    logger.info(f"Found {len(universal_genes)} universally high-expressed genes")
    
    # Get per-tissue gene sets with anatomical IDs
    high_exp_per_tissue_with_ids = {}
    
    for tissue_name in all_tissues:
        tissue_df = df[df["Anatomical entity name"] == tissue_name]
        tissue_genes = set(tissue_df["Gene ID"].unique())
        
        # Get the anatomical entity ID for this tissue
        anatomical_id = tissue_df["Anatomical entity ID"].iloc[0] if not tissue_df.empty else None
        
        # Store with both name and ID
        high_exp_per_tissue_with_ids[tissue_name] = {
            "genes": tissue_genes,
            "anatomical_id": anatomical_id,
        }
    
    # For backward compatibility, also return just the gene sets
    high_exp_per_tissue = {
        tissue: data["genes"] for tissue, data in high_exp_per_tissue_with_ids.items()
    }
    
    return universal_genes, high_exp_per_tissue


__all__ = [
    'load_anatomical_mapping',
    'read_expr_simple',
    'get_universal_high_expression_genes_bgee_with_ids',
]