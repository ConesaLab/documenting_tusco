"""GTEx expression data processing."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

from tusco_selector.logging_utils import get_logger

logger = get_logger(__name__)


def read_gct_gz(file_path: Path) -> pd.DataFrame:
    """Read a GTEx .gct.gz file into a DataFrame.
    
    GCT format:
    - Line 1: version (e.g., "#1.2")
    - Line 2: dimensions (e.g., "57820   17382")
    - Line 3: header with "Name", "Description", followed by sample IDs
    - Rest: data rows with gene ID, description, followed by expression values
    
    Args:
        file_path: Path to .gct.gz file
        
    Returns:
        DataFrame with genes as rows and samples as columns
    """
    file_path = Path(file_path)
    logger.info(f"Reading GTEx GCT file: {file_path}")
    
    with gzip.open(file_path, 'rt') as f:
        # Skip version line
        version = f.readline().strip()
        logger.debug(f"GCT version: {version}")
        
        # Read dimensions
        dims = f.readline().strip().split('\t')
        n_genes, n_samples = int(dims[0]), int(dims[1])
        logger.info(f"Expected dimensions: {n_genes} genes x {n_samples} samples")
        
        # Read the rest into DataFrame
        df = pd.read_csv(f, sep='\t', index_col=0)
        
        # Drop the Description column if present
        if 'Description' in df.columns:
            df = df.drop(columns=['Description'])
        
        logger.info(f"Loaded expression matrix: {df.shape[0]} genes x {df.shape[1]} samples")
        
        return df


def load_gtex_tissue_mapping(gtex_mapping_file_path: str | Path) -> pd.DataFrame:
    """Load GTEx to UBERON tissue mapping from TSV file.
    
    Expected format:
    - Columns: tissue_name, uberon_code
    - Tab-separated values
    
    Args:
        gtex_mapping_file_path: Path to GTEx tissue mapping file
        
    Returns:
        DataFrame with tissue mapping
    """
    mapping_path = Path(gtex_mapping_file_path)
    
    if not mapping_path.exists():
        logger.warning(f"GTEx tissue mapping file not found: {mapping_path}")
        return pd.DataFrame()
    
    logger.info(f"Loading GTEx tissue mapping from: {mapping_path}")
    
    try:
        df = pd.read_csv(
            mapping_path,
            sep='\t',
            names=['tissue_name', 'uberon_code'],
            dtype=str,
        )
        logger.info(f"Loaded {len(df)} GTEx tissue mappings")
        return df
    except Exception as e:
        logger.error(f"Error loading GTEx tissue mapping: {e}")
        return pd.DataFrame()


def get_universal_high_expression_genes(
    gct_files: list[Path | str],
    *,
    prevalence_expression_cutoff: float = 0.1,
    prevalence_ubiquitous_cutoff: float = 0.9,
    min_tissues_for_high_expression: int = 3,
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """Identify universally high-expressed genes from GTEx data.
    
    A gene is considered:
    - High-expressed in a tissue if median TPM >= prevalence_expression_cutoff
    - Universally high-expressed if present in >= prevalence_ubiquitous_cutoff of tissues
    
    Args:
        gct_files: List of GTEx GCT file paths
        prevalence_expression_cutoff: Minimum median TPM for high expression
        prevalence_ubiquitous_cutoff: Minimum fraction of tissues for universal expression
        min_tissues_for_high_expression: Minimum tissues required
        
    Returns:
        Tuple of (universal_genes, high_exp_per_tissue)
    """
    if not gct_files:
        logger.warning("No GTEx files provided")
        return set(), {}
    
    all_tissue_data = {}
    tissue_names = []
    
    # Load all GTEx files
    for gct_file in gct_files:
        gct_path = Path(gct_file)
        
        if not gct_path.exists():
            logger.warning(f"GTEx file not found: {gct_path}")
            continue
        
        # Extract tissue name from filename
        # Expected format: gtex_v8_rnaseq_<tissue_name>_tpm.gct.gz
        tissue_name = gct_path.stem.replace('.gct', '').replace('gtex_v8_rnaseq_', '').replace('_tpm', '')
        tissue_names.append(tissue_name)
        
        logger.info(f"Processing GTEx tissue: {tissue_name}")
        
        try:
            df = read_gct_gz(gct_path)
            
            # Calculate median expression per gene
            median_expr = df.median(axis=1)
            
            # Find high-expressed genes
            high_exp_genes = set(median_expr[median_expr >= prevalence_expression_cutoff].index)
            all_tissue_data[tissue_name] = high_exp_genes
            
            logger.info(f"  {tissue_name}: {len(high_exp_genes)} high-expressed genes")
            
        except Exception as e:
            logger.error(f"Error processing {gct_path}: {e}")
            continue
    
    if not all_tissue_data:
        logger.warning("No valid GTEx data loaded")
        return set(), {}
    
    # Find universally high-expressed genes
    total_tissues = len(all_tissue_data)
    min_tissues_required = max(
        min_tissues_for_high_expression,
        int(total_tissues * prevalence_ubiquitous_cutoff)
    )
    
    logger.info(f"Requiring presence in {min_tissues_required}/{total_tissues} tissues for universal expression")
    
    # Count tissue occurrences for each gene
    gene_tissue_counts = {}
    for tissue_genes in all_tissue_data.values():
        for gene in tissue_genes:
            gene_tissue_counts[gene] = gene_tissue_counts.get(gene, 0) + 1
    
    # Select universally expressed genes
    universal_genes = {
        gene for gene, count in gene_tissue_counts.items()
        if count >= min_tissues_required
    }
    
    logger.info(f"Found {len(universal_genes)} universally high-expressed genes")
    
    # Also return per-tissue high expression for tissue-specific analysis
    high_exp_per_tissue = {
        tissue: genes for tissue, genes in all_tissue_data.items()
    }
    
    return universal_genes, high_exp_per_tissue


__all__ = [
    'read_gct_gz',
    'load_gtex_tissue_mapping',
    'get_universal_high_expression_genes',
]