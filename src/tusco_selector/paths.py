"""Path utilities and standardized file naming."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Tuple

from tusco_selector.logging_utils import get_logger

logger = get_logger(__name__)


def project_root() -> Path:
    """Return the absolute project root directory.
    
    This resolves to the repository root by walking up from this module's file.
    
    Returns:
        Path to project root directory
    """
    # src/tusco_selector/paths.py -> src -> project root
    return Path(__file__).parent.parent.parent.resolve()


def get_cache_dir(species: str = None) -> Path:
    """Get the cache directory for downloaded resources.
    
    Args:
        species: Optional species code to create species-specific subdirectory
        
    Returns:
        Path to cache directory
    """
    root = project_root()
    cache_dir = root / "cache"
    
    if species:
        cache_dir = cache_dir / species
    
    # Create if doesn't exist
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    return cache_dir


def get_output_dir(species: str = None, custom_dir: str = None) -> Path:
    """Get the output directory for pipeline results.
    
    Args:
        species: Species code (e.g., 'hsa', 'mmu', 'dre')
        custom_dir: Custom output directory path
        
    Returns:
        Path to output directory
    """
    if custom_dir:
        output_dir = Path(custom_dir)
    else:
        # Default to species-specific output directory
        root = project_root()
        if species:
            output_dir = root / f"{species}_output_gtfs"
        else:
            output_dir = root / "output_gtfs"
    
    # Create if doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    return output_dir


def get_step_file_names(step_num: int, slug: str) -> Tuple[str, str]:
    """Generate standard file names for a pipeline step.
    
    Args:
        step_num: Step number (1-based)
        slug: Step identifier slug (e.g., 'single_isoform', 'splice_tss_filtered')
        
    Returns:
        Tuple of (gtf_filename, mapping_filename)
    """
    gtf_filename = f"step{step_num}_{slug}.gtf.gz"
    mapping_filename = f"step{step_num}_{slug}_mapping.tsv"
    
    return gtf_filename, mapping_filename


def get_step_paths(output_dir: Path | str, step_num: int, slug: str) -> Tuple[Path, Path]:
    """Generate full paths for pipeline step output files.
    
    Args:
        output_dir: Output directory path
        step_num: Step number (1-based)
        slug: Step identifier slug
        
    Returns:
        Tuple of (gtf_path, mapping_path) as Path objects
    """
    output_dir = Path(output_dir)
    gtf_filename, mapping_filename = get_step_file_names(step_num, slug)
    
    return output_dir / gtf_filename, output_dir / mapping_filename


def get_final_output_names(species: str) -> Tuple[str, str, str, str]:
    """Get standardized final output file names for a species.
    
    Args:
        species: Species code
        
    Returns:
        Tuple of (gtf_name, tsv_name, single_exon_name, multi_exon_name)
    """
    if species == "hsa":
        prefix = "tusco_human"
    elif species == "mmu":
        prefix = "tusco_mouse"
    elif species == "dre":
        prefix = "tusco_zebrafish"
    else:
        prefix = f"tusco_{species}"
    
    return (
        f"{prefix}.gtf",
        f"{prefix}.tsv",
        f"{prefix}_single_exon.tsv",
        f"{prefix}_multi_exon.tsv",
    )


def get_tissue_output_names(species: str, tissue_name: str) -> Tuple[str, str, str]:
    """Get standardized tissue-specific output file names.
    
    Args:
        species: Species code
        tissue_name: Tissue name (will be sanitized)
        
    Returns:
        Tuple of (gtf_name, tsv_name, detailed_tsv_name)
    """
    # Sanitize tissue name for file naming
    safe_tissue = tissue_name.lower().replace(" ", "_").replace("/", "_")
    
    prefix = f"tusco_{species}_{safe_tissue}"
    
    return (
        f"{prefix}.gtf.gz",
        f"{prefix}.tsv",
        f"{prefix}_detailed.tsv",
    )


def ensure_dir(path: Path | str) -> Path:
    """Ensure a directory exists, creating it if necessary.
    
    Args:
        path: Directory path
        
    Returns:
        Path object for the directory
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_path(path: str | Path, base_dir: Path = None) -> Path:
    """Resolve a path, optionally relative to a base directory.
    
    Args:
        path: Path to resolve
        base_dir: Optional base directory for relative paths
        
    Returns:
        Resolved absolute Path
    """
    path = Path(path)
    
    if path.is_absolute():
        return path
    
    if base_dir:
        return (Path(base_dir) / path).resolve()
    
    return path.resolve()


__all__ = [
    'project_root',
    'get_cache_dir',
    'get_output_dir',
    'get_step_file_names',
    'get_step_paths',
    'get_final_output_names',
    'get_tissue_output_names',
    'ensure_dir',
    'resolve_path',
]