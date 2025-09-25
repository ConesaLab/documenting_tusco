"""Expression filtering helpers for universal and tissue-specific genes.

This module provides utilities to:
- Read GTEx ``.gct.gz`` matrices
- Read and pre-filter Bgee ``*_expr_simple.tsv.gz``
- Derive universally expressed genes and per-tissue high-expression sets
- Load housekeeping genes and map transcripts to genes
"""

import gzip
import time
import requests

# Optional logging
import logging
from pathlib import Path
from typing import Dict, Set, Tuple
from collections import Counter

import pandas as pd

logger = logging.getLogger(__name__)

# Optional dependencies  
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

# -----------------------------------------------------------------------------
# Expression-specific helpers (GCT handling, high-expression gene selection)
# -----------------------------------------------------------------------------

# GTEx tissue mapping (for TSV interest filtering)
TSV_FILES_OF_INTEREST: Set[str] = set()  # Can be customized as needed

def _extract_tissue_name(gct_path: Path) -> str:
    """Extract tissue name from GCT file path."""
    # Extract filename without extensions 
    stem = gct_path.name.replace('.gct.gz', '').replace('.gct', '')
    # Remove any prefix patterns if needed
    return stem

def read_gct_gz(file_path: Path) -> pd.DataFrame:
    """Read a GTEx .gct.gz file into a DataFrame.
    
    GCT format:
    - Line 1: #1.2 (version)
    - Line 2: nrows ncols (dimensions)
    - Line 3: column headers (Name, Description, sample1, sample2, ...)
    - Lines 4+: data rows
    """
    with gzip.open(file_path, 'rt') as f:
        # Skip version line
        version_line = f.readline().strip()
        if not version_line.startswith('#'):
            raise ValueError(f"Expected GCT version line, got: {version_line}")
        
        # Read dimensions
        _nrows, _ncols = map(int, f.readline().strip().split())
        
        # Read the rest as a DataFrame
        df = pd.read_csv(f, sep='\t', dtype={'Name': str})
        
    return df

def read_encode_expression_tsv(file_path: Path) -> pd.DataFrame:
    """Read an ENCODE expression TSV file."""
    # Placeholder function - not currently used for GTEx workflow
    return pd.read_csv(file_path, sep='\t', dtype={'Name': str})


# FANTOM5 (wide TPM matrix) helpers
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Bgee (long TSV format) helpers
# -----------------------------------------------------------------------------


def load_bgee_anatomical_mapping(mapping_file_path: str | Path) -> pd.DataFrame:
    """
    Loads Bgee anatomical mapping from the RNA-Seq libraries TSV file.

    Parameters
    ----------
    mapping_file_path
        Path to the Bgee 'Mus_musculus_RNA-Seq_libraries.tsv' (or similar) file.

    Returns
    -------
    pd.DataFrame
        DataFrame with "Anatomical entity ID" and "Anatomical entity name",
        with duplicates dropped to ensure one name per ID.
    """
    mapping_file_path = Path(mapping_file_path)
    logger.info("Loading Bgee anatomical mapping from: %s", mapping_file_path)
    try:
        # Assuming the file is a standard TSV
        df_map = pd.read_csv(
            mapping_file_path,
            sep="\t",
            usecols=["Anatomical entity ID", "Anatomical entity name"],
            dtype=str,
        )
        # Remove quotes from names if present, e.g., "leukocyte" -> leukocyte
        df_map["Anatomical entity name"] = df_map["Anatomical entity name"].str.replace(
            '"', ""
        )
        # Drop duplicates: if an ID has multiple names, keep the first one encountered.
        df_map.drop_duplicates(
            subset=["Anatomical entity ID"], keep="first", inplace=True
        )
        logger.info("Loaded %d unique anatomical ID to name mappings.", len(df_map))
    except Exception as e:
        logger.error("Error reading Bgee anatomical mapping file %s: %s", mapping_file_path, e)
        # Return an empty DataFrame with expected columns on error
        return pd.DataFrame(columns=["Anatomical entity ID", "Anatomical entity name"])
    return df_map


def read_bgee_expr_simple_tsv(
    file_path: str | Path,
    anatomical_mapping_df: pd.DataFrame,
    min_genes_per_tissue: int = 10000,
    qual_ok: Tuple[str, ...] = ("gold quality", "silver quality"),
) -> pd.DataFrame:
    """Read and pre-filter a Bgee ``*_expr_simple.tsv.gz`` file.

    Simple filtering logic:
    - Call quality in `qual_ok` (gold/silver quality)
    - Expression == "present"
    - Keep tissues with at least `min_genes_per_tissue` distinct genes
    """
    file_path = Path(file_path)
    logger.info("Reading Bgee expression file: %s", file_path)

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
        logger.error("Error reading Bgee file %s: %s", file_path, e)
        return pd.DataFrame()

    logger.info("Initial Bgee records: %d", len(df))

    # Simple filtering: good quality and filter expression == "present"
    df_filtered = df[df["Call quality"].isin(qual_ok) & (df["Expression"] == "present")]
    logger.info("Records after quality and expression filters: %d", len(df_filtered))

    if df_filtered.empty:
        logger.warning("No data left after Bgee quality filters.")
        return pd.DataFrame()

    # Merge with anatomical names (keep simple - drop if no mapping)
    df_with_names = pd.merge(
        df_filtered, anatomical_mapping_df, on="Anatomical entity ID", how="left"
    )
    logger.info("Records after anatomical name mapping (no drop): %d", len(df_with_names))

    if df_with_names.empty:
        logger.warning("No data left after anatomical name mapping.")
        return pd.DataFrame()

    # Count genes per tissue and keep only tissues with enough genes
    tissue_sizes = (
        df_with_names.groupby("Anatomical entity name")["Gene ID"]
        .nunique()
        .reset_index(name="n_genes")
    )

    good_tissue_names = tissue_sizes[tissue_sizes["n_genes"] >= min_genes_per_tissue][
        "Anatomical entity name"
    ].tolist()

    if not good_tissue_names:
        logger.warning("No tissues found with at least %d genes.", min_genes_per_tissue)
        return pd.DataFrame()

    logger.info(
        "Retaining %d tissues with >= %d genes.", len(good_tissue_names), min_genes_per_tissue
    )

    # Final filtering: keep only records from good tissues
    df_final = df_with_names[
        df_with_names["Anatomical entity name"].isin(good_tissue_names)
    ]
    logger.info("Final records after tissue size filtering: %d", len(df_final))

    return df_final


def get_universal_high_expression_genes_bgee_with_ids(
    file_path: str | Path,
    bgee_mapping_file_path: str | Path,
    *,
    min_genes_per_tissue: int = 10000,
    qual_ok: Tuple[str, ...] = ("gold quality", "silver quality"),
    prevalence_threshold: float = 1.0,
) -> Tuple[Set[str], Dict[Tuple[str, str], Set[str]]]:
    """Identify universally expressed genes and tissue-specific genes from Bgee.

    This variant returns tissue information as tuples of (anatomical_entity_id, tissue_name).

    Parameters
    ----------
    file_path
        Path to the Bgee ``*_expr_simple.tsv.gz`` file.
    bgee_mapping_file_path
        Path to the Bgee anatomical mapping file (e.g., 'Mus_musculus_RNA-Seq_libraries.tsv').
    min_genes_per_tissue
        Minimum number of distinct genes a tissue must have to be considered.
    qual_ok
        Tuple of "Call quality" values to accept (e.g., "gold quality").
    prevalence_threshold
        Minimum fraction of tissues where a gene must be expressed to be considered universal.
        Default 1.0 means gene must be present in all tissues.

    Returns
    -------
    Tuple[Set[str], Dict[Tuple[str, str], Set[str]]]
        - Set with universally expressed gene IDs (present in >= prevalence_threshold fraction of tissues).
        - Mapping ``(anatomical_entity_id, tissue_name) → set(gene_id)`` for tissue-specific genes.
    """
    anatomical_mapping_df = load_bgee_anatomical_mapping(bgee_mapping_file_path)
    if anatomical_mapping_df.empty:
        logger.error("Bgee anatomical mapping failed to load or is empty. Cannot proceed with Bgee name mapping.")
        return set(), {}

    df_good = read_bgee_expr_simple_tsv(
        file_path,
        anatomical_mapping_df=anatomical_mapping_df,
        min_genes_per_tissue=min_genes_per_tissue,
        qual_ok=qual_ok,
    )

    if df_good.empty:
        logger.warning("Bgee pre-processing returned empty DataFrame. No genes will be identified.")
        return set(), {}

    # Helper to drop version suffix from Ensembl IDs (if any, though Bgee usually doesn't have them for mouse)
    def _strip(gid: str) -> str:
        return gid.split(".")[0] if "." in gid else gid

    # D. Genes that appear in **all** remaining tissues
    # Count distinct tissues (by name) each gene appears in
    gene_tissue_counts = (
        df_good.groupby("Gene ID")["Anatomical entity name"]
        .nunique()
        .reset_index(name="n_tissues_with_gene")
    )

    n_good_tissues = df_good["Anatomical entity name"].nunique()
    if n_good_tissues == 0:  # Should be caught by df_good.empty earlier, but defensive
        logger.warning(
            "No good tissues (by name) available for Bgee universal gene calculation."
        )
        return set(), {}

    # Calculate minimum number of tissues required based on prevalence threshold
    min_tissues_required = max(1, int(n_good_tissues * prevalence_threshold))
    
    # Get genes that appear in at least the minimum number of tissues required
    universal_gene_ids: Set[str] = set(
        gene_tissue_counts[gene_tissue_counts["n_tissues_with_gene"] >= min_tissues_required][
            "Gene ID"
        ].apply(_strip)
    )

    logger.info(
        "Bgee universal genes (in >= %d/%d tissues, prevalence >= %.2f): %d",
        min_tissues_required,
        n_good_tissues,
        prevalence_threshold,
        len(universal_gene_ids),
    )

    # E. Tissue-specific genes (using (anatomical_entity_id, tissue_name) tuples as keys)
    high_exp_per_tissue: Dict[Tuple[str, str], Set[str]] = {}

    # Get unique anatomical entity ID and name combinations
    unique_tissues = df_good[
        ["Anatomical entity ID", "Anatomical entity name"]
    ].drop_duplicates()

    for _, row in unique_tissues.iterrows():
        anatomical_id = row["Anatomical entity ID"]
        tissue_name = row["Anatomical entity name"]

        # Get genes for this tissue
        tissue_genes = set(
            df_good[df_good["Anatomical entity ID"] == anatomical_id]["Gene ID"].apply(
                _strip
            )
        )

        # Use tuple of (anatomical_entity_id, tissue_name) as key
        high_exp_per_tissue[(anatomical_id, tissue_name)] = tissue_genes

    logger.info(
        "Bgee tissue-specific gene sets created for %d tissues (with anatomical IDs).",
        len(high_exp_per_tissue),
    )

    return universal_gene_ids, high_exp_per_tissue


def get_universal_high_expression_genes(
    gct_files: list[Path | str],
    *,
    prevalence_expression_cutoff: float = 0.1,
    median_expression_cutoff: float = 1,
    prevalence_threshold: float = 0.99,
    tsv_files_of_interest: Set[str] | None = None,
    universal_tissue_fraction: float = 1.0,
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """Identify *universally* and tissue‑specific expressed genes using prevalence
    and median expression criteria.

    A gene is considered expressed in a tissue if **both** of the following hold:

    • The fraction of samples whose TPM exceeds *expression_cutoff* is at least
      *prevalence_threshold*; and  
    • The **median** TPM across all samples in the tissue is at least
      *median_cutoff*.

    A gene is classified as *universal* if it satisfies both criteria in every
    tissue processed.

    Parameters
    ----------
    gct_files
        Sequence of paths to ``*.gct.gz`` files – each corresponding to one
        tissue.  Usually obtained from :func:`process_gct_files`.
    prevalence_expression_cutoff
        TPM threshold used for the prevalence test.
    median_expression_cutoff
        TPM threshold that the median expression must meet or exceed.
    prevalence_threshold
        Minimum fraction of samples that must be expressed (TPM >
        *prevalence_expression_cutoff*).
    tsv_files_of_interest
        Optional set with the **file names** (e.g. ``encode.heart.expression.tsv``)
        of ENCODE ``*.expression.tsv`` files that should be considered when
        computing the *universal* gene intersection.  Other ENCODE files will
        still be analysed for tissue-specific genes but will be ignored for the
        universal gene calculation.  If *None* (default) the predefined
        :data:`TSV_FILES_OF_INTEREST` constant is used.
    universal_tissue_fraction
        Minimum fraction of tissues that must satisfy the expression criteria for a gene to be called universal.

    Returns
    -------
    Tuple[Set[str], Dict[str, Set[str]]]
        1. Set with universally expressed gene IDs.
        2. Mapping ``tissue → set(gene_id)`` with the genes that meet the
           criteria in that tissue.
    """

    # Helper to remove version suffixes (e.g. ENSG00000123456.3 → ENSG00000123456)
    def _strip_version(gid: str) -> str:
        return gid.split(".")[0] if "." in gid else gid

    high_exp_per_tissue: Dict[str, Set[str]] = {}
    include_for_universal: Dict[str, bool] = {}  # tissue → bool flag
    gene_ids: list[str] | None = None

    # Use provided interest list or fall back to the global default
    if tsv_files_of_interest is None:
        tsv_files_of_interest = TSV_FILES_OF_INTEREST

    iterator = tqdm(gct_files, desc="Reading expression files", unit="file") if HAS_TQDM else gct_files

    for gct_path in iterator:
        gct_path = Path(gct_path)
        tissue = _extract_tissue_name(gct_path)

        # Select appropriate reader based on file suffix
        if gct_path.suffix == ".gz" and gct_path.name.endswith(".gct.gz"):
            df = read_gct_gz(gct_path)
            include_in_universal = True  # Always include GTEx / GCT files
        else:
            df = read_encode_expression_tsv(gct_path)
            # Include ENCODE TSVs only if explicitly listed
            include_in_universal = gct_path.name in tsv_files_of_interest

        if gene_ids is None:
            gene_ids = df["Name"].astype(str).tolist()

        expr = df.iloc[:, 2:].astype(float)   # numeric expression matrix
        n_samples = expr.shape[1]

        prevalence = (expr > prevalence_expression_cutoff).sum(axis=1) / n_samples
        median_expr = expr.median(axis=1)

        genes_pass = (prevalence >= prevalence_threshold) & (median_expr >= median_expression_cutoff)
        genes_passing = {
            _strip_version(gene_ids[i]) for i, keep in enumerate(genes_pass) if keep
        }
        high_exp_per_tissue[tissue] = genes_passing
        include_for_universal[tissue] = include_in_universal

        logger.info(
            "%s: %d/%d genes pass prevalence ≥ %.2f (TPM > %.1f) and median ≥ %.1f TPM (n=%d)",
            tissue, len(genes_passing), len(gene_ids),
            prevalence_threshold, prevalence_expression_cutoff, median_expression_cutoff, n_samples,
        )

    # Determine universal housekeeping genes as the intersection across tissues
    if high_exp_per_tissue:
        # Select only tissues flagged for inclusion
        tissues_for_universal = [t for t, ok in include_for_universal.items() if ok]

        if not tissues_for_universal:
            logger.warning(
                "None of the provided tissues qualified for universal gene calculation; falling back to all tissues"
            )
            tissues_for_universal = list(high_exp_per_tissue.keys())

        universal_sets = [high_exp_per_tissue[t] for t in tissues_for_universal]
        tissue_count = len(tissues_for_universal)

        # Decide strategy: strict intersection (default) vs. relaxed fraction-based inclusion
        # The relaxed strategy is activated only when a TSV interest list is provided *and*
        # `universal_tissue_fraction` is < 1.0 (e.g. 0.9 for 90 % of tissues).
        if universal_tissue_fraction < 1.0 and tsv_files_of_interest:
            gene_counter: Counter[str] = Counter()
            for gene_set in universal_sets:
                gene_counter.update(gene_set)

            # Build the "table": gene → #tissues where it passes the criteria
            # and keep genes present in at least the requested fraction of tissues.
            universal_genes = {
                gene for gene, cnt in gene_counter.items()
                if (cnt / tissue_count) >= universal_tissue_fraction
            }

            logger.info(
                "Universal set: %d genes present in ≥ %.2f of %d selected tissues (%d total tissues processed)",
                len(universal_genes), universal_tissue_fraction, tissue_count, len(high_exp_per_tissue),
            )
        else:
            # Original behaviour: strict intersection across tissues
            if tissue_count == 1:
                universal_genes = universal_sets[0].copy()
            else:
                universal_genes = set.intersection(*universal_sets)

            logger.info(
                "Universal set: %d genes across %d selected tissues (%d total tissues processed)",
                len(universal_genes), tissue_count, len(high_exp_per_tissue),
            )
    else:
        universal_genes = set()
        logger.warning("No genes passed in any tissue")

    return universal_genes, high_exp_per_tissue


def load_housekeeping_genes(file_path: str | Path) -> Set[str]:
    """Load housekeeping genes from CSV file and convert transcript IDs to gene IDs.
    
    The housekeeping genes file is expected to be in a CSV format with semicolons as delimiters,
    where the first column contains Ensembl transcript IDs (e.g., ENST00000123456).
    
    Parameters
    ----------
    file_path
        Path to CSV file containing housekeeping genes with transcript IDs.
        
    Returns
    -------
    Set[str]
        Set of Ensembl gene IDs for housekeeping genes.
    """
    # Try to read the file with various possible delimiters
    try:
        df = pd.read_csv(file_path, sep=';')
        logger.info("Read housekeeping genes file with semicolon delimiter")
    except Exception:
        try:
            df = pd.read_csv(file_path)
            logger.info("Read housekeeping genes file with comma delimiter")
        except Exception:
            try:
                df = pd.read_csv(file_path, sep='\t')
                logger.info("Read housekeeping genes file with tab delimiter")
            except Exception as e:
                logger.error(f"Failed to read housekeeping genes file: {e}")
                return set()
    
    # Get the first column which should contain Ensembl transcript IDs
    if len(df.columns) == 0:
        logger.error("Housekeeping genes file has no columns")
        return set()
    
    transcript_ids = df.iloc[:, 0].tolist()
    
    # Log the number of transcript IDs found
    logger.info("Found %d transcript IDs in housekeeping genes file", len(transcript_ids))
    
    # Extract gene IDs directly from transcript IDs (for Ensembl format)
    gene_ids = transcript_ids_to_gene_ids(transcript_ids)
    
    logger.info("Extracted %d gene IDs from housekeeping genes file", len(gene_ids))
    
    return gene_ids


def transcript_ids_to_gene_ids(transcript_ids: list[str]) -> Set[str]:
    """Convert Ensembl transcript IDs to gene IDs using batched POST calls (1000 IDs per batch)."""
    gene_ids: Set[str] = set()
    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    session = requests.Session()
    batch_size = 1000

    # Use tqdm for progress if available
    iterator = tqdm(range(0, len(transcript_ids), batch_size), desc="Converting transcript IDs") if HAS_TQDM else range(0, len(transcript_ids), batch_size)

    for i in iterator:
        batch = transcript_ids[i:i + batch_size]
        try:
            response = session.post(
                url,
                headers=headers,
                params={"object_type": "transcript"},
                json={"ids": batch},
            )
            if response.status_code == 429:
                retry_after = int(response.headers.get("Retry-After", 1))
                logger.warning("Rate limited by Ensembl API. Waiting %d seconds.", retry_after)
                time.sleep(retry_after)
                response = session.post(
                    url,
                    headers=headers,
                    params={"object_type": "transcript"},
                    json={"ids": batch},
                )

            response.raise_for_status()
            data = response.json()
            for tid, entry in data.items():
                if not entry:
                    continue
                parent = entry.get("Parent") or entry.get("parent")
                if not parent:
                    continue
                gene_id = parent.split(".")[0]
                gene_ids.add(gene_id)
        except Exception as e:
            logger.error("Error fetching transcript batch starting at index %d: %s", i, e)
            continue # Continue with the next batch
        # Optional small delay between batches
        time.sleep(0.1)
        
    return gene_ids


# Re-export for upstream modules (CLI) ---------------------------------------
__all__ = [
    "get_universal_high_expression_genes_bgee_with_ids",
    "get_universal_high_expression_genes",
    "load_housekeeping_genes",
]
