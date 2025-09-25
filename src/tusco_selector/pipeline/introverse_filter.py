from __future__ import annotations

"""Introverse filter helper functions.

This module provides utilities to load the *IntroVerse* mis-splicing catalogue
(`introverse_10.csv`) and use the genomic coordinates to identify overlapping
Ensembl gene IDs via the Ensembl REST API.  The resulting gene-ID set can be
used to remove genes that show widespread mis-splicing (>10 % individuals)
from the TUSCO candidate set.

The heavy lifting (HTTP requests) is performed lazily and cached in-memory so
repeat invocations during a single run do not issue duplicate requests.
"""

import json
import logging
import os
from pathlib import Path
from typing import Dict, Set, Tuple

# Third-party dependencies are optional at *import* time so that test discovery
# does not require the full stack.  We import them lazily inside helper
# functions.

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helper – cache Ensembl look-ups to avoid duplicate HTTP requests
# ---------------------------------------------------------------------------

# Simple in-memory cache:  ``location -> {gene_id, …}``
_LOCATION_CACHE: Dict[str, Set[str]] = {}


# ---------------------------------------------------------------------------
# Core helper – map genomic location to Ensembl gene IDs
# ---------------------------------------------------------------------------


def _ensembl_genes_for_location(location: str) -> Set[str]:
    """Return *set* of Ensembl gene IDs overlapping *location*.

    The *location* string must follow the pattern ``"chr:start-end:strand"`` as
    used by the IntroVerse CSV.  Chromosome names are passed to the Ensembl
    endpoint without the leading ``chr`` as the REST server expects bare
    chromosome labels (e.g. ``1`` not ``chr1``).
    """

    # Use cached result if available ------------------------------------------------
    cached = _LOCATION_CACHE.get(location)
    if cached is not None:
        return cached

    # Import *requests* lazily to keep dependency footprint light ---------------
    import requests  # pylint: disable=import-error

    chrom, positions, _strand = location.split(":")
    # Strip optional leading "chr"
    chrom = chrom.removeprefix("chr")
    start, end = positions.split("-")

    server = "https://rest.ensembl.org"
    ext = f"/overlap/region/human/{chrom}:{start}-{end}?feature=gene"
    headers = {"Content-Type": "application/json"}

    try:
        resp = requests.get(f"{server}{ext}", headers=headers, timeout=30)
        resp.raise_for_status()
    except Exception as exc:  # pragma: no cover – network errors
        logger.warning("Ensembl REST query failed for %s – %s", location, exc)
        return set()

    genes_json = resp.json()
    gene_ids: Set[str] = {
        g["id"].split(".")[0] for g in genes_json if g.get("feature_type") == "gene"
    }

    # Cache and return -----------------------------------------------------------
    _LOCATION_CACHE[location] = gene_ids
    return gene_ids


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def read_csv_and_get_ensembl_ids(csv_path: str | Path) -> Set[str]:  # noqa: D401
    """Return *set* of Ensembl gene IDs present in *csv_path*.

    The function loads the IntroVerse CSV and queries the Ensembl REST API for
    every *unique* genomic location listed in the ``ID`` column.  Results are
    cached to minimise redundant HTTP requests.
    """

    csv_path = Path(csv_path)
    if not csv_path.is_file():
        logger.warning("IntroVerse CSV not found: %s", csv_path)
        return set()

    # Lazy import of pandas + tqdm -------------------------------------------------
    import pandas as pd  # pylint: disable=import-error

    from tusco_selector.optional_deps import tqdm

    def _progress(iterable, **kwargs):  # noqa: D401
        return tqdm(iterable, **kwargs)

    logger.info("Loading IntroVerse CSV (%s)", csv_path.name)
    df = pd.read_csv(csv_path)
    if "ID" not in df.columns:
        logger.warning("IntroVerse CSV missing 'ID' column – %s", csv_path)
        return set()

    unique_locations = df["ID"].dropna().unique()
    logger.info("%d unique genomic locations in IntroVerse CSV", len(unique_locations))

    gene_ids: Set[str] = set()
    for loc in _progress(
        unique_locations, desc="Querying Ensembl for overlapping genes"
    ):
        # Ensure consistent formatting (strip quotes, whitespace)
        loc = str(loc).replace('"', "").strip()
        gene_ids.update(_ensembl_genes_for_location(loc))

    logger.info(
        "IntroVerse mapping identified %d unique Ensembl gene IDs", len(gene_ids)
    )
    return gene_ids


def apply_introverse_filter(
    transcripts: Dict[str, dict], introverse_genes: Set[str]
) -> Tuple[Dict[str, dict], dict, Set[str]]:
    """Return *transcripts* excluding genes present in *introverse_genes*."""

    if not introverse_genes:
        logger.info("IntroVerse gene set empty – no genes removed")
        return (
            transcripts,
            {
                "before": {
                    "total": len(transcripts),
                    "single_exon": 0,
                    "multi_exon": 0,
                },
                "after": {"total": len(transcripts), "single_exon": 0, "multi_exon": 0},
                "removed": {"total": 0, "single_exon": 0, "multi_exon": 0},
            },
            set(),
        )

    # Count single-exon vs multi-exon genes before filtering
    single_exon_before = 0
    multi_exon_before = 0

    for gene_id, info in transcripts.items():
        exons = info.get("exons", [])
        if len(exons) == 1:
            single_exon_before += 1
        else:
            multi_exon_before += 1

    # Apply filtering
    kept = {
        gid: info for gid, info in transcripts.items() if gid not in introverse_genes
    }

    # Count after filtering
    single_exon_after = 0
    multi_exon_after = 0

    # Count removed genes by type
    single_exon_removed = 0
    multi_exon_removed = 0

    for gene_id, info in kept.items():
        exons = info.get("exons", [])
        if len(exons) == 1:
            single_exon_after += 1
        else:
            multi_exon_after += 1

    single_exon_removed = single_exon_before - single_exon_after
    multi_exon_removed = multi_exon_before - multi_exon_after
    total_removed = len(transcripts) - len(kept)

    logger.info(
        "Before IntroVerse filtering: %d genes total (%d single-exon, %d multi-exon)",
        len(transcripts),
        single_exon_before,
        multi_exon_before,
    )
    logger.info(
        "After IntroVerse filtering: %d genes remain (%d single-exon, %d multi-exon)",
        len(kept),
        single_exon_after,
        multi_exon_after,
    )
    logger.info(
        "IntroVerse filter removed %d genes (%d single-exon, %d multi-exon)",
        total_removed,
        single_exon_removed,
        multi_exon_removed,
    )

    removed_ids = set(transcripts.keys()) - set(kept.keys())

    stats = {
        "before": {
            "total": len(transcripts),
            "single_exon": single_exon_before,
            "multi_exon": multi_exon_before,
        },
        "after": {
            "total": len(kept),
            "single_exon": single_exon_after,
            "multi_exon": multi_exon_after,
        },
        "removed": {
            "total": total_removed,
            "single_exon": single_exon_removed,
            "multi_exon": multi_exon_removed,
        },
    }
    return kept, stats, removed_ids
