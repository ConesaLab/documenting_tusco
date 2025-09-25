"""
Step 4: Apply manual filtering rules.
"""

import logging
from pathlib import Path
from typing import Dict, Set, Tuple

logger = logging.getLogger(__name__)


def load_manual_filter(filter_path: Path) -> Set[str]:
    """Return set of gene IDs listed in *filter_path*.

    The file may contain one gene ID per line; leading/trailing whitespace and
    empty/comment lines (#) are ignored.  Version suffixes (".NN") are removed
    for consistency with pipeline gene IDs.
    """

    genes: Set[str] = set()
    if not filter_path.is_file():
        logger.warning("Manual filter file %s does not exist", filter_path)
        return genes

    with filter_path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genes.add(line.split(".")[0])

    logger.info("Loaded %d gene IDs from manual filter %s", len(genes), filter_path)
    return genes


def apply_manual_filter(
    transcripts: Dict[str, dict], filter_genes: Set[str]
) -> Tuple[Dict[str, dict], Set[str]]:
    """Return *transcripts* excluding *filter_genes*.

    Parameters
    ----------
    transcripts
        Current mapping ``gene_id → info``.
    filter_genes
        Set of gene IDs (version-stripped) to remove.

    Returns
    -------
    Tuple[Dict[str, dict], Set[str]]
        1. Dictionary of transcripts that passed the filter.
        2. Set of gene IDs that were removed by this filter.
    """

    if not filter_genes:
        logger.info("Manual filter set empty – no genes removed")
        return transcripts, set()

    kept: Dict[str, dict] = {
        gid: info for gid, info in transcripts.items() if gid not in filter_genes
    }
    removed_count = len(transcripts) - len(kept)
    removed_ids = filter_genes & set(transcripts.keys())
    logger.info(
        "Manual filter removed %d genes (remaining %d)", removed_count, len(kept)
    )
    return kept, removed_ids
