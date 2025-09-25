"""
Step 3: Remove alternative splice and TSS transcripts.
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

# Setup module logger
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helper functions – splice junction handling
# ---------------------------------------------------------------------------

JunctionKey = Tuple[str, int, int, str]  # (chrom, start, end, strand)


def _get_annotated_splice_junctions(
    exons_sorted: List[Tuple[str, int, int]], strand: str
) -> List[JunctionKey]:
    """Return list of annotated splice junctions for *exons_sorted*.

    Junction coordinates follow recount3 style (0-based *start*, 1-based *end*).
    We simply take *exon_end* .. *next_exon_start-1* on the same chromosome.
    """

    annotated: List[JunctionKey] = []
    if len(exons_sorted) < 2:
        return annotated

    chrom = exons_sorted[0][0]

    for i in range(len(exons_sorted) - 1):
        _chr, exon_end = exons_sorted[i][0], exons_sorted[i][2]
        next_chr, next_start = exons_sorted[i + 1][0], exons_sorted[i + 1][1]
        if _chr != next_chr:
            # Unexpected: exons on different chromosomes – skip junction
            continue
        junc_start = exon_end  # 0-based inclusive (end of upstream exon)
        junc_end = next_start - 1  # inclusive, recount3 convention
        if junc_start < junc_end:
            annotated.append((chrom, junc_start, junc_end, strand))
    return annotated


def _average_coverage(
    annotated: List[JunctionKey], bed_cov: dict[JunctionKey, int]
) -> float:
    """Return average coverage across *annotated* junctions using *bed_cov*."""

    if not annotated:
        return 0.0
    values = [bed_cov.get(j, 0) for j in annotated]
    return sum(values) / len(values) if values else 0.0


# ---------------------------------------------------------------------------
# Fast junction lookup helpers
# ---------------------------------------------------------------------------


def _build_junction_index(
    bed_cov: dict[JunctionKey, int]
) -> dict[Tuple[str, str], list[Tuple[int, int, int]]]:
    """Return mapping *(chrom, strand) → sorted list of (start, end, cov)* for quick region queries."""
    index: dict[Tuple[str, str], list[Tuple[int, int, int]]] = defaultdict(list)
    for (chrom, s, e, strand), cov in bed_cov.items():
        index[(chrom, strand)].append((s, e, cov))
    # Sort each list by start coordinate so we can break early when *start* exceeds region end
    for lst in index.values():
        lst.sort(key=lambda x: x[0])
    return index


def _has_high_coverage_novel_junction(
    annotated: List[JunctionKey],
    junction_index: dict[Tuple[str, str], list[Tuple[int, int, int]]],
    gene_region: Tuple[int, int],
    coverage_cut: float,
    min_length: int = 0,
) -> Tuple[bool, Dict[JunctionKey, int]]:
    """Return ``True`` if *junction_index* contains a junction within *gene_region*
    that is *not* in *annotated* and has coverage >= *coverage_cut*.
    Also returns a dictionary of the novel junctions found and their coverage.
    The lookup is limited to the (chrom,strand) pairs present in *annotated*
    for speed.

    min_length : int
        Minimum junction length (bp) required to be considered; 0 disables length filtering.
    """
    gene_start, gene_end = gene_region
    found_novel = False
    novel_junctions_found = {}  # Store novel junctions { (chrom, s, e, strand): cov }
    if not annotated or coverage_cut <= 0:
        return False, novel_junctions_found

    annotated_set = set(annotated)
    chrom_strand_pairs = {(c, s) for (c, _s, _e, s) in annotated}

    for chrom, strand in chrom_strand_pairs:
        for s, e, cov in junction_index.get((chrom, strand), []):
            if s > gene_end:
                break  # lists are sorted by start – nothing further can overlap
            # Require the novel junction to be completely inside the gene region
            if s < gene_start or e > gene_end or cov < coverage_cut:
                continue
            if (e - s + 1) < min_length:
                continue  # junction too short
            junction_key = (chrom, s, e, strand)
            if junction_key in annotated_set:
                continue  # annotated junction – ignore

            # Found a high-coverage novel junction
            found_novel = True
            novel_junctions_found[junction_key] = cov
            # Continue checking other junctions in case multiple novel ones exist
            # We return all novel junctions exceeding the threshold.
            # The calling function will decide based on whether *any* novel junction was found.

    return found_novel, novel_junctions_found


# ---------------------------------------------------------------------------
# Additional helper: highest-coverage novel junction
# ---------------------------------------------------------------------------


def _highest_novel_junction_coverage(
    annotated: List[JunctionKey],
    junction_index: dict[Tuple[str, str], list[Tuple[int, int, int]]],
    gene_region: Tuple[int, int],
    min_length: int = 0,
) -> Tuple[int, int, int]:
    """Return a triple (max_cov, max_start, max_end) for the novel junction with the highest coverage
    inside gene_region, considering only junctions with length ≥ min_length. If none, return (0, -1, -1).
    """
    gene_start, gene_end = gene_region
    if not annotated:
        return 0, -1, -1

    annotated_set = set(annotated)
    chrom_strand_pairs = {(c, s) for (c, _s, _e, s) in annotated}

    max_cov = 0
    max_start = -1
    max_end = -1
    for chrom, strand in chrom_strand_pairs:
        for s, e, cov in junction_index.get((chrom, strand), []):
            if s > gene_end:
                break  # lists are sorted by start – nothing further can overlap
            # Only consider junctions fully contained within the gene region
            if s < gene_start or e > gene_end:
                continue
            if (e - s + 1) < min_length:
                continue  # skip short junctions
            if (chrom, s, e, strand) in annotated_set:
                continue  # annotated junction – ignore
            if cov > max_cov:
                max_cov = cov
                max_start = s
                max_end = e
    return max_cov, max_start, max_end


# ---------------------------------------------------------------------------
# Helper functions – TSS handling
# ---------------------------------------------------------------------------


def _load_tss_bed(tss_bed: str | Path) -> pd.DataFrame:
    """Utility to load a BED file with refTSS data and return dataframe."""

    df = pd.read_csv(
        tss_bed,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "refTSS_ID", "score", "strand"],
    )
    # Ensure chromosome naming consistent (chrX)
    df["chrom"] = df["chrom"].apply(
        lambda c: c if str(c).startswith("chr") else f"chr{c}"
    )
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["start", "end"])
    return df


def _normalize_chromosome_label(label: str) -> str:
    """Return chromosome label with 'chr' prefix to match refTSS BED.

    Ensures consistency between gene annotation (which may be '1') and refTSS
    (which is normalized to 'chr1' in :func:`_load_tss_bed`).
    """
    s = str(label)
    return s if s.startswith("chr") else f"chr{s}"


# ---------------------------------------------------------------------------
# Core filtering functions
# ---------------------------------------------------------------------------


def check_splice_junctions(
    transcripts: Dict[str, dict],
    junctions: dict[JunctionKey, int] | None,
    novel_threshold: float = 0.01,
    min_novel_length: int = 0,
) -> Tuple[Dict[str, dict], Dict[str, str]]:
    """Filter *transcripts* by removing multi-exon genes with high-coverage novel junctions.

    Parameters
    ----------
    transcripts
        Mapping *gene_id → info* as produced by :pyfunc:`single_isoform.select_matched_single_isoforms`.
        Each info dict must contain keys *exons* (List[Tuple[chrom,start,end]]) and *strand*.
    junctions
        Dictionary mapping *(chrom, start, end, strand)* to *coverage*.  If *None* → skip filtering.
    novel_threshold
        Novel junction coverage threshold relative to average annotated coverage (default 1 %).
    min_novel_length : int
        Minimum length (bp) a non-annotated junction must have to be considered (default 0).

    Returns
    -------
    Tuple[Dict[str, dict], Dict[str, str]]
        1. Dictionary of transcripts that passed the filter.
        2. Dictionary mapping removed gene IDs to their filter reason.
    """

    if junctions is None:
        logger.warning(
            "No junctions coverage supplied – skipping splice-junction filtering"
        )
        return transcripts, {}

    keep: Dict[str, dict] = {}
    # Store removed gene IDs and the reason
    removed_genes_with_reason: Dict[str, str] = {}

    # Build an index once for fast (chrom,strand) region queries
    junction_index = _build_junction_index(junctions)

    for gene_id, info in transcripts.items():
        exons = info.get("exons", [])
        strand = info.get("strand", "+")
        if len(exons) <= 1:
            keep[gene_id] = info  # single-exon – handled later by TSS
            continue

        # Multi-exon – compute annotated junctions and coverage (sort by start)
        exons_sorted = sorted(exons, key=lambda x: (x[1], x[2]))
        annotated = _get_annotated_splice_junctions(exons_sorted, strand)
        avg_cov = _average_coverage(annotated, junctions)
        coverage_cut = avg_cov * novel_threshold if avg_cov > 0 else 0
        # If coverage_cut is effectively zero (e.g., avg_cov is tiny), set a minimum absolute cutoff
        # to avoid filtering based on negligible novel coverage. Let's say 1 read minimum.
        effective_coverage_cut = max(coverage_cut, 1.0) if avg_cov > 0 else 1.0

        logger.debug(
            "Gene %s: %d annotated junctions, avg_cov=%.3f, novel_threshold=%.2f, effective_coverage_cut=%.3f",
            gene_id,
            len(annotated),
            avg_cov,
            novel_threshold,
            effective_coverage_cut,  # Use effective cut for logging
        )

        gene_start = min(e[1] for e in exons)
        gene_end = max(e[2] for e in exons)

        # Use the updated _has_high_coverage_novel_junction
        has_novel, novel_juncs = _has_high_coverage_novel_junction(
            annotated,
            junction_index,
            (gene_start, gene_end),
            effective_coverage_cut,
            min_novel_length,
        )

        if has_novel:
            # Log details of the first high-coverage novel junction found (for simplicity)
            first_novel_key = next(iter(novel_juncs.keys()))
            novel_cov = novel_juncs[first_novel_key]
            logger.debug(
                "Gene %s removed due to high-coverage novel junction "
                "(e.g., %s cov=%d, cut=%.3f)",
                gene_id,
                str(first_novel_key),
                novel_cov,
                effective_coverage_cut,
            )
            removed_genes_with_reason[gene_id] = "Novel Splice Junction"
        else:
            max_novel_cov, max_novel_start, max_novel_end = (
                _highest_novel_junction_coverage(
                    annotated, junction_index, (gene_start, gene_end), min_novel_length
                )
            )
            logger.debug(
                "Gene %s kept after splice-junction check (max_novel_cov=%d, cut=%.3f)",
                gene_id,
                max_novel_cov,
                effective_coverage_cut,
            )
            keep[gene_id] = info

    logger.info(
        "Splice-junction filtering removed %d of %d multi-exon genes",
        len(removed_genes_with_reason),
        len(transcripts),
    )
    return keep, removed_genes_with_reason


def check_tss(
    transcripts: Dict[str, dict],
    tss_df: pd.DataFrame | None,
    *,
    tss_region_bp: int | None = None,
) -> Tuple[Dict[str, dict], Dict[str, str]]:
    """Filter transcripts that have evidence for alternative or inconsistent TSS.

    Behaviour (applies to both single- and multi-exon genes):
    - Determine the transcript TSS exon and TSS coordinate:
      - Single-exon: the sole exon region; TSS is exon start ("+") or exon end ("-").
      - Multi-exon: use the first exon for "+" strand or the last exon for "-" strand;
        TSS is the exon boundary nearest the transcript start (start for "+", end for "-").
    - Find refTSS entries overlapping that TSS exon region (BED-style overlap).
      - If >1 overlaps: remove as "Alternative TSS (Multiple Overlaps)".
      - If exactly 1 overlap: require the TSS coordinate to fall inside that refTSS interval;
        otherwise remove as "TSS Mismatch".
      - If no overlaps: keep.

    Returns
    -------
    Tuple[Dict[str, dict], Dict[str, str]]
        1. Dictionary of transcripts that passed the filter.
        2. Dictionary mapping removed gene IDs to the reason ("Alternative TSS" or "TSS Mismatch").
    """

    if tss_df is None:
        logger.warning("No TSS dataframe supplied – skipping TSS filtering")
        return transcripts, {}

    keep: Dict[str, dict] = {}
    removed_genes_with_reason: Dict[str, str] = {}

    for gene_id, info in transcripts.items():
        exons: List[Tuple[str, int, int]] = info.get("exons", [])
        if not exons:
            # No exon information – cannot determine TSS; keep as-is
            keep[gene_id] = info
            continue

        strand = info.get("strand", "+")
        chrom = info.get("chrom", exons[0][0] if exons else "")
        chrom = _normalize_chromosome_label(chrom)

        # Identify the TSS exon
        exons_sorted = sorted(exons, key=lambda x: (x[1], x[2]))
        if strand == "+":
            tss_exon = exons_sorted[0]
            tss_pos = tss_exon[1]  # 1-based start
        else:
            tss_exon = exons_sorted[-1]
            tss_pos = tss_exon[2]  # 1-based end

        exon_start = tss_exon[1]
        exon_end = tss_exon[2]

        # Convert to 0-based half-open for comparison (BED style)
        # GTF/GFF exon coordinates are typically 1-based inclusive. BED is 0-based half-open.
        exon_start_0 = exon_start - 1
        exon_end_0 = exon_end
        tss_pos_0 = tss_pos - 1

        logger.debug(
            "Gene %s: TSS exon [%d, %d) on %s strand; TSS pos (1-based) = %d",
            gene_id,
            exon_start_0,
            exon_end_0,
            strand,
            tss_pos,
        )

        # Two modes:
        #  - Window mode (preferred): count CAGE peaks within ±tss_region_bp around TSS.
        #  - Legacy mode: use TSS exon overlap with refTSS and require exact containment when single overlap.
        if tss_region_bp and tss_region_bp > 0:
            window_start_0 = max(0, tss_pos_0 - tss_region_bp)
            window_end_0 = tss_pos_0 + tss_region_bp + 1  # half-open end
            overlaps = tss_df[
                (tss_df["chrom"] == chrom)
                & (tss_df["strand"] == strand)
                & (tss_df["start"] < window_end_0)
                & (tss_df["end"] > window_start_0)
            ]
            logger.debug(
                "Gene %s: CAGE peaks within ±%dbp window [%d, %d) = %d",
                gene_id,
                tss_region_bp,
                window_start_0,
                window_end_0,
                len(overlaps),
            )

            if len(overlaps) > 1:
                removed_genes_with_reason[gene_id] = (
                    f"Alternative TSS (>1 CAGE peak within ±{tss_region_bp}bp)"
                )
                continue
            else:
                keep[gene_id] = info
                continue
        else:
            overlaps = tss_df[
                (tss_df["chrom"] == chrom)
                & (tss_df["strand"] == strand)
                & (tss_df["start"] < exon_end_0)
                & (tss_df["end"] > exon_start_0)
            ]
            logger.debug(
                "Gene %s: overlapping TSS entries (exon-level) = %d", gene_id, len(overlaps)
            )
            if len(overlaps) > 1:
                logger.debug(
                    "Gene %s removed due to multiple TSS overlaps (%d overlaps)",
                    gene_id,
                    len(overlaps),
                )
                removed_genes_with_reason[gene_id] = "Alternative TSS (Multiple Overlaps)"
                continue
            if len(overlaps) == 1:
                ref_tss = overlaps.iloc[0]
                if ref_tss["start"] <= tss_pos_0 < ref_tss["end"]:
                    keep[gene_id] = info
                else:
                    logger.debug(
                        "Gene %s removed: TSS position %d (0-based: %d) not within refTSS [%d, %d)",
                        gene_id,
                        tss_pos,
                        tss_pos_0,
                        ref_tss["start"],
                        ref_tss["end"],
                    )
                    removed_genes_with_reason[gene_id] = "TSS Mismatch"
                continue
            # No overlapping TSS entries
            keep[gene_id] = info

    logger.info("TSS filtering removed %d genes", len(removed_genes_with_reason))
    return keep, removed_genes_with_reason


# ---------------------------------------------------------------------------
# Combined convenience function
# ---------------------------------------------------------------------------


def check_splice_junctions_tss(
    transcripts: Dict[str, dict],
    *,
    junctions: dict[JunctionKey, int] | None = None,
    tss_df: pd.DataFrame | None = None,
    novel_threshold: float = 0.01,
    min_novel_length: int = 0,
    tss_scope: str = "both",
    tss_region_bp: int | None = None,
):
    """Run TSS checks (configurable scope), then splice-junction checks for multi-exon genes.

    Order of operations:
    1. Apply :func:`check_tss` to selected transcripts depending on ``tss_scope``:
       - "single": only single-exon
       - "multi": only multi-exon
       - "both": all transcripts
    2. From the TSS-passing set, apply :func:`check_splice_junctions` to multi-exon genes.
    3. Merge results and return kept/removed sets with reasons.

    Parameters
    ----------
    min_novel_length : int
        Minimum length (bp) a non-annotated junction must have to be considered (default 0).

    Returns
    -------
    Tuple[Dict[str, dict], Dict[str, str]]
        1. Dictionary of transcripts that passed the filters.
        2. Dictionary mapping removed gene IDs to their filter reason.
    """

    logger.info(
        "Starting combined TSS + splice‑junction filtering (TSS scope: %s)", tss_scope
    )

    # Split upfront
    single_group = {gid: info for gid, info in transcripts.items() if len(info.get("exons", [])) == 1}
    multi_group = {gid: info for gid, info in transcripts.items() if len(info.get("exons", [])) > 1}

    kept_after_tss: Dict[str, dict] = {}
    removed_tss: Dict[str, str] = {}

    # Run TSS checks based on scope
    scope = (tss_scope or "both").lower()
    if scope not in {"single", "multi", "both"}:
        scope = "both"

    if scope in {"single", "both"}:
        # First check single-exon genes without using the window mode
        kept_s, removed_s = check_tss(single_group, tss_df, tss_region_bp=0)
        logger.info("[Exon-overlap mode] TSS filtering retained %d of %d single-exon genes", len(kept_s), len(single_group))
        kept_after_tss.update(kept_s)
        removed_tss.update(removed_s)
        
        # Then check multi-exon genes using the window mode
        kept_s, removed_s = check_tss(single_group, tss_df, tss_region_bp=tss_region_bp)
        logger.info("[Window mode] TSS filtering retained %d of %d single-exon genes", len(kept_s), len(single_group))
        kept_after_tss.update(kept_s)
        removed_tss.update(removed_s)
    else:
        kept_after_tss.update(single_group)

    if scope in {"multi", "both"}:
        kept_m, removed_m = check_tss(multi_group, tss_df, tss_region_bp=tss_region_bp)
        logger.info("[Window mode] TSS filtering retained %d of %d multi-exon genes", len(kept_m), len(multi_group))
        kept_after_tss.update(kept_m)
        removed_tss.update(removed_m)
    else:
        kept_after_tss.update(multi_group)

    # Now split the TSS-passing set into single- and multi-exon
    single_kept_after_tss: Dict[str, dict] = {}
    multi_candidates: Dict[str, dict] = {}
    for gid, info in kept_after_tss.items():
        if len(info.get("exons", [])) > 1:
            multi_candidates[gid] = info
        else:
            single_kept_after_tss[gid] = info

    # Run splice-junction checks on multi-exon genes that passed TSS
    filtered_multi, removed_multi = check_splice_junctions(
        multi_candidates, junctions, novel_threshold, min_novel_length
    )

    # Merge results
    combined_kept = {**single_kept_after_tss, **filtered_multi}
    combined_removed = {**removed_tss, **removed_multi}

    logger.info(
        "Filtering retained %d of %d transcripts (after TSS + splice checks)",
        len(combined_kept),
        len(transcripts),
    )

    logger.info("Combined filtering removed %d genes", len(combined_removed))

    return combined_kept, combined_removed
