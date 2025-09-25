import gzip
import logging
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

# Optional progress bar -----------------------------------------------------
from tusco_selector.optional_deps import HAS_TQDM, tqdm

# ----------------------------------------------------------------------------
# Module-level logger
# ----------------------------------------------------------------------------
logger = logging.getLogger(__name__)

# Public symbols intentionally exported for use across the package
__all__ = [
    # Parsing & mapping utilities
    "parse_gtf",
    "load_assembly_mapping",
    "load_junctions_bed",
    # GTF subsetting / reconstruction
    "extract_gene_id",
    "subset_gtf_by_gene_ids",
    "build_transcripts_dict_from_gtf",
    # Mapping subset / summaries
    "subset_mapping_file",
    "count_exons_in_gtf",
    # TUSCO deliverables
    "generate_tusco_tables",
]

# ----------------------------------------------------------------------------
# Generic helpers
# ----------------------------------------------------------------------------


def _open_maybe_gzip(path: Path | str):
    """Return a *text* handle regardless whether *path* is gzipped or not."""

    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")  # type: ignore[arg-type]
    return open(path, "r")


def _parse_attributes(attr_field: str) -> Dict[str, str]:
    """Return a mapping of GTF attributes from column 9 of a GTF line."""

    attrs: Dict[str, str] = {}
    # Attributes are separated by semicolons; each key/value pair uses space after key
    for raw in attr_field.strip().split(";"):
        raw = raw.strip()
        if not raw or " " not in raw:
            # Skip malformed attributes gracefully
            continue
        key, value = raw.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def _normalise_chr(raw_chr: str, mapping: Dict[str, str] | None = None) -> str:
    """Return a *consistent* chromosome label using *mapping* if provided.

    Rules
    -----
    1. If *raw_chr* matches a key in *mapping* → use mapped value.
    2. Else if it already starts with 'chr' / 'CHR' (case-insensitive) → ensure lower-case 'chr' prefix.
    3. Otherwise prefix with 'chr'.
    """

    if mapping and raw_chr in mapping:
        return mapping[raw_chr]

    if raw_chr.lower().startswith("chr"):
        return "chr" + raw_chr[3:]

    return "chr" + raw_chr


# ----------------------------------------------------------------------------
# Public parsing helpers  ----------------------------------------------------
# ----------------------------------------------------------------------------

GeneTranscripts = Dict[str, Set[str]]
TranscriptExons = Dict[str, List[Tuple[str, int, int]]]
TranscriptStrand = Dict[str, str]


def parse_gtf(
    gtf_path: Path | str,
    chr_mapping: Dict[str, str] | None = None,
) -> Tuple[GeneTranscripts, TranscriptExons, TranscriptStrand]:
    """Parse *gtf_path* and extract exon coordinates for every transcript.

    Parameters
    ----------
    gtf_path
        Path to the GTF (optionally gzipped).
    chr_mapping
        Optional mapping from RefSeq → UCSC chromosome labels (only required for
        human GRCh38).  If *None* an empty mapping is assumed.

    Returns
    -------
    gene_tx
        Mapping *gene_id → set(transcript_id)*.
    tx_exons
        Mapping *transcript_id → List[(chrom, start, end)]*.
    tx_strand
        Mapping *transcript_id → strand ('+'|'-')*.
    """

    chr_mapping = chr_mapping or {}

    gene_tx: GeneTranscripts = {}
    tx_exons: TranscriptExons = {}
    tx_strand: TranscriptStrand = {}

    try:
        # ── Two-pass approach: First count data lines -> progress bar length ──
        with _open_maybe_gzip(gtf_path) as handle:
            total_lines = sum(1 for line in handle if line and not line.startswith("#"))

        gene_tx = {}
        tx_exons = {}
        tx_strand = {}

        with _open_maybe_gzip(gtf_path) as handle:
            iterator = tqdm(
                handle,
                total=total_lines,
                desc=f"Parsing {Path(gtf_path).name}",
                unit="lines",
                disable=not HAS_TQDM,
            )
            for line in iterator:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue  # malformed – skip
                seqname, _source, feature, start, end, _score, strand, _frame, attrs = (
                    parts
                )
                if feature != "exon":
                    continue

                attr_map = _parse_attributes(attrs)
                gene_id = attr_map.get("gene_id")
                tx_id = attr_map.get("transcript_id")
                if not gene_id or not tx_id:
                    continue  # missing essential identifiers

                # Store exon coordinates
                chr_name = _normalise_chr(seqname, chr_mapping)
                tx_exons.setdefault(tx_id, []).append((chr_name, int(start), int(end)))
                tx_strand[tx_id] = strand
                gene_tx.setdefault(gene_id, set()).add(tx_id)

        logger.debug(
            "Extracted %d genes with %d transcripts from %s",
            len(gene_tx),
            len(tx_exons),
            Path(gtf_path).name,
        )
    except FileNotFoundError:
        logger.warning("GTF not found: %s", gtf_path)
    except Exception as exc:  # pragma: no cover – catch-all for robustness
        logger.error("Failed to parse %s: %s", gtf_path, exc)

    return gene_tx, tx_exons, tx_strand


# ----------------------------------------------------------------------------
# Assembly report helper (RefSeq ↔ UCSC chrom mapping) -----------------------
# ----------------------------------------------------------------------------


def load_assembly_mapping(report_file: Path | str | None) -> Dict[str, str]:
    """Parse an *assembly_report.txt* file and return RefSeq → UCSC mapping."""

    mapping: Dict[str, str] = {}
    if not report_file:
        return mapping
    report_file = Path(report_file)
    if not report_file.is_file():
        return mapping  # silently ignore – best effort

    with report_file.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue  # unexpected – skip
            refseq_accn = parts[6]
            ucsc_name = parts[9] if len(parts) > 9 else ""
            if ucsc_name:
                mapping[refseq_accn] = ucsc_name
    return mapping


# ----------------------------------------------------------------------------
# Splice-junction BED loader --------------------------------------------------
# ----------------------------------------------------------------------------

JunctionKey = Tuple[str, int, int, str]  # (chrom, start, end, strand)


def load_junctions_bed(
    junctions_bed_path: Path | str, species: str = ""
) -> Dict[JunctionKey, int]:
    """Load junction coverage from a *recount3*-style BED file.

    Only junctions with a canonical splice site (GT donor / AG acceptor, i.e. ``GT:AG``)
    are retained, unless the species is mmu.

    Parameters
    ----------
    junctions_bed_path
        Path to the junctions BED file.
    species
        Species identifier; if "mmu", splice site filtering is bypassed.

    Returns
    -------
    Dict[JunctionKey, int]
        Mapping ``(chrom, start, end, strand) → coverage`` for canonical junctions.
    """

    junctions: Dict[JunctionKey, int] = {}

    try:
        opener = gzip.open if str(junctions_bed_path).endswith(".gz") else open
        with opener(junctions_bed_path, "rt") as handle:  # type: ignore[arg-type]
            iterator = tqdm(handle, desc="Loading junction data", disable=not HAS_TQDM)
            for line in iterator:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue

                chrom, start, end, name_field, score, strand = parts[:6]

                # Skip splice site check for mouse
                if species.lower() != "mmu":
                    # The splice-site motif is stored as the last semicolon‑separated
                    # token in column 4, e.g. "…;GT:AG"
                    splice_site = name_field.split(";")[-1].upper()

                    # Keep only canonical GT/AG splice junctions
                    if splice_site != "GT:AG":
                        continue

                try:
                    cov = int(score)
                except ValueError:
                    cov = 1  # fallback

                junctions[(chrom, int(start), int(end), strand)] = cov

        log_msg = "Loaded %d junctions from %s"
        # For non-mouse species, we filtered to canonical GT:AG only
        if species.lower() != "mmu":
            log_msg = "Loaded %d canonical GT:AG junctions from %s"
        logger.info(log_msg, len(junctions), junctions_bed_path)
    except FileNotFoundError:
        logger.warning("Junction BED not found: %s", junctions_bed_path)
    except Exception as exc:
        logger.error("Error loading junction BED %s: %s", junctions_bed_path, exc)

    return junctions


# ----------------------------------------------------------------------------
# GTF helper functions for subset and reconstruction -------------------------
# ----------------------------------------------------------------------------


def extract_gene_id(attr_field: str, *, strip_version: bool = True) -> str | None:
    """Extract the *gene_id* value from a raw GTF *attr_field* (column 9).

    Parameters
    ----------
    attr_field
        The raw attribute string from a GTF line (column 9).
    strip_version
        Whether to remove a trailing Ensembl version suffix (".NN").

    Returns
    -------
    str | None
        The extracted gene ID or *None* if not present.
    """

    try:
        attrs = _parse_attributes(attr_field)
        gene_id = attrs.get("gene_id")
        if gene_id and strip_version:
            gene_id = gene_id.split(".")[0]
        return gene_id
    except Exception:
        # Be lenient – return None on any error
        return None


def subset_gtf_by_gene_ids(
    source_gtf: Path | str,
    gene_ids: set[str],
    dest_gtf: Path | str,
) -> int:
    """Write *dest_gtf* containing only records from *gene_ids* based on *source_gtf*.

    The function preserves header/comment lines (starting with "#") *as is* and
    copies all feature lines where the *gene_id* attribute (version-stripped)
    matches one of *gene_ids*.

    The destination file format (plain vs. gzipped) is inferred from the file
    extension: if it ends with ".gz" output is gzip-compressed.

    Returns the number of feature lines written (excluding header lines).
    """

    source_gtf = Path(source_gtf)
    dest_gtf = Path(dest_gtf)

    if dest_gtf.suffix == ".gz":
        opener_out = lambda p: gzip.open(p, "wt")  # type: ignore[arg-type]
    else:
        opener_out = lambda p: open(p, "w")  # type: ignore[arg-type]

    written = 0

    opener_in = gzip.open if str(source_gtf).endswith(".gz") else open  # type: ignore[arg-type]
    with opener_in(source_gtf, "rt") as ih, opener_out(dest_gtf) as oh:  # type: ignore[misc]
        for line in ih:
            if line.startswith("#"):
                oh.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue  # malformed – skip
            gene_id = extract_gene_id(parts[8])
            if gene_id and gene_id in gene_ids:
                oh.write(line)
                written += 1

    logger.info("Wrote %d records to subset GTF %s", written, dest_gtf)
    return written


def build_transcripts_dict_from_gtf(gtf_path: Path | str) -> Dict[str, Dict[str, Any]]:
    """Reconstruct a *transcripts* mapping from a (subset) GTF file.

    The returned structure mirrors the output of
    :func:`tusco_selector.pipeline.single_isoform.select_matched_single_isoforms`,
    i.e. a mapping ``gene_id → {transcript_id, exons, strand, chrom}``.
    """

    gene_tx, tx_exons, tx_strand = parse_gtf(gtf_path, chr_mapping=None)

    transcripts: Dict[str, Dict[str, Any]] = {}
    for gene_id, tx_ids in gene_tx.items():
        tx_id = next(iter(tx_ids))  # single-isoform by construction
        exons = sorted(tx_exons[tx_id], key=lambda x: (x[0], x[1]))
        transcripts[gene_id.split(".")[0]] = {
            "transcript_id": tx_id,
            "exons": exons,
            "strand": tx_strand.get(tx_id, "+"),
            "chrom": exons[0][0] if exons else "",
        }

    logger.info(
        "Reconstructed transcript mapping for %d genes from %s",
        len(transcripts),
        Path(gtf_path).name,
    )
    return transcripts


def subset_mapping_file(gene_ids, mapping_file, output_file, description):
    """
    Subset a mapping file (TSV) by filtering rows where the first column (Ensembl gene ID)
    is in the provided list of gene IDs.

    If a single Ensembl gene ID has multiple rows, they will be merged into one row
    with unique values in each column separated by commas.

    Args:
        gene_ids (list): List of Ensembl gene IDs to keep
        mapping_file (str): Path to the input mapping TSV file
        output_file (str): Path to save the filtered TSV file
        description (str): Description of the filter for the header

    Returns:
        str: Path to the output file
    """
    import csv
    import gzip
    from collections import defaultdict

    gene_ids_set = set(gene_ids)

    # Determine if input file is gzipped
    is_gzipped = mapping_file.endswith(".gz")
    opener = gzip.open if is_gzipped else open
    mode = "rt" if is_gzipped else "r"

    # Group rows by gene ID
    gene_rows = defaultdict(list)
    header = None

    with opener(mapping_file, mode) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)  # Try to get header if it exists

        # If there's no header, use the first row as data and create a default header
        if not header:
            f.seek(0)  # Reset file pointer
            reader = csv.reader(f, delimiter="\t")

        for row in reader:
            if row and row[0] in gene_ids_set:
                gene_id = row[0]
                gene_rows[gene_id].append(row)

    # Merge rows for each gene ID
    merged_rows = []
    merged_count = 0

    for gene_id, rows in gene_rows.items():
        if len(rows) == 1:
            # Only one row for this gene ID, no need to merge
            merged_rows.append(rows[0])
        else:
            # Multiple rows for this gene ID, merge them
            merged_row = [gene_id]  # Start with the gene ID

            # For each column (except the first one which is the gene ID)
            for col_idx in range(1, len(rows[0])):
                # Get unique values for this column across all rows
                unique_values = set()
                for row in rows:
                    if (
                        col_idx < len(row) and row[col_idx]
                    ):  # Check if column exists and not empty
                        unique_values.add(row[col_idx])

                # Join unique values with commas
                merged_row.append(",".join(sorted(unique_values)))

            merged_rows.append(merged_row)
            merged_count += 1

    # Write merged data to output file
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Write header comments
        f.write(f"# Gene IDs passing {description} filter\n")
        f.write(f"# Total: {len(merged_rows)}\n")
        if merged_count > 0:
            f.write(
                f"# Note: {merged_count} gene IDs had multiple rows that were merged\n"
            )

        # Write data rows (skip header)
        for row in merged_rows:
            writer.writerow(row)

    logger.info(
        "Saved %d merged gene entries to %s (%d were merged from duplicates)",
        len(merged_rows),
        output_file,
        merged_count,
    )

    return output_file


def count_exons_in_gtf(gtf_path: Path | str) -> Tuple[int, int, int]:
    """Parses a GTF file and counts total, single-exon, and multi-exon genes.

    Args:
        gtf_path: Path to the GTF file (can be gzipped).

    Returns:
        A tuple: (total_genes, single_exon_genes, multi_exon_genes).
    """
    gtf_path = Path(gtf_path)
    if not gtf_path.exists():
        logger.warning(f"GTF file not found for exon counting: {gtf_path}")
        return 0, 0, 0

    gene_tx, tx_exons, _ = parse_gtf(gtf_path, chr_mapping=None)

    total_genes = len(gene_tx)
    single_exon_genes = 0
    multi_exon_genes = 0

    # Create a map from transcript_id back to gene_id for easier lookup
    tx_to_gene = {}
    for gene_id, tx_ids in gene_tx.items():
        for tx_id in tx_ids:
            tx_to_gene[tx_id] = gene_id

    # Use gene_tx which has unique gene IDs as keys
    processed_genes = set()
    for gene_id, transcript_ids in gene_tx.items():
        if gene_id in processed_genes:
            continue  # Should not happen with typical GTFs but good practice

        # Since previous steps aim for single isoform, we usually expect one tx_id here.
        # However, we count exons based on *any* transcript associated with the gene found in the file.
        max_exons_for_gene = 0
        has_transcript_info = False
        for tx_id in transcript_ids:
            if tx_id in tx_exons:
                has_transcript_info = True
                num_exons = len(tx_exons[tx_id])
                if num_exons > max_exons_for_gene:
                    max_exons_for_gene = num_exons

        if has_transcript_info:  # Only count genes for which we found exon info
            if max_exons_for_gene == 1:
                single_exon_genes += 1
            elif max_exons_for_gene > 1:
                multi_exon_genes += 1
            # If max_exons_for_gene is 0 (e.g., gene feature without exons), it's not counted as single or multi.
            processed_genes.add(gene_id)
        else:
            logger.debug(
                f"Gene {gene_id} found but no corresponding exon features in the GTF."
            )

    # Recalculate total based on genes with actual exon counts
    counted_genes = single_exon_genes + multi_exon_genes
    if counted_genes != total_genes:
        logger.warning(
            f"Discrepancy in gene counts for {gtf_path}. Initial parse found {total_genes}, but only {counted_genes} had exon features."
        )
        # Decide how to report: Either return counted_genes or total_genes. Let's return the count of genes we could classify.
        total_genes = counted_genes

    logger.info(
        f"Counted exons in {gtf_path}: {total_genes} genes ({single_exon_genes} single, {multi_exon_genes} multi)"
    )
    return total_genes, single_exon_genes, multi_exon_genes


# ----------------------------------------------------------------------------
# TUSCO table generators ------------------------------------------------------
# ----------------------------------------------------------------------------


def _classify_genes_by_exons(
    gtf_path: Path | str,
) -> Tuple[Set[str], Set[str]]:
    """Return sets of single- and multi-exon gene IDs from a GTF.

    Gene IDs are returned without version suffixes.
    """

    gene_tx, tx_exons, _ = parse_gtf(gtf_path, chr_mapping=None)

    single: Set[str] = set()
    multi: Set[str] = set()

    for raw_gene_id, transcript_ids in gene_tx.items():
        gene_id = raw_gene_id.split(".")[0]
        max_exons = 0
        for tx_id in transcript_ids:
            if tx_id in tx_exons:
                max_exons = max(max_exons, len(tx_exons[tx_id]))
        if max_exons == 1:
            single.add(gene_id)
        elif max_exons > 1:
            multi.add(gene_id)

    return single, multi


def generate_tusco_tables(
    mapping_tsv_path: Path | str,
    final_gtf_path: Path | str,
    out_dir: Path | str | None = None,
    human: str = "hsa",
    mouse: str = "mmu",
    *,
    species_label: str | None = None,
) -> Tuple[Path, Path, Path]:
    """Generate consolidated TUSCO TSVs from final mapping and GTF.

    Creates three files in ``out_dir`` (defaults to the mapping file's dir):
      - If ``species_label`` is provided: ``tusco_{species_label}.tsv`` and
        corresponding ``..._multi_exon.tsv`` and ``..._single_exon.tsv``.
      - Otherwise (backwards compatibility): ``tusco_{human}{mouse}.tsv``,
        ``tusco_{human}{mouse}_multi_exon.tsv``, and
        ``tusco_{human}{mouse}_single_exon.tsv``.

    The rows are copied directly from ``mapping_tsv_path`` (data-only; lines
    beginning with "#" are ignored). Single/multi files are split by exon counts
    computed from ``final_gtf_path``.
    """

    import csv

    mapping_tsv_path = Path(mapping_tsv_path)
    final_gtf_path = Path(final_gtf_path)
    if out_dir is None:
        out_dir = mapping_tsv_path.parent
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Classify genes
    single_genes, multi_genes = _classify_genes_by_exons(final_gtf_path)
    all_genes = single_genes | multi_genes

    # Read mapping rows (skip comment lines)
    rows: List[List[str]] = []
    if not mapping_tsv_path.exists():
        logger.warning("Mapping TSV not found: %s", mapping_tsv_path)
    else:
        with mapping_tsv_path.open("r") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for raw in reader:
                if not raw:
                    continue
                if raw[0].startswith("#"):
                    continue
                rows.append(raw)

    # Filter rows by presence in GTF-derived gene set
    def _filter_by_gene_set(target: Set[str]) -> List[List[str]]:
        kept: List[List[str]] = []
        for row in rows:
            gene_id = row[0].split(".")[0] if row else ""
            if gene_id in target:
                kept.append(row)
        return kept

    rows_all = _filter_by_gene_set(all_genes) if rows else []
    rows_multi = _filter_by_gene_set(multi_genes) if rows else []
    rows_single = _filter_by_gene_set(single_genes) if rows else []

    # Output paths
    base = f"tusco_{species_label}" if species_label else f"tusco_{human}{mouse}"
    out_all = out_dir / f"{base}.tsv"
    out_multi = out_dir / f"{base}_multi_exon.tsv"
    out_single = out_dir / f"{base}_single_exon.tsv"

    def _write(path: Path, data_rows: List[List[str]], label: str) -> None:
        with path.open("w", newline="") as oh:
            writer = csv.writer(oh, delimiter="\t")
            oh.write(f"# Generated from {mapping_tsv_path.name} and {final_gtf_path.name}\n")
            oh.write(f"# {label} genes: {len(data_rows)}\n")
            for r in data_rows:
                writer.writerow(r)
        logger.info("Wrote %d rows to %s", len(data_rows), path)

    _write(out_all, rows_all, "All")
    _write(out_multi, rows_multi, "Multi-exon")
    _write(out_single, rows_single, "Single-exon")

    return out_all, out_multi, out_single
