import sys
import os
import random
import argparse
import re
from tqdm import tqdm
from collections import defaultdict
import logging

# make random picks reproducible when needed (debugging)
random.seed(42)

# Constants for canonical splice sites
CANONICAL_DONOR = 'GT'
CANONICAL_ACCEPTOR = 'AG'


# Reverse‑complements for minus‑strand genes
REVERSE_ACCEPTOR = 'CT'    # complement of AG
REVERSE_DONOR = 'AC'       # complement of GT

# Default window size for splice site search
DEFAULT_WINDOW = 100  # bases


def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate fake GTF annotations for selected genes')
    parser.add_argument('--gtf', required=True,
                        help='Input reference GTF file (RefSeq or GENCODE)')
    parser.add_argument('--target_genes', required=True,
                        help='List of target gene names (one per line) or TSV file')
    parser.add_argument('--input_format', choices=['auto', 'simple', 'tsv'], default='auto',
                        help='Format of the input genes file: "simple" (one gene per line), '
                             '"tsv" (TSV with gene symbols in column 3), or "auto" to detect (default)')
    parser.add_argument('--symbol_column', type=int, default=3,
                        help='Column number (1-based) containing gene symbols in TSV format (default: 3)')
    parser.add_argument('--genome', required=True,
                        help='Reference genome FASTA (for splice site sequences)')
    parser.add_argument('--output_fake', required=True,
                        help='Output GTF containing only fake target transcripts')
    parser.add_argument('--output_full', required=True,
                        help='Full GTF with fake target genes replaced')
    parser.add_argument('--log', required=True,
                        help='Detailed modification log file')
    return parser.parse_args()


def load_target_gene_list(target_file, input_format='auto', symbol_column=3):
    """
    Read target gene symbols and return them in *upper‑case* so we can
    match case‑insensitively against different GTF conventions.
    
    Supports two formats:
    1. Simple text file with one gene name per line
    2. TSV file with gene symbols in the specified column (default: column 3)
       If a gene symbol is empty, falls back to gene ID in the next column
    
    Parameters
    ----------
    target_file : str
        Path to the input file
    input_format : str
        'auto' (detect format), 'simple' (one gene per line), or 'tsv' 
        (TSV with gene symbols in specified column)
    symbol_column : int
        1-based column index containing gene symbols in TSV format (default: 3)
    """
    # Convert to 0-based index for processing
    symbol_idx = symbol_column - 1
    fallback_idx = symbol_idx + 1
    
    gene_symbols = set()
    with open(target_file) as fh:
        lines = [line.strip() for line in fh if line.strip() and not line.startswith('#')]
        
        # Determine the file format
        is_tsv = False
        if input_format == 'auto':
            # Auto-detect format by looking for tabs in the first non-comment line
            is_tsv = lines and '\t' in lines[0]
            logging.info(f"Auto-detected format: {'TSV' if is_tsv else 'simple list'}")
        else:
            is_tsv = input_format == 'tsv'
            logging.info(f"Using specified format: {input_format}")
            
        if is_tsv:
            # TSV format - extract gene symbols from the specified column
            logging.info(f"Using column {symbol_column} for gene symbols")
            for line in lines:
                columns = line.split('\t')
                if len(columns) > symbol_idx:  # Ensure we have enough columns
                    gene_symbol = columns[symbol_idx].strip()
                    
                    # If gene symbol is empty but we have a fallback column, use that instead
                    if not gene_symbol and len(columns) > fallback_idx:
                        gene_id = columns[fallback_idx].strip()
                        if gene_id:
                            logging.debug(f"Using fallback ID '{gene_id}' in place of empty gene symbol")
                            gene_symbol = gene_id
                    
                    if gene_symbol:  # Only add non-empty gene symbols
                        gene_symbols.add(gene_symbol.upper())
        else:
            # Simple list format - one gene per line
            gene_symbols = set(line.upper() for line in lines)
    
    logging.info(f"Loaded {len(gene_symbols)} unique gene symbols")
    return gene_symbols


def setup_logging(log_file):
    """
    Configure root logger to write DEBUG information to *log_file* and
    INFO‑level messages to the console, giving us a detailed trace that
    can be inspected later.
    """
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(message)s')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
    logging.getLogger('').addHandler(console)


def is_canonical(seq, donor_pos, acceptor_pos):
    donor = seq[donor_pos:donor_pos+2].upper()
    acceptor = seq[acceptor_pos-2:acceptor_pos].upper()
    return donor == CANONICAL_DONOR and acceptor == CANONICAL_ACCEPTOR


def adjust_site(seq, pos, motif, max_shift, orientation, min_bound=None, max_bound=None):
    """
    Shift *pos* along *orientation* ( +1 = downstream, -1 = upstream on the
    reference strand) up to *max_shift* bp, limited to [min_bound, max_bound]
    if those are given (relative coordinates to the same origin as *pos*).
    Returns the new coordinate where *motif* matches or *None* if none found.
    """
    if 0 <= pos <= len(seq) - len(motif) and seq[pos:pos+len(motif)].upper() == motif:
        return pos
    for delta in range(1, max_shift + 1):
        cand = pos + orientation * delta
        if cand < 0 or cand + len(motif) > len(seq):
            break
        if min_bound is not None and cand < min_bound:
            continue
        if max_bound is not None and cand + len(motif) - 1 > max_bound:
            continue
        if seq[cand:cand+len(motif)].upper() == motif:
            return cand
    return None


def motif_at(seq, pos):
    """Return the two–base motif at *pos* (upper‑case) or '??' if out of range."""
    if 0 <= pos <= len(seq) - 2:
        return seq[pos:pos+2].upper()
    return '??'

def label_site(motif, canonical_motif):
    canonical = motif == canonical_motif
    return 'CANONICAL' if canonical else 'NON‑CANONICAL'

# Helper to find canonical splice site
def find_canonical_site(seq, rel_pos, motif, orientation,
                        intron_min, intron_max, min_shift):
    """
    Search unidirectionally along *orientation* (±1) beginning *min_shift*
    bases away from *rel_pos* until the end of the intron, looking for *motif*.
    Returns the new *relative* coordinate if found, else None.
    """
    step = orientation
    for delta in range(min_shift, intron_max - intron_min - 1):
        cand = rel_pos + step * delta
        if cand < intron_min or cand + len(motif) - 1 > intron_max:
            break
        if seq[cand:cand + len(motif)].upper() == motif:
            return cand
    return None


def process_single_exon(chrom, start, end, strand, seq, gene_start, min_intron=50):
    # pick donor and acceptor positions within exon such that intron >= min_intron
    L = end - start + 1
    if L < min_intron + 10:
        raise ValueError('Exon too short for fake intron')
    # try up to 100 random splits until we get canonical motifs
    for attempt in range(100):
        p1 = random.randint(start + 5, end - min_intron - 5)
        p2 = random.randint(p1 + min_intron, end - 5)
        donor_motif    = CANONICAL_DONOR    if strand == '+' else REVERSE_DONOR
        acceptor_motif = CANONICAL_ACCEPTOR if strand == '+' else REVERSE_ACCEPTOR
        donor_orientation    = +1 if strand == '+' else -1
        acceptor_orientation = -1 if strand == '+' else +1
        d_rel = p1 - gene_start
        a_rel = p2 - gene_start
        d = adjust_site(seq, d_rel, donor_motif, 20, donor_orientation)
        a = adjust_site(seq, a_rel, acceptor_motif, 20, acceptor_orientation)
        if d is not None and a is not None:
            donor = gene_start + d
            acceptor = gene_start + a
            if abs(donor - p1) >= 5 and abs(acceptor - p2) >= 5:
                break
    else:
        # fallback to original positions if canonical motifs truly absent
        donor, acceptor = p1, p2
    # return two exons, skipping the 2-nt splice motifs
    if strand == '+':
        # donor motif at intron start, acceptor motif at intron end
        ex1_end    = donor - 1
        ex2_start  = acceptor + 2
    else:
        # on minus strand, acceptor motif at intron start, donor motif at intron end
        ex1_end    = acceptor - 1
        ex2_start  = donor + 2

    ex1_start = start
    ex2_end   = end
    return [
        (chrom, ex1_start, ex1_end, strand),
        (chrom, ex2_start, ex2_end, strand)
    ], donor, acceptor


def process_multi_exon(exons, seq, gene_start, min_shift=5, window=DEFAULT_WINDOW):
    """
    Adjust splice junctions for a multi‑exon transcript.

    Parameters
    ----------
    exons : list[tuple]
        [(chrom, start, end, strand), ...] **ascending genomic order**.
    seq : str
        Genomic sequence of the whole gene (same orientation as reference).
    gene_start : int
        Absolute start coordinate of *seq* on the reference genome (1‑based).
    Returns
    -------
    new_exons : list[tuple]
        Same structure as *exons* but with modified coordinates.
    modified   : list[tuple]
        (orig_donor, new_donor, orig_acceptor, new_acceptor,
         orig_donor_motif, new_donor_motif,
         orig_acceptor_motif, new_acceptor_motif)
        for every intron.
    """
    strand = exons[0][3]
    donor_motif    = CANONICAL_DONOR    if strand == '+' else REVERSE_DONOR
    acceptor_motif = CANONICAL_ACCEPTOR if strand == '+' else REVERSE_ACCEPTOR

    # work‑copy of coordinates so we can patch both ends of the intron
    coords = [[start, end] for (_, start, end, _) in exons]
    modified = []

    for i in range(len(exons) - 1):
        chrom, su, eu, _ = exons[i]     # upstream exon
        chrom, sd, ed, _ = exons[i + 1] # downstream exon

        # intron bounds
        intron_start_abs = eu + 1
        intron_end_abs   = sd - 1

        # original donor / acceptor positions (genomic)
        if strand == '+':
            orig_donor_abs    = intron_start_abs
            orig_acceptor_abs = intron_end_abs - 1  # first base of AG
            donor_orientation    = +1
            acceptor_orientation = -1
        else:
            orig_donor_abs    = intron_end_abs - 1  # first base of AC
            orig_acceptor_abs = intron_start_abs
            donor_orientation    = -1
            acceptor_orientation = +1

        # relative coords to gene_start
        rel_donor    = orig_donor_abs    - gene_start
        rel_acceptor = orig_acceptor_abs - gene_start
        rel_intr_start = intron_start_abs - gene_start
        rel_intr_end   = intron_end_abs   - gene_start

        # search first within ±window
        d_new_rel = adjust_site(seq, rel_donor, donor_motif,
                                window, donor_orientation,
                                rel_intr_start, rel_intr_end)
        a_new_rel = adjust_site(seq, rel_acceptor, acceptor_motif,
                                window, acceptor_orientation,
                                rel_intr_start, rel_intr_end)

        # if no shift or same position, widen scan over entire intron
        if d_new_rel is None or abs((d_new_rel + gene_start) - orig_donor_abs) < min_shift:
            d_new_rel = find_canonical_site(seq, rel_donor, donor_motif,
                                            donor_orientation,
                                            rel_intr_start, rel_intr_end,
                                            min_shift) or rel_donor
        if a_new_rel is None or abs((a_new_rel + gene_start) - orig_acceptor_abs) < min_shift:
            a_new_rel = find_canonical_site(seq, rel_acceptor, acceptor_motif,
                                            acceptor_orientation,
                                            rel_intr_start, rel_intr_end,
                                            min_shift) or rel_acceptor

        d_new_abs = gene_start + d_new_rel
        a_new_abs = gene_start + a_new_rel

        # patch exon boundaries (first/last exon untouched)
        if strand == '+':
            # upstream exon end at one base before new donor motif
            coords[i][1]      = d_new_abs - 1
            # downstream exon start just after the end of the new acceptor motif (2 nt)
            coords[i + 1][0]  = a_new_abs + 2
        else:
            # upstream exon end just before the start of the new acceptor motif (2 nt)
            coords[i][1]      = a_new_abs - 1
            # downstream exon start just after the end of the new donor motif (2 nt)
            coords[i + 1][0]  = d_new_abs + 2

        # log motifs
        modified.append((
            orig_donor_abs, d_new_abs,
            orig_acceptor_abs, a_new_abs,
            motif_at(seq, rel_donor),
            motif_at(seq, d_new_abs - gene_start),
            motif_at(seq, rel_acceptor),
            motif_at(seq, a_new_abs - gene_start)
        ))

    # rebuild exon list with updated coords
    new_exons = []
    for (chrom, _, __, strand), (new_start, new_end) in zip(exons, coords):
        new_exons.append((chrom, new_start, new_end, strand))

    return new_exons, modified


def main():
    args = parse_args()
    # initialise detailed logging
    setup_logging(args.log)
    logging.info("Script arguments received: %s", sys.argv)
    print("Script arguments received:", sys.argv)
    print("Loading target gene list...")
    target_genes = load_target_gene_list(args.target_genes, args.input_format, args.symbol_column)
    
    # Build gffutils DB
    import gffutils
    db_fn = args.gtf + '.db'
    if not os.path.exists(db_fn):
        print("Creating GFF database...")
        gffutils.create_db(args.gtf, db_fn, disable_infer_transcripts=True,
                           disable_infer_genes=True)
    print("Loading GFF database...")
    db = gffutils.FeatureDB(db_fn)
    
    fake_feats = []
    full_feats = []
    logs = []
    
    # Get total number of genes for progress bar
    total_genes = len(list(db.features_of_type('gene')))
    print(f"Processing {total_genes} genes...")
    
    for gene in tqdm(db.features_of_type('gene'), total=total_genes, desc="Processing genes"):
        target_tx_feats = []
        target_exon_feats = []

        # gather possible identifiers and match target gene list case‑insensitively
        gene_name_fields = []
        for key in ('gene_name', 'gene_id', 'Name'):
            if key in gene.attributes:
                gene_name_fields.append(gene.attributes[key][0])
                # Also add version-stripped gene_id if applicable
                if key == 'gene_id' and '.' in gene.attributes[key][0]:
                    gene_id_no_version = gene.attributes[key][0].split('.')[0]
                    gene_name_fields.append(gene_id_no_version)
        
        # Also check gene.id and its version-stripped form
        if gene.id not in gene_name_fields:
            gene_name_fields.append(gene.id)
        if '.' in gene.id:
            gene_id_no_version = gene.id.split('.')[0]
            if gene_id_no_version not in gene_name_fields:
                gene_name_fields.append(gene_id_no_version)
                
        gene_name = gene_name_fields[0] if gene_name_fields else ''
        is_target = any(name.upper() in target_genes for name in gene_name_fields)
        logging.debug("Processing gene %s (%s) – target=%s",
                      gene.id, gene_name, is_target)
        
        # Debug for specific gene ENSG00000165572
        if gene.id.startswith('ENSG00000165572') or 'ENSG00000165572' in str(gene.attributes):
            print("\n\nFOUND TARGET GENE FOR DEBUGGING:", gene.id)
            print(f"Gene attributes: {gene.attributes}")
            print(f"Gene name fields: {gene_name_fields}")
            print(f"Is recognized as target: {is_target}")
            print(f"Target genes list (sample): {list(target_genes)[:5]}")
            print(f"Gene name fields uppercase: {[name.upper() for name in gene_name_fields]}")
            print(f"Gene ID without version: {gene.id.split('.')[0]}")
            
            # Try checking if the gene ID without version is in the target genes
            gene_id_no_version = gene.id.split('.')[0]
            is_target_no_version = gene_id_no_version.upper() in target_genes
            print(f"Is target when ignoring version: {is_target_no_version}")
        
        # Get all transcripts for this gene
        transcripts = list(db.children(gene, featuretype='transcript'))
        modified_transcripts = []
        
        for tx in tqdm(transcripts, desc=f"Processing transcripts for {gene.id}", leave=False):
            # Get all features for this transcript
            tx_features = list(db.children(tx, order_by='start'))
            
            if is_target:
                # For target genes, we'll modify exon coordinates
                modified_exons = []
                seq = fetch_sequence(args.genome, gene.seqid, gene.start-1, gene.end)
                gene_start = gene.start
                
                # Process exons
                exons = [f for f in tx_features if f.featuretype == 'exon']
                try:
                    if len(exons) == 1:
                        new_exons, donor, acceptor = process_single_exon(
                            gene.seqid, exons[0].start, exons[0].end, gene.strand, seq, gene_start)
                        donor_motif    = CANONICAL_DONOR    if gene.strand == '+' else REVERSE_DONOR
                        acceptor_motif = CANONICAL_ACCEPTOR if gene.strand == '+' else REVERSE_ACCEPTOR
                        donor_motif_str = motif_at(seq, donor - gene_start)
                        acceptor_motif_str = motif_at(seq, acceptor - gene_start)
                        donor_lab = label_site(donor_motif_str, donor_motif) + '_DONOR'
                        accept_lab = label_site(acceptor_motif_str, acceptor_motif) + '_ACCEPTOR'
                        log = (f"{gene.id}\t{tx.id}\tsingle-exon\t"
                               f"new_donor={donor}({donor_lab} = '{donor_motif_str}'),"
                               f"new_acceptor={acceptor}({accept_lab} = '{acceptor_motif_str}')")
                    else:
                        exlist = [(ex.chrom, ex.start, ex.end, ex.strand) for ex in exons]
                        new_exons, modified = process_multi_exon(exlist, seq, gene_start)
                        donor_motif    = CANONICAL_DONOR    if gene.strand == '+' else REVERSE_DONOR
                        acceptor_motif = CANONICAL_ACCEPTOR if gene.strand == '+' else REVERSE_ACCEPTOR
                        mods_parts = []
                        for (od, nd, oa, na, od_m, nd_m, oa_m, na_m) in modified:
                            od_label = label_site(od_m, donor_motif) + '_DONOR'
                            nd_label = label_site(nd_m, donor_motif) + '_DONOR'
                            oa_label = label_site(oa_m, acceptor_motif) + '_ACCEPTOR'
                            na_label = label_site(na_m, acceptor_motif) + '_ACCEPTOR'
                            mods_parts.append(
                                f"{od}({od_label}, '{od_m}')->{nd}({nd_label}, '{nd_m}')|"
                                f"{oa}({oa_label}, '{oa_m}')->{na}({na_label}, '{na_m}')"
                            )
                        mods_str = ';'.join(mods_parts)
                        log = f"{gene.id}\t{tx.id}\tmulti-exon\t{mods_str}"
                    logs.append(log)
                    
                    # Create modified exon features
                    for i, (c, s, e, strnd) in enumerate(new_exons, 1):
                        # Keep all original attributes from the original exon
                        attrs = {}
                        for key, value in exons[min(i-1, len(exons)-1)].attributes.items():
                            attrs[key] = value
                        attrs['exon_number'] = [str(i)]
                        # Create a new exon feature
                        feat = gffutils.Feature(
                            seqid=c, source=tx.source, featuretype='exon',
                            start=s, end=e, strand=strnd, frame='.',
                            attributes=attrs)
                        modified_exons.append(feat)
                        target_exon_feats.append(feat)
                    
                    # Update transcript coordinates based on modified exons
                    if modified_exons:
                        tx_start = min(exon.start for exon in modified_exons)
                        tx_end = max(exon.end for exon in modified_exons)
                        tx_attrs = dict(tx.attributes)
                        modified_tx = gffutils.Feature(
                            seqid=tx.seqid,
                            source=tx.source,
                            featuretype='transcript',
                            start=tx_start,
                            end=tx_end,
                            strand=tx.strand,
                            frame='.',
                            attributes=tx_attrs
                        )
                        modified_transcripts.append(modified_tx)
                        target_tx_feats.append(modified_tx)
                except ValueError as e:
                    # Log the error and skip this transcript
                    error_msg = f"{gene.id}\t{tx.id}\tERROR\t{str(e)}"
                    logs.append(error_msg)
                    logging.warning(f"Skipping transcript {tx.id} of gene {gene.id}: {str(e)}")
                    
                    # Use the original features for this transcript
                    for feat in tx_features:
                        full_feats.append(feat)
                
                # Keep original ancillary features only for non‑target genes
                if not is_target:
                    for feat in tx_features:
                        if feat.featuretype not in ('exon', 'transcript'):
                            full_feats.append(feat)
            else:
                # For non-target genes, copy all features as is
                for feat in tx_features:
                    full_feats.append(feat)
        
        # Update gene coordinates based on all transcripts (modified or original)
        if is_target:
            # Make a fresh gene feature after transcripts/exons are defined
            if target_tx_feats or target_exon_feats:
                gene_start = min(feat.start for feat in target_exon_feats + target_tx_feats)
                gene_end   = max(feat.end   for feat in target_exon_feats + target_tx_feats)
                gene_attrs = dict(gene.attributes)
                modified_gene = gffutils.Feature(
                    seqid=gene.seqid,
                    source=gene.source,
                    featuretype='gene',
                    start=gene_start,
                    end=gene_end,
                    strand=gene.strand,
                    frame='.',
                    attributes=gene_attrs
                )
                ordered_feats = [modified_gene] + target_tx_feats + target_exon_feats
                fake_feats.extend(ordered_feats)
                full_feats.extend(ordered_feats)
            else:
                # No exons found – fall back to original record
                full_feats.append(gene)
        else:
            # Non‑target gene: keep original order/features
            full_feats.append(gene)

    print("Writing output files...")
    # Write outputs
    write_gtf(args.output_fake, fake_feats)
    write_gtf(args.output_full, full_feats)
    
    # Write log file
    with open(args.log, 'a') as lh:
        lh.write("#gene\ttranscript\tmodification\n")
        lh.write("\n".join(logs))
    print("Done!")


def fetch_sequence(fasta, chrom, start, end):
    # simple FASTA index-based fetcher (requires pysam or pyfaidx)
    from pyfaidx import Fasta
    fa = Fasta(fasta)
    return str(fa[chrom][start:end].seq)


def _format_attributes(attrs):
    """
    Convert a gffutils attributes dict to a valid GTF attribute string.
    """
    parts = []
    for key, values in attrs.items():
        for v in values:
            parts.append(f'{key} "{v}";')
    return ' '.join(parts)

def write_gtf(path, features):
    """
    Write a coordinate‑sorted GTF file (gene → transcript → exon ordering),
    guaranteeing a valid structure for IGV.
    """
    def _sort_key(f):
        feature_rank = 0 if f.featuretype == 'gene' else 1 if f.featuretype == 'transcript' else 2
        return (f.seqid, f.start, feature_rank, f.end)

    with open(path, 'w') as fh:
        for f in sorted(features, key=_sort_key):
            score = f.score if f.score not in (None, '') else '.'
            frame = f.frame if f.frame not in (None, '') else '.'
            attr_str = _format_attributes(f.attributes)
            fh.write('\t'.join([
                str(f.seqid),
                str(f.source),
                str(f.featuretype),
                str(f.start),
                str(f.end),
                str(score),
                str(f.strand),
                str(frame),
                attr_str
            ]) + '\n')

def write_feature(path, feature):
    # append a single feature
    with open(path, 'a') as fh:
        fh.write(str(feature) + '\n')

if __name__ == '__main__':
    main()
