#!/usr/bin/env python3
"""
Extract TUSCO filtering thresholds from log files.

This script parses all TUSCO log files in the expression directory and extracts
key filtering thresholds used in each configuration.
"""

import re
import os
from pathlib import Path
import pandas as pd
from typing import Dict, Optional, Tuple


def extract_bgee_thresholds(log_content: str) -> Tuple[Optional[float], Optional[int], Optional[int]]:
    """
    Extract Bgee prevalence and tissue counts.
    
    Returns: (prevalence, tissues_required, tissues_total)
    """
    # Pattern: "Bgee universal genes (in >= 97/130 universal tissues, prevalence >= 0.75): 31607"
    pattern = r"Bgee universal genes \(in >= (\d+)/(\d+) universal tissues, prevalence >= ([\d.]+)\)"
    match = re.search(pattern, log_content)
    if match:
        tissues_required = int(match.group(1))
        tissues_total = int(match.group(2))
        prevalence = float(match.group(3))
        return prevalence, tissues_required, tissues_total
    return None, None, None


def extract_gtex_prevalence(log_content: str) -> Optional[float]:
    """
    Extract GTEx prevalence threshold (hsa only).
    
    Looks for the first occurrence of "genes pass prevalence ≥ X (TPM > 0.1)"
    """
    # Pattern: "UBERON:0001159: 16422/59033 genes pass prevalence ≥ 0.75 (TPM > 0.1) and median ≥ 1.0 TPM (n=419)"
    pattern = r"genes pass prevalence ≥ ([\d.]+) \(TPM > 0\.1\)"
    match = re.search(pattern, log_content)
    if match:
        return float(match.group(1))
    return None


def extract_final_universal_genes(log_content: str) -> Optional[int]:
    """
    Extract final universal genes count (intersection).
    """
    # Pattern: "Final universal genes (intersection): 66"
    pattern = r"Final universal genes \(intersection\): (\d+)"
    match = re.search(pattern, log_content)
    if match:
        return int(match.group(1))
    return None


def extract_final_transcripts(log_content: str) -> Optional[int]:
    """
    Extract final transcript count.
    """
    # Pattern: "TUSCO selection completed successfully – 64 final transcripts"
    pattern = r"TUSCO selection completed.*? – (\d+) final transcripts"
    match = re.search(pattern, log_content)
    if match:
        return int(match.group(1))
    return None


def extract_genes_removed_step5(log_content: str) -> Optional[int]:
    """
    Extract number of genes removed in Step 5.
    """
    # Pattern: "Total genes removed in Step 5: 717"
    pattern = r"Total genes removed in Step 5: (\d+)"
    match = re.search(pattern, log_content)
    if match:
        return int(match.group(1))
    return None


def parse_log_file(log_path: Path) -> Dict[str, any]:
    """
    Parse a single log file and extract all thresholds.
    
    Returns a dictionary with extracted values.
    """
    with open(log_path, 'r') as f:
        content = f.read()
    
    # Determine config name and species from path
    # e.g., "hsa_prev_lenient/tusco_hsa.log" -> config="hsa_prev_lenient", species="hsa"
    config_name = log_path.parent.name
    if config_name.startswith('hsa'):
        species = 'hsa'
    elif config_name.startswith('mmu'):
        species = 'mmu'
    else:
        species = 'unknown'
    
    # Extract thresholds
    bgee_prevalence, bgee_tissues_required, bgee_tissues_total = extract_bgee_thresholds(content)
    gtex_prevalence = extract_gtex_prevalence(content) if species == 'hsa' else None
    final_universal_genes = extract_final_universal_genes(content)
    final_transcripts = extract_final_transcripts(content)
    genes_removed_step5 = extract_genes_removed_step5(content)
    
    return {
        'config': config_name,
        'species': species,
        'bgee_prevalence': bgee_prevalence,
        'bgee_tissues_required': bgee_tissues_required,
        'bgee_tissues_total': bgee_tissues_total,
        'gtex_prevalence': gtex_prevalence,
        'final_universal_genes': final_universal_genes,
        'final_transcripts': final_transcripts,
        'genes_removed_step5': genes_removed_step5,
    }


def main():
    """Main function to process all log files and generate summary."""
    # Get the directory containing this script
    script_dir = Path(__file__).parent
    
    # Find all log files matching pattern: */tusco_*.log
    log_files = []
    for log_file in script_dir.glob('*/tusco_*.log'):
        log_files.append(log_file)
    
    # Sort for consistent output
    log_files.sort()
    
    if not log_files:
        print(f"Error: No log files found in {script_dir}")
        return
    
    print(f"Found {len(log_files)} log files to process:")
    for log_file in log_files:
        print(f"  - {log_file}")
    
    # Parse all log files
    results = []
    for log_file in log_files:
        print(f"\nProcessing: {log_file.name}")
        try:
            data = parse_log_file(log_file)
            results.append(data)
            print(f"  ✓ Extracted thresholds for {data['config']}")
        except Exception as e:
            print(f"  ✗ Error processing {log_file}: {e}")
            continue
    
    if not results:
        print("\nError: No data extracted from log files")
        return
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Reorder columns
    column_order = [
        'config',
        'species',
        'bgee_prevalence',
        'bgee_tissues_required',
        'bgee_tissues_total',
        'gtex_prevalence',
        'final_universal_genes',
        'final_transcripts',
        'genes_removed_step5',
    ]
    df = df[column_order]
    
    # Sort by species, then config name
    df = df.sort_values(['species', 'config']).reset_index(drop=True)
    
    # Save to TSV
    output_file = script_dir / 'thresholds_summary.tsv'
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\n✓ Successfully extracted thresholds from {len(results)} configurations")
    print(f"✓ Output saved to: {output_file}")
    print(f"\nSummary:")
    print(df.to_string(index=False))


if __name__ == '__main__':
    main()

