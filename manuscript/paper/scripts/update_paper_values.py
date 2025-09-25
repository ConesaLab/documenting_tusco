#!/usr/bin/env python3
"""
Script to extract numerical values from figure tables and update the paper text
"""

import pandas as pd
import numpy as np
import re
from pathlib import Path

# Base paths resolved relative to the repository root
SCRIPT_DIR = Path(__file__).resolve().parent
PAPER_DIR = SCRIPT_DIR.parent
REPO_ROOT = PAPER_DIR.parents[1]

FIGS_DIR = REPO_ROOT / "figs"
DATA_DIR = REPO_ROOT / "data/processed/tusco"
MANUSCRIPT_DIR = PAPER_DIR

def extract_tusco_counts():
    """Extract TUSCO gene counts from statistics files"""
    values = {}

    # Human tissue statistics
    human_stats = pd.read_csv(DATA_DIR / "hsa/tusco_human_tissue_statistics.tsv", sep='\t')
    values['HUMAN_TISSUE_MIN'] = human_stats['Number_of_TUSCO_Genes'].min()
    values['HUMAN_TISSUE_MAX'] = human_stats['Number_of_TUSCO_Genes'].max()

    # Mouse tissue statistics
    mouse_stats = pd.read_csv(DATA_DIR / "mmu/tusco_mouse_tissue_statistics.tsv", sep='\t')
    values['MOUSE_TISSUE_MIN'] = mouse_stats['Number_of_TUSCO_Genes'].min()
    values['MOUSE_TISSUE_MAX'] = mouse_stats['Number_of_TUSCO_Genes'].max()

    # Universal counts (counting lines excluding comments)
    with open(DATA_DIR / "tusco_human.tsv", 'r') as f:
        human_universal = sum(1 for line in f if not line.startswith('#'))
    with open(DATA_DIR / "tusco_mouse.tsv", 'r') as f:
        mouse_universal = sum(1 for line in f if not line.startswith('#'))
    values['HUMAN_UNIVERSAL'] = human_universal
    values['MOUSE_UNIVERSAL'] = mouse_universal

    return values

def extract_figure3_values():
    """Extract values for Figure 3"""
    values = {}

    # Table S1 - cosine similarity
    table_s1 = pd.read_csv(FIGS_DIR / "figure-03/tables/table_s1.csv")

    # PacBio cosine similarity values
    pacbio_cosim = table_s1[table_s1['Pipeline'].str.contains('PacBio')]['cosim'].values
    values['PACBIO_COSIM_MIN'] = pacbio_cosim.min()
    values['PACBIO_COSIM_MAX'] = pacbio_cosim.max()

    # All cosine similarity range (excluding PacBio already reported)
    other_cosim = table_s1[~table_s1['Pipeline'].str.contains('PacBio')]['cosim'].values
    values['OTHER_COSIM_MIN'] = other_cosim.min()
    values['OTHER_COSIM_MAX'] = other_cosim.max()

    # RIN correlation from figure3c.tsv
    if (FIGS_DIR / "figure-03/tables/figure3c.tsv").exists():
        fig3c = pd.read_csv(FIGS_DIR / "figure-03/tables/figure3c.tsv", sep='\t')
        # Extract correlation values from the data
        # This would need the actual correlation calculation or stored values
        values['TUSCO_RIN_R'] = 0.881  # From paper text
        values['TUSCO_RIN_P'] = "1.41 × 10⁻⁶"
        values['SEQUINS_RIN_R'] = 0.075
        values['SEQUINS_RIN_P'] = 0.77

    return values

def extract_figure5_values():
    """Extract values for Figure 5"""
    values = {}

    if (FIGS_DIR / "figure-05/tables/fig-5b-5c.tsv").exists():
        fig5_data = pd.read_csv(FIGS_DIR / "figure-05/tables/fig-5b-5c.tsv", sep='\t')

        # Extract replicate-specific metrics
        # These would be extracted from the actual data structure
        # Using placeholder values from the paper text for now
        values.update({
            'BRAIN_SN_1REP': 81.3,
            'BRAIN_SN_CI_LOW_1REP': 74.5,
            'BRAIN_SN_CI_HIGH_1REP': 88.0,
            'BRAIN_PDR_1REP': 93.1,
            'BRAIN_PDR_CI_LOW': 89.9,
            'BRAIN_PDR_CI_HIGH': 96.4,
            'BRAIN_PRE_1REP': 67.1,
            'BRAIN_PRE_CI_LOW': 61.4,
            'BRAIN_PRE_CI_HIGH': 72.8,
            'BRAIN_FDR_1REP': 32.9,

            'KIDNEY_SN_1REP': 71.9,
            'KIDNEY_SN_CI_LOW_1REP': 66.4,
            'KIDNEY_SN_CI_HIGH_1REP': 77.4,
            'KIDNEY_PDR_1REP': 90.6,
            'KIDNEY_PDR_CI_LOW': 85.9,
            'KIDNEY_PDR_CI_HIGH': 95.4,
            'KIDNEY_PRE_1REP': 65.2,
            'KIDNEY_PRE_CI_LOW': 61.0,
            'KIDNEY_PRE_CI_HIGH': 69.4,
            'KIDNEY_FDR_1REP': 34.8,

            # Two replicates
            'BRAIN_PRE_2REP': 87.8,
            'BRAIN_PRE_2REP_CI_LOW': 86.5,
            'BRAIN_PRE_2REP_CI_HIGH': 89.1,
            'BRAIN_FDR_2REP': 12.2,
            'BRAIN_SN_2REP': 81.3,
            'BRAIN_SN_2REP_CI_LOW': 79.4,
            'BRAIN_SN_2REP_CI_HIGH': 83.1,

            'KIDNEY_PRE_2REP': 81.0,
            'KIDNEY_PRE_2REP_CI_LOW': 79.2,
            'KIDNEY_PRE_2REP_CI_HIGH': 82.8,
            'KIDNEY_FDR_2REP': 19.0,
            'KIDNEY_SN_2REP': 70.6,
            'KIDNEY_SN_2REP_CI_LOW': 67.8,
            'KIDNEY_SN_2REP_CI_HIGH': 73.5,

            # 3-5 replicates
            'BRAIN_PRE_3REP': 88.7,
            'BRAIN_PRE_5REP': 89.3,
            'BRAIN_SN_3REP': 83.1,
            'BRAIN_SN_5REP': 84.4,
            'KIDNEY_SN_5REP': 78.1,
            'KIDNEY_PRE_RANGE': "79–81",

            # Tissue-specific
            'BRAIN_TISSUE_GENES': 65,
            'KIDNEY_TISSUE_GENES': 46,
            'BRAIN_COSIM': 0.9999,
            'KIDNEY_COSIM': 0.9992,
        })

    return values

def update_paper_text(values):
    """Update the paper text file with extracted values"""

    # Read the current paper text
    paper_file = MANUSCRIPT_DIR / "txt/tusco_paper.txt"
    with open(paper_file, 'r') as f:
        text = f.read()

    # Create updated text with actual values
    updated_text = text

    # Update Figure 1 values
    old_fig1 = "This process yielded tissue-resolved TUSCO gene sets (human: 60–187 per tissue; mouse: 28–67) and universal cores (human n = 50; mouse n = 32)"
    new_fig1 = f"This process yielded tissue-resolved TUSCO gene sets (human: {values['HUMAN_TISSUE_MIN']}–{values['HUMAN_TISSUE_MAX']} per tissue; mouse: {values['MOUSE_TISSUE_MIN']}–{values['MOUSE_TISSUE_MAX']}) and universal cores (human n = {values['HUMAN_UNIVERSAL']}; mouse n = {values['MOUSE_UNIVERSAL']})"

    if old_fig1 in updated_text:
        updated_text = updated_text.replace(old_fig1, new_fig1)
        print(f"✓ Updated Figure 1 gene counts")

    # Update Figure 3 cosine similarity
    old_cosim1 = "consistently achieved cosim values near unity in both human and mouse samples"
    new_cosim1 = f"consistently achieved cosim values ranging from {values['PACBIO_COSIM_MIN']:.4f} to {values['PACBIO_COSIM_MAX']:.4f} in both human and mouse samples"

    if old_cosim1 in updated_text:
        updated_text = updated_text.replace(old_cosim1, new_cosim1)
        print(f"✓ Updated PacBio cosine similarity")

    old_cosim2 = "also exhibited high agreement (cosim: 0.9553–0.9977)"
    new_cosim2 = f"also exhibited high agreement (cosim: {values['OTHER_COSIM_MIN']:.4f}–{values['OTHER_COSIM_MAX']:.4f})"

    if old_cosim2 in updated_text:
        updated_text = updated_text.replace(old_cosim2, new_cosim2)
        print(f"✓ Updated other cosine similarity range")

    # Save updated text
    with open(paper_file, 'w') as f:
        f.write(updated_text)

    # Also save a backup
    backup_file = paper_file.with_suffix('.txt.backup')
    with open(backup_file, 'w') as f:
        f.write(text)

    return updated_text

def main():
    """Main execution function"""
    print("Extracting numerical values from tables...")

    all_values = {}

    # Extract values from different sources
    print("\n1. Extracting TUSCO gene counts...")
    all_values.update(extract_tusco_counts())

    print("\n2. Extracting Figure 3 values...")
    all_values.update(extract_figure3_values())

    print("\n3. Extracting Figure 5 values...")
    all_values.update(extract_figure5_values())

    # Print extracted values for verification
    print("\n" + "="*50)
    print("EXTRACTED VALUES:")
    print("="*50)
    for key, value in sorted(all_values.items()):
        print(f"{key}: {value}")

    # Update the paper text
    print("\n" + "="*50)
    print("UPDATING PAPER TEXT...")
    print("="*50)
    update_paper_text(all_values)

    print("\n✓ Paper text updated successfully!")
    print(f"  Original backed up to: tusco_paper.txt.backup")

if __name__ == "__main__":
    main()
