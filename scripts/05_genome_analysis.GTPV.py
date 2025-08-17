#!/usr/bin/env python3
#Genome Region Pairwise Analysis Script
#Processes MAFFT-aligned FASTA files and creates pairwise similarity matrices
#for specified genomic regions
import sys
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import subprocess

def parse_fasta(fasta_file):
    """Parse FASTA file and return sequences dictionary"""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract accession from header (handle #1#NC_004003.1 format)
        seq_id = record.id
        if '#' in seq_id:
            parts = seq_id.split('#')
            seq_id = parts[-1] if parts[-1] else parts[0]

        sequences[seq_id] = str(record.seq).upper()

    return sequences

def get_clade_mapping():
    """Return the clade mapping for samples"""
    return {
        'KC951854.1': '2.1',
        'MH381810.1': '2.1',
        'MN072620.1': '2.1',
        'MN072621.1': '2.1',
        'MW020570.1': '2.1',
        'PV167794.1': '2.1',
        'AY077835.1': '2.2',
        'AY077836.1': '2.2',
        'KX576657.1': '2.2',
        'MN072622.1': '2.2',
        'NC_004003.1': '2.2',
        'MN072623.1': '2.3',
        'MN072624.1': '2.3',
        'MN072625.1': '2.3'
    }

def get_regions():
    """Return regions of interest"""
    return {
#         'LSDV009': (5583, 5869),
#         'LSDV013': (8395, 9465),
#         'LSDV026': (17410, 18341)
        'LSDV132': (119053, 119613),
        'LSDV136': (128526, 128998)
    }

def get_clade_colors():
    """Return color mapping for clades"""
    return {
        '2.1': '#404040',  # Dark grey
        '2.2': '#8E44AD',  # Purple
        '2.3': '#E67E22'   # Orange
    }

def extract_region_from_alignment(aligned_seq, reference_seq, start, end):
    """
    Extract region from aligned sequence based on reference coordinates.
    Find alignment columns corresponding to reference coordinates, then extract same columns from all sequences.
    """
    # Find the alignment columns that correspond to the reference coordinates
    ref_pos = 0  # Position in unaligned reference sequence (1-based)
    align_start = None
    align_end = None

    for i, char in enumerate(reference_seq):
        if char != '-':  # Only count non-gap characters
            ref_pos += 1
            if ref_pos == start and align_start is None:
                align_start = i
            if ref_pos == end:
                align_end = i + 1
                break

    if align_start is None or align_end is None:
        return ""

    # Extract the same alignment columns from the target sequence
    if align_end > len(aligned_seq):
        align_end = len(aligned_seq)

    return aligned_seq[align_start:align_end]

def extract_region_direct(sequence, start, end):
    """
    Extract region directly from sequence using 1-based coordinates
    """
    return sequence[start-1:end]

def calculate_pairwise_similarity(seq1, seq2):
    """
    Calculate pairwise similarity excluding gaps.
    Only counts positions where both sequences have non-gap characters.
    """
    if len(seq1) != len(seq2):
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]

    matches = 0
    valid_positions = 0

    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            valid_positions += 1
            if seq1[i] == seq2[i]:
                matches += 1

    if valid_positions == 0:
        return 0.0

    return (matches / valid_positions) * 100

def calculate_sliding_window_differences(sequences, reference_id, start, end, window_size=100, step_size=10):
    """
    Calculate differences per 100 bp sliding windows across a region
    """
    clade_map = get_clade_mapping()
    reference_seq = sequences[reference_id]

    # Find alignment columns for the reference coordinates
    ref_pos = 0
    align_start = None
    align_end = None

    for i, char in enumerate(reference_seq):
        if char != '-':
            ref_pos += 1
            if ref_pos == start and align_start is None:
                align_start = i
            if ref_pos == end:
                align_end = i + 1
                break

    if align_start is None or align_end is None:
        return {}

    # Extract regions for all sequences using the same alignment columns
    region_seqs = {}
    for seq_id, seq in sequences.items():
        if align_end > len(seq):
            region_seqs[seq_id] = seq[align_start:]
        else:
            region_seqs[seq_id] = seq[align_start:align_end]

    # Group sequences by clade
    clade_seqs = {}
    for seq_id, seq in region_seqs.items():
        if seq_id == reference_id:
            continue
        clade = clade_map.get(seq_id, 'Unknown')
        if clade not in clade_seqs:
            clade_seqs[clade] = []
        clade_seqs[clade].append(seq)

    # Calculate sliding window differences using the reference region
    ref_region = region_seqs[reference_id]
    results = {}
    region_length = len(ref_region.replace('-', ''))

    for clade, seqs in clade_seqs.items():
        positions = []
        differences = []

        for window_start in range(0, region_length - window_size + 1, step_size):
            window_end = window_start + window_size

            # Get reference window directly from aligned region
            ref_window = ref_region[window_start:window_end]

            total_diffs = 0
            valid_comparisons = 0

            for seq in seqs:
                seq_window = seq[window_start:window_end] if window_end <= len(seq) else seq[window_start:]

                # Count differences in this window
                diffs = 0
                valid_positions = 0

                min_len = min(len(ref_window), len(seq_window))
                for j in range(min_len):
                    if ref_window[j] != '-' and seq_window[j] != '-':
                        valid_positions += 1
                        if ref_window[j] != seq_window[j]:
                            diffs += 1

                if valid_positions > 0:
                    total_diffs += (diffs / valid_positions) * 100
                    valid_comparisons += 1

            if valid_comparisons > 0:
                avg_diff = total_diffs / valid_comparisons
                positions.append(start + window_start + window_size // 2)
                differences.append(avg_diff)

        results[clade] = (positions, differences)

    return results

def create_sliding_window_plot(sequences, region_name, start, end, reference_id='NC_004003.1'):
    """Create sliding window plot for a region"""
    # Calculate sliding window differences
    results = calculate_sliding_window_differences(sequences, reference_id, start, end)

    if not results:
        return

    # Create plot
    plt.figure(figsize=(12, 8))
    colors = get_clade_colors()

    for clade, (positions, differences) in results.items():
        if positions and differences:
            plt.plot(positions, differences, color=colors.get(clade, 'black'),
                     linewidth=2, label=f'Clade {clade}')

    plt.xlabel('NC_004003.1 Coordinates (bp)', fontsize=12)
    plt.ylabel('Differences per 100 bp (%)', fontsize=12)
    plt.title(f'{region_name} - Sliding Window Analysis (100 bp window, 10 bp step)', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save plot
    plot_filename = f"{region_name}_sliding_window.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()

def create_pairwise_matrix(sequences, region_name, start, end, reference_id='NC_004003.1'):
    """Create pairwise similarity matrix for a specific region"""
    clade_map = get_clade_mapping()

    if reference_id not in sequences:
        reference_id = list(sequences.keys())[0]

    reference_seq = sequences[reference_id]

    # Find alignment columns for the reference coordinates
    ref_pos = 0
    align_start = None
    align_end = None

    for i, char in enumerate(reference_seq):
        if char != '-':
            ref_pos += 1
            if ref_pos == start and align_start is None:
                align_start = i
            if ref_pos == end:
                align_end = i + 1
                break

    if align_start is None or align_end is None:
        return None, None

    # Extract the same alignment columns from all sequences
    region_seqs = {}
    for seq_id, seq in sequences.items():
        if align_end > len(seq):
            region_seqs[seq_id] = seq[align_start:]
        else:
            region_seqs[seq_id] = seq[align_start:align_end]

    # Create similarity matrix
    seq_ids = list(sequences.keys())
    seq_ids.sort(key=lambda x: (clade_map.get(x, 'Unknown'), x))

    # Create DataFrame with clade information
    index_labels = [f"{clade_map.get(seq_id, 'Unknown')}_{seq_id}" for seq_id in seq_ids]
    matrix = pd.DataFrame(index=index_labels, columns=index_labels, dtype=float)

    for i, seq_id1 in enumerate(seq_ids):
        for j, seq_id2 in enumerate(seq_ids):
            if i != j:
                similarity = calculate_pairwise_similarity(
                    region_seqs[seq_id1],
                    region_seqs[seq_id2]
                )
                matrix.iloc[i, j] = round(similarity, 1)
            else:
                matrix.iloc[i, j] = np.nan

    return matrix, region_seqs

def print_matrix_lower_triangle(matrix, region_name):
    """Print lower triangle of matrix to console"""
    print(f"\n{region_name} - Pairwise Similarity Matrix (%)")
    print("=" * (len(region_name) + 35))

    # Print header
    header = "Sample".ljust(20)
    for col in matrix.columns:
        header += col.split('_')[1][:10].ljust(12)
    print(header)

    # Print rows (lower triangle only)
    for i, (idx, row) in enumerate(matrix.iterrows()):
        row_str = idx.ljust(20)
        for j, val in enumerate(row):
            if i == j:
                row_str += "-".ljust(12)
            elif j > i:
                row_str += "".ljust(12)
            else:
                if pd.isna(val):
                    row_str += "-".ljust(12)
                else:
                    row_str += f"{val:.1f}%".ljust(12)
        print(row_str)

# --- CORRECTED FUNCTIONS START HERE ---

def save_region_aligned_fasta(region_seqs, region_name, output_prefix="region"):
    """
    Save region aligned sequences to a new FASTA alignment file.
    This function correctly preserves the gaps ('-') to maintain the alignment.
    """
    # The output is an alignment, so we use a .aln extension
    filename = f"{output_prefix}_{region_name}.aln"

    with open(filename, 'w') as f:
        for seq_id, aligned_seq in region_seqs.items():
            if aligned_seq:  # Only write if sequence is not empty
                f.write(f">{seq_id}\n")
                f.write(f"{aligned_seq}\n") # Write the sequence WITH gaps

    return filename

def run_raxml_on_alignment(alignment_file, region_name):
    """
    Run RAxML analysis on a pre-aligned file.
    This function skips the redundant MAFFT step.
    """
    # RAxML analysis - it now takes the alignment file directly
    raxml_cmd = [
        "raxml-ng", "--all", "--msa", alignment_file,
        "--model", "GTR+G", "--prefix", region_name,
        "--seed", "123", "--bs-metric", "fbp,tbe", "--redo"
    ]

    try:
        # Using capture_output=True to hide stdout/stderr unless there's an error
        result = subprocess.run(raxml_cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        print(f"Error: raxml-ng not found. Please ensure it's in your system's PATH.", file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error running RAxML for {region_name}:", file=sys.stderr)
        print(f"Command: {' '.join(e.cmd)}", file=sys.stderr)
        print(f"Stderr:\n{e.stderr}", file=sys.stderr)

# --- CORRECTED FUNCTIONS END HERE ---

def save_matrices_csv(matrices, output_prefix="similarity_matrices"):
    """Save all matrices to CSV files"""
    for region_name, matrix in matrices.items():
        filename = f"{output_prefix}_{region_name}.csv"
        matrix.to_csv(filename)
        print(f"Saved matrix to {filename}")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze pairwise similarity of genomic regions from MAFFT alignment"
    )
    parser.add_argument("fasta_file",
                        help="MAFFT-aligned FASTA file")
    parser.add_argument("--output-prefix", "-o",
                        default="similarity_matrices",
                        help="Output file prefix for CSV files (default: similarity_matrices)")
    parser.add_argument("--reference", "-r",
                        default="NC_004003.1",
                        help="Reference sequence ID for coordinates (default: NC_004003.1)")
    parser.add_argument("--csv-only", "-c",
                        action="store_true",
                        help="Only output CSV files, don't print to console")
    parser.add_argument("--skip-phylo", "-s",
                        action="store_true",
                        help="Skip RAxML analysis")
    parser.add_argument("--skip-plots", "-p",
                        action="store_true",
                        help="Skip sliding window plots")
    parser.add_argument("--region-prefix",
                        default="region",
                        help="Prefix for region FASTA files (default: region)")

    args = parser.parse_args()

    # Check if file exists
    if not os.path.exists(args.fasta_file):
        print(f"Error: Input FASTA file not found at '{args.fasta_file}'", file=sys.stderr)
        sys.exit(1)

    # Parse sequences
    try:
        sequences = parse_fasta(args.fasta_file)
    except Exception as e:
        print(f"Error parsing FASTA file: {e}", file=sys.stderr)
        sys.exit(1)

    # Get regions
    regions = get_regions()

    # Calculate matrices for each region
    matrices = {}
    region_sequences = {}

    for region_name, (start, end) in regions.items():
        try:
            matrix, region_seqs = create_pairwise_matrix(sequences, region_name, start, end, args.reference)
            if matrix is not None and not matrix.empty:
                matrices[region_name] = matrix
                region_sequences[region_name] = region_seqs

                if not args.csv_only:
                    print_matrix_lower_triangle(matrix, region_name)
            else:
                print(f"Warning: Could not generate matrix for region {region_name}. It might be outside the reference sequence range.", file=sys.stderr)
        except Exception as e:
            print(f"An error occurred while processing region {region_name}: {e}", file=sys.stderr)
            pass

    # Save to CSV files
    if matrices:
        save_matrices_csv(matrices, args.output_prefix)

    # Create sliding window plots
    if not args.skip_plots:
        print("\nGenerating sliding window plots...")
        for region_name, (start, end) in regions.items():
            try:
                create_sliding_window_plot(sequences, region_name, start, end, args.reference)
            except Exception as e:
                print(f"Could not create plot for {region_name}: {e}", file=sys.stderr)
                pass

    # Save region alignment files and run phylogenetic analysis
    if region_sequences and not args.skip_phylo:
        print("\nRunning phylogenetic analysis...")
        for region_name, region_seqs in region_sequences.items():
            # Save the already-aligned regions to a new .aln file
            alignment_file = save_region_aligned_fasta(region_seqs, region_name, args.region_prefix)
            # Run RAxML directly on this new alignment file
            run_raxml_on_alignment(alignment_file, region_name)

    # Print summary only if matrices were generated and not csv_only mode
    if matrices:
        print(f"\n--- Analysis Summary ---")
        print(f"Processed {len(sequences)} sequences from {args.fasta_file}")
        print(f"Analyzed {len(matrices)} regions: {', '.join(matrices.keys())}")
        print(f"Generated {len(matrices)} similarity matrices (CSV).")

        if not args.skip_plots:
            print(f"Created sliding window plots (.png).")

        if not args.skip_phylo and region_sequences:
            print(f"Created region alignment files (.aln) and ran RAxML analysis.")
            print(f"RAxML results have prefixes: {', '.join(region_sequences.keys())}")

        # Show clade distribution
        clade_map = get_clade_mapping()
        clade_counts = {}
        for seq_id in sequences.keys():
            clade = clade_map.get(seq_id, 'Unknown')
            clade_counts[clade] = clade_counts.get(clade, 0) + 1

        print("\nClade Distribution:")
        for clade, count in sorted(clade_counts.items()):
            print(f"  Clade {clade}: {count} samples")
        print("------------------------")

    else:
        print("\nNo valid matrices were generated. Please check region coordinates and input file.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
