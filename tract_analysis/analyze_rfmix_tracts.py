#!/usr/bin/env python3
"""
Script to analyze RFMix output and calculate ancestry tract lengths
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip
import os
import argparse
from collections import defaultdict
import seaborn as sns

def read_hmm_file(hmm_file_path):
    """
    Read the HMM output file to get SNP positions and Morgan distances
    """
    print(f"Reading HMM file: {hmm_file_path}")
    
    # Read the HMM file without pre-defining column names
    hmm_data = pd.read_csv(hmm_file_path, sep='\s+', header=None)
    
    print(f"HMM file has {hmm_data.shape[1]} columns")
    print(f"Loaded {len(hmm_data)} SNPs from HMM file")
    
    # Add meaningful column names for the columns we need
    # Column 0: chromosome, Column 1: position, Column 6: Morgan distance (0-indexed)
    hmm_data.rename(columns={0: 'chr', 1: 'pos', 6: 'morgan_dist'}, inplace=True)
    
    return hmm_data

def read_rfmix_file(rfmix_file_path):
    """
    Read RFMix output file
    """
    print(f"Reading RFMix file: {rfmix_file_path}")
    
    if rfmix_file_path.endswith('.gz'):
        with gzip.open(rfmix_file_path, 'rt') as f:
            lines = f.readlines()
    else:
        with open(rfmix_file_path, 'r') as f:
            lines = f.readlines()
    
    # Parse the data - each line is a SNP, each column is an individual
    rfmix_data = []
    for line in lines:
        if line.strip():
            values = [int(x) for x in line.strip().split()]
            rfmix_data.append(values)
    
    rfmix_array = np.array(rfmix_data)
    print(f"Loaded RFMix data: {rfmix_array.shape[0]} SNPs x {rfmix_array.shape[1]} individuals")
    
    return rfmix_array

def calculate_tract_lengths(rfmix_data, morgan_distances):
    """
    Calculate tract lengths for each individual and ancestry
    """
    n_snps, n_individuals = rfmix_data.shape
    
    if len(morgan_distances) != n_snps:
        print(f"Warning: Mismatch between SNPs in RFMix ({n_snps}) and Morgan distances ({len(morgan_distances)})")
        min_len = min(n_snps, len(morgan_distances))
        rfmix_data = rfmix_data[:min_len, :]
        morgan_distances = morgan_distances[:min_len]
        n_snps = min_len
    
    all_tract_lengths = {'ancestry_1': [], 'ancestry_2': []}
    
    # Calculate cumulative Morgan distances
    cumulative_morgan = np.cumsum([0] + list(morgan_distances[1:]))
    
    for individual in range(n_individuals):
        ancestry_sequence = rfmix_data[:, individual]
        
        # Find tract boundaries
        tract_starts = [0]
        current_ancestry = ancestry_sequence[0]
        
        for snp in range(1, n_snps):
            if ancestry_sequence[snp] != current_ancestry:
                tract_starts.append(snp)
                current_ancestry = ancestry_sequence[snp]
        
        # Calculate tract lengths
        for i in range(len(tract_starts)):
            start_snp = tract_starts[i]
            end_snp = tract_starts[i + 1] if i + 1 < len(tract_starts) else n_snps - 1
            
            # Tract length in Morgans
            tract_length_morgan = cumulative_morgan[end_snp] - cumulative_morgan[start_snp]
            
            # Convert to centiMorgans
            tract_length_cm = tract_length_morgan * 100
            
            # Get ancestry for this tract
            tract_ancestry = ancestry_sequence[start_snp]
            
            if tract_ancestry == 1:
                all_tract_lengths['ancestry_1'].append(tract_length_cm)
            elif tract_ancestry == 2:
                all_tract_lengths['ancestry_2'].append(tract_length_cm)
    
    return all_tract_lengths

def plot_tract_length_distribution(all_tract_lengths, output_file=None):
    """
    Create a plot similar to the reference figure
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Define bins (logarithmic spacing)
    bins = np.logspace(0, 3, 50)  # From 1 to 1000 cM
    
    # Plot histograms for each ancestry
    colors = ['blue', 'red']
    ancestries = ['ancestry_1', 'ancestry_2']
    labels = ['Ancestry 1', 'Ancestry 2']
    
    for i, (ancestry, color, label) in enumerate(zip(ancestries, colors, labels)):
        if len(all_tract_lengths[ancestry]) > 0:
            counts, bin_edges = np.histogram(all_tract_lengths[ancestry], bins=bins)
            
            # Calculate bin centers
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            
            # Plot as scatter points
            ax.scatter(bin_centers[counts > 0], counts[counts > 0], 
                      color=color, alpha=0.7, s=30, label=label)
            
            # Fit exponential decay (optional)
            valid_idx = counts > 0
            if np.sum(valid_idx) > 3:
                x_fit = bin_centers[valid_idx]
                y_fit = counts[valid_idx]
                
                # Fit exponential: y = a * exp(-b*x)
                log_y = np.log(y_fit)
                coeffs = np.polyfit(x_fit, log_y, 1)
                
                # Plot fitted line
                x_smooth = np.linspace(min(x_fit), max(x_fit), 100)
                y_smooth = np.exp(coeffs[1]) * np.exp(coeffs[0] * x_smooth)
                ax.plot(x_smooth, y_smooth, color=color, linewidth=2, alpha=0.8)
    
    # Formatting
    ax.set_xscale('linear')
    ax.set_yscale('log')
    ax.set_xlabel('Tract Length (cM)', fontsize=12)
    ax.set_ylabel('Number of tracts', fontsize=12)
    ax.set_title('Ancestry Tract Length Distribution', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    
    plt.show()
    
    return fig, ax

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Analyze RFMix output and calculate ancestry tract lengths')
    parser.add_argument('--hmm-file', required=True, 
                       help='Path to the HMM file with Morgan distances')
    parser.add_argument('--rfmix-pattern', required=True,
                       help='Path pattern for RFMix files with ${chr} placeholder (e.g., /path/to/output_rfmix_${chr}.txt.gz)')
    parser.add_argument('--output-prefix', default='tract_analysis',
                       help='Prefix for output files (default: tract_analysis)')
    parser.add_argument('--chr-start', type=int, default=1,
                       help='Start chromosome number (default: 1)')
    parser.add_argument('--chr-end', type=int, default=22,
                       help='End chromosome number (default: 22)')
    
    args = parser.parse_args()
    
    # Read HMM file
    hmm_data = read_hmm_file(args.hmm_file)
    
    # Group HMM data by chromosome
    hmm_by_chr = {}
    for chr_num in range(args.chr_start, args.chr_end + 1):
        chr_data = hmm_data[hmm_data['chr'] == chr_num]
        if len(chr_data) > 0:
            hmm_by_chr[chr_num] = chr_data
            print(f"Chromosome {chr_num}: {len(chr_data)} SNPs")
    
    # Collect all tract lengths across all chromosomes
    all_tract_lengths_combined = {'ancestry_1': [], 'ancestry_2': []}
    
    # Process each chromosome
    for chr_num in range(args.chr_start, args.chr_end + 1):
        # Replace ${chr} placeholder with chromosome number
        rfmix_file = args.rfmix_pattern.replace('${chr}', str(chr_num))
        
        if not os.path.exists(rfmix_file):
            print(f"Warning: RFMix file for chromosome {chr_num} not found: {rfmix_file}")
            continue
        
        if chr_num not in hmm_by_chr:
            print(f"Warning: No HMM data for chromosome {chr_num}")
            continue
        
        # Read RFMix data
        rfmix_data = read_rfmix_file(rfmix_file)
        
        # Get Morgan distances for this chromosome
        morgan_distances = hmm_by_chr[chr_num]['morgan_dist'].values
        
        # Calculate tract lengths for this chromosome
        tract_lengths = calculate_tract_lengths(rfmix_data, morgan_distances)
        
        # Add to combined results
        all_tract_lengths_combined['ancestry_1'].extend(tract_lengths['ancestry_1'])
        all_tract_lengths_combined['ancestry_2'].extend(tract_lengths['ancestry_2'])
        
        print(f"Chromosome {chr_num}: Found {len(tract_lengths['ancestry_1'])} ancestry 1 tracts, "
              f"{len(tract_lengths['ancestry_2'])} ancestry 2 tracts")
    
    # Print summary statistics
    print("\n=== SUMMARY STATISTICS ===")
    for ancestry in ['ancestry_1', 'ancestry_2']:
        tracts = all_tract_lengths_combined[ancestry]
        if len(tracts) > 0:
            print(f"{ancestry.replace('_', ' ').title()}:")
            print(f"  Total tracts: {len(tracts)}")
            print(f"  Mean length: {np.mean(tracts):.2f} cM")
            print(f"  Median length: {np.median(tracts):.2f} cM")
            print(f"  Min length: {np.min(tracts):.2f} cM")
            print(f"  Max length: {np.max(tracts):.2f} cM")
            print()
    
    # Create the plot
    output_plot = f"{args.output_prefix}_tract_length_distribution.png"
    plot_tract_length_distribution(all_tract_lengths_combined, output_plot)
    
    # Save tract length data to file
    output_data = f"{args.output_prefix}_tract_lengths_data.txt"
    with open(output_data, 'w') as f:
        f.write("ancestry\ttract_length_cm\n")
        for ancestry in ['ancestry_1', 'ancestry_2']:
            for length in all_tract_lengths_combined[ancestry]:
                f.write(f"{ancestry}\t{length:.6f}\n")
    
    print(f"Tract length data saved to: {output_data}")

if __name__ == "__main__":
    main()
