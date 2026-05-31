import csv
import glob
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Only define RUES_MES path since that's all we need
RUES_MES_PATH = "/VR/RUES_MES"
CHRNUM = 14  # Only chromosome 14 is used

def readPersisFile(filename):
    persisdiag = []
    with open(filename,'rt') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            if int(line[0]) != 1: continue
            persisdiag.append((float(line[1]),float(line[2])))
    return persisdiag

def find_rues_mes_files():
    """Find RUES_MES files for chromosome 14 only"""
    files = {}
    
    # Find original persistence diagram
    orig_pattern = os.path.join(RUES_MES_PATH, f"*{CHRNUM}_persisdiagram.txt")
    orig_files = glob.glob(orig_pattern)
    # Filter out permutation files
    orig_files = [f for f in orig_files if 'distperm' not in f and 'randperm' not in f and 'edgeperm' not in f]
    if orig_files:
        files['original'] = orig_files[0]
    
    # Find distance permutation
    dist_pattern = os.path.join(RUES_MES_PATH, f"*{CHRNUM}_distperm*_persisdiagram.txt")
    dist_files = glob.glob(dist_pattern)
    if dist_files:
        files['distperm'] = dist_files[0]
    
    # Find random permutation
    rand_pattern = os.path.join(RUES_MES_PATH, f"*{CHRNUM}_randperm*_persisdiagram.txt")
    rand_files = glob.glob(rand_pattern)
    if rand_files:
        files['randperm'] = rand_files[0]
    
    # Find edge permutation
    edge_pattern = os.path.join(RUES_MES_PATH, f"*{CHRNUM}_edgeperm*_persisdiagram.txt")
    edge_files = glob.glob(edge_pattern)
    if edge_files:
        files['edgeperm'] = edge_files[0]
    
    return files

def load_rues_mes_data():
    """Load RUES_MES chromosome 14 data only"""
    files = find_rues_mes_files()
    data = {}
    
    for data_type, filepath in files.items():
        print(f"Loading {data_type}: {os.path.basename(filepath)}")
        data[data_type] = readPersisFile(filepath)
    
    return data

def plot_barcode_figure3a(data, filename):
    """Create barcode plots for Figure 3a - chromosome 14 RUES2 MES cells"""
    
    # Create figure with 2x2 subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Colors matching the original figure
    colors = {
        'original': '#003f5c',    # Dark blue
        'distperm': '#ef5675',    # Pink
        'edgeperm': '#7a5195',    # Purple
        'randperm': '#ffa600'     # Orange
    }
    
    # Plot original data (RUES2 MES, chr 14) - REVERSED ARRAY
    if 'original' in data:
        reversed_orig_data = data['original'][::-1]
        for i, (birth, death) in enumerate(reversed_orig_data):
            ax1.plot([birth, death], [i, i], color=colors['original'], linewidth=1, alpha=0.8)
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, len(reversed_orig_data))
        ax1.set_title('RUES2 MES, chr 14', fontsize=22, fontweight='bold')
        ax1.set_xlabel('Radius', fontsize=20)
        ax1.set_ylabel('Index', fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=20)
    
    # Plot linear dependence - REVERSED ARRAY
    if 'distperm' in data:
        reversed_dist_data = data['distperm'][::-1]
        for i, (birth, death) in enumerate(reversed_dist_data):
            ax2.plot([birth, death], [i, i], color=colors['distperm'], linewidth=1, alpha=0.8)
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, len(reversed_dist_data))
        ax2.set_title('Linear dependence', fontsize=22, fontweight='bold')
        ax2.set_xlabel('Radius', fontsize=20)
        ax2.set_ylabel('Index', fontsize=20)
        ax2.tick_params(axis='both', which='major', labelsize=20)
    
    # Plot edge permutation - REVERSED ARRAY
    if 'edgeperm' in data:
        reversed_edge_data = data['edgeperm'][::-1]
        for i, (birth, death) in enumerate(reversed_edge_data):
            ax3.plot([birth, death], [i, i], color=colors['edgeperm'], linewidth=1, alpha=0.8)
        ax3.set_xlim(0, 1)
        ax3.set_ylim(0, len(reversed_edge_data))
        ax3.set_title('Edge permutation', fontsize=22, fontweight='bold')
        ax3.set_xlabel('Radius', fontsize=20)
        ax3.set_ylabel('Index', fontsize=20)
        ax3.tick_params(axis='both', which='major', labelsize=20)
    
    # Plot random permutation - REVERSED ARRAY
    if 'randperm' in data:
        reversed_rand_data = data['randperm'][::-1]
        for i, (birth, death) in enumerate(reversed_rand_data):
            ax4.plot([birth, death], [i, i], color=colors['randperm'], linewidth=1, alpha=0.8)
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, len(reversed_rand_data))
        ax4.set_title('Random permutation', fontsize=22, fontweight='bold')
        ax4.set_xlabel('Radius', fontsize=20)
        ax4.set_ylabel('Index', fontsize=20)
        ax4.tick_params(axis='both', which='major', labelsize=20)
    
    # Adjust layout
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight', dpi=600, format='png')
    print(f'Saved Figure 3a barcode plots to {filename}')
    plt.close()

# Main execution
if __name__ == "__main__":
    outputloc = "your_output_path" #write your own
    
    print(f"Loading RUES_MES chromosome {CHRNUM} data...")
    data = load_rues_mes_data()
    
    print(f"Found {len(data)} data types:")
    for data_type, persis_data in data.items():
        print(f"  {data_type}: {len(persis_data)} persistence pairs")
    
    print("Creating Figure 3a barcode plots...")
    plot_barcode_figure3a(data, outputloc + 'figure3a_barcode.png')
    
    print("Analysis complete!")