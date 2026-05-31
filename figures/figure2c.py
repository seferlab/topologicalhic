import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import re

# Set up the paths
input_directory = "/figure2_outputs" # you should update here
output_directory = "your_output_file" # you should update here

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Function to extract chromosome number from filename
def extract_chromosome_number(filename):
    """Extract chromosome number from filename"""
    match = re.search(r'chr(\d+)', filename)
    if match:
        return int(match.group(1))
    return None

# Read all loop size analysis files
def read_all_chromosome_data():
    """Read all chromosome data files and combine them"""
    all_data = []
    
    # Get all files in the input directory
    file_pattern = os.path.join(input_directory, "*loop_size_analysis.txt")
    files = glob.glob(file_pattern)
    
    for file_path in files:
        filename = os.path.basename(file_path)
        chr_num = extract_chromosome_number(filename)
        
        if chr_num is not None and 1 <= chr_num <= 22:
            try:
                # Read the file
                df = pd.read_csv(file_path, sep='\t', comment='#')
                
                # Filter for lifespan >= 0.1 as mentioned in the figure caption
                df_filtered = df[df['lifespan'] >= 0.1].copy()
                
                # Add chromosome information
                df_filtered['chromosome'] = chr_num
                
                all_data.append(df_filtered)
                
            except Exception as e:
                pass  # Silent error handling
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        return combined_df
    else:
        return None

# Create the boxplot
def create_chromosome_boxplot(data):
    """Create customized boxplot with square border and adjusted origin"""
    
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(14,12))
    
    chromosomes = sorted(data['chromosome'].unique())
    plot_data = [data[data['chromosome'] == chr_num]['loop_size_bins'].values for chr_num in chromosomes]
    
    bp = ax.boxplot(plot_data, 
                    positions=chromosomes,
                    widths=0.6,
                    patch_artist=True,
                    showfliers=True,
                    flierprops=dict(marker='d', markersize=3, alpha=0.6, color='gray'),
                    medianprops=dict(color='black', linewidth=1.5),
                    boxprops=dict(facecolor='#7FB3A3', alpha=0.7, edgecolor='black', linewidth=1.2),
                    whiskerprops=dict(color='black', linewidth=1.2),
                    capprops=dict(color='black', linewidth=1.2))

    # Log scale and manual origin adjustment
    ax.set_yscale('log')
    
    min_val = data['loop_size_bins'].min()
    y_min = max(1, min_val * 0.8) 

    y_max = 1000 #you can change here and make it more generic
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0.5, 22.5)

    # Labels and ticks
    ax.set_xlabel('Chromosome number', fontsize=26, fontweight='bold')
    ax.set_xticks(range(1, 23))
    ax.set_xticklabels(range(1, 23), fontsize=22)
    ax.tick_params(axis='y', labelsize=22)

    # Remove grid
    ax.grid(False)

    # Make all borders visible and equal width
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
        spine.set_edgecolor('black')

    # Make the plot area square (equal aspect)
    ax.set_box_aspect(1)

    # Save
    plt.tight_layout()
    output_path = os.path.join(output_directory, 'genomic_loop_sizes_by_chromosome.png')
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    output_path_pdf = os.path.join(output_directory, 'genomic_loop_sizes_by_chromosome.pdf')
    plt.savefig(output_path_pdf, bbox_inches='tight')

    plt.close()  # Close figure to avoid showing

    print(f"Plot saved to: {output_path}")
    print(f"PDF saved to: {output_path_pdf}")
    
    return fig

# Generate summary statistics
def generate_summary_stats(data):
    """Generate summary statistics for each chromosome"""
    summary_stats = []
    
    for chr_num in sorted(data['chromosome'].unique()):
        chr_data = data[data['chromosome'] == chr_num]['loop_size_bins']
        
        stats = {
            'chromosome': chr_num,
            'count': len(chr_data),
            'mean': chr_data.mean(),
            'median': chr_data.median(),
            'std': chr_data.std(),
            'min': chr_data.min(),
            'max': chr_data.max(),
            'q25': chr_data.quantile(0.25),
            'q75': chr_data.quantile(0.75)
        }
        summary_stats.append(stats)
    
    summary_df = pd.DataFrame(summary_stats)
    
    # Save summary statistics
    summary_path = os.path.join(output_directory, 'chromosome_loop_size_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    
    print(f"Summary statistics saved to: {summary_path}")
    
    return summary_df

# Main execution
def main():
    """Main function to execute the analysis"""
    # Read all data
    combined_data = read_all_chromosome_data()
    
    if combined_data is not None:
        # Generate summary statistics
        summary_stats = generate_summary_stats(combined_data)
        
        # Create the plot
        fig = create_chromosome_boxplot(combined_data)

# Run the analysis
if __name__ == "__main__":
    main()