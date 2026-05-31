import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

def read_all_loop_analysis_files(input_folder):
    pattern = os.path.join(input_folder, "*_loop_size_analysis.txt")
    files = glob.glob(pattern)
    
    if len(files) == 0:
        raise ValueError(f"No *_loop_size_analysis.txt files found in {input_folder}")
    
    all_data = []
    file_info = []
    
    for file_path in files:
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#')
            filename = os.path.basename(file_path)
            df['source_file'] = filename
            sample_name = filename.replace('_loop_size_analysis.txt', '')
            df['sample'] = sample_name
            
            all_data.append(df)
            file_info.append((filename, len(df)))
            
        except Exception as e:
            print(f"  Error reading {file_path}: {e}")
            continue
    
    if len(all_data) == 0:
        raise ValueError("No valid data found in any files")
    
    combined_df = pd.concat(all_data, ignore_index=True)
    
    return combined_df, file_info

def create_combined_loop_distribution_plot(combined_df, output_folder, n_bins=9):
    os.makedirs(output_folder, exist_ok=True)
    
    lifespan_min = combined_df['lifespan'].min()
    lifespan_max = combined_df['lifespan'].max()
    bin_edges = np.linspace(lifespan_min, lifespan_max, n_bins + 1)
    
    bin_labels = [f"({bin_edges[i]:.2f}, {bin_edges[i+1]:.2f})" for i in range(n_bins)]
    
    grouped_data = []
    bin_stats = []
    
    for i in range(n_bins):
        low, high = bin_edges[i], bin_edges[i+1]
        if i == n_bins - 1:
            mask = (combined_df['lifespan'] >= low) & (combined_df['lifespan'] <= high)
        else:
            mask = (combined_df['lifespan'] >= low) & (combined_df['lifespan'] < high)
        
        bin_data = combined_df[mask]['loop_size_bins'].values
        grouped_data.append(bin_data)
        
        if len(bin_data) > 0:
            bin_stats.append({
                'range': f"({low:.2f}, {high:.2f})",
                'count': len(bin_data),
                'median': np.median(bin_data),
                'mean': np.mean(bin_data),
                'q25': np.percentile(bin_data, 25),
                'q75': np.percentile(bin_data, 75)
            })
        else:
            bin_stats.append({
                'range': f"({low:.2f}, {high:.2f})",
                'count': 0,
                'median': 0,
                'mean': 0,
                'q25': 0,
                'q75': 0
            })
    
    # Changed figure size to (10, 8)
    fig, ax = plt.subplots(figsize=(14, 12))
    
    box_plot = ax.boxplot(grouped_data, 
                         labels=bin_labels,
                         patch_artist=True,
                         showfliers=True,
                         flierprops=dict(marker='D', markerfacecolor='gray', 
                                       markeredgecolor='gray', markersize=3, alpha=0.5))
    
    box_color = '#7FB3A3'
    for patch in box_plot['boxes']:
        patch.set_facecolor(box_color)
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')
        patch.set_linewidth(1)
    
    for median in box_plot['medians']:
        median.set_color('black')
        median.set_linewidth(2)
    
    for whisker in box_plot['whiskers']:
        whisker.set_color('black')
        whisker.set_linewidth(1)
    
    for cap in box_plot['caps']:
        cap.set_color('black')
        cap.set_linewidth(1)
    
    ax.set_yscale('log')
    
    # Updated label font sizes to 20
    ax.set_xlabel('Lifespan range', fontsize=26, fontweight='bold')
    ax.set_ylabel('Genomic size of loop (log scale)', fontsize=26, fontweight='bold')
    
    # Removed grid
    # ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # Updated tick labels: rotation=90 for vertical, fontsize=18
    plt.xticks(rotation=90, ha='center', fontsize=24)
    ax.tick_params(axis='y', labelsize=24)
    
    # Remove minor ticks
    ax.minorticks_off()
    
    # Yeni y_min ve y_max hesaplama
    true_min = combined_df['loop_size_bins'].min()
    y_min = max(1, true_min * 0.8)
    y_max = combined_df['loop_size_bins'].max() * 2

    if y_max <= 100:
        ax.set_ylim(y_min, 100)
        ax.set_yticks([10, 100])  
        ax.set_yticklabels(['10', '10²'])
    elif y_max <= 1000:
        ax.set_ylim(y_min, 1000)
        ax.set_yticks([10, 100, 1000])
        ax.set_yticklabels(['10', '10²', '10³'])
    else:
        ax.set_ylim(y_min, y_max)
        log_ticks = [10**i for i in range(int(np.log10(y_min)), int(np.log10(y_max))+1)]
        ax.set_yticks(log_ticks)
        ax.set_yticklabels([f'10^{int(np.log10(tick))}' if tick > 1 else '1' for tick in log_ticks])

    plt.tight_layout()
    
    output_file = os.path.join(output_folder, "figure2a_genomic_loops_distribution.png")
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    print(f"Plot saved as: {output_file}")
    
    output_pdf = os.path.join(output_folder, "figure2a_genomic_loops_distribution.pdf")
    plt.savefig(output_pdf, dpi=600, bbox_inches='tight')
    print(f"Plot saved as: {output_pdf}")
    
    plt.show()
    
    stats_file = os.path.join(output_folder, "combined_loop_statistics.txt")
    with open(stats_file, 'w') as f:
        f.write("Combined Genomic Loop Analysis Statistics\n")
        f.write("="*50 + "\n\n")
        f.write(f"Total loops analyzed: {len(combined_df)}\n")
        f.write(f"Number of files combined: {len(combined_df['sample'].unique())}\n")
        f.write(f"Lifespan range: {combined_df['lifespan'].min():.4f} - {combined_df['lifespan'].max():.4f}\n")
        f.write(f"Loop size range: {combined_df['loop_size_bins'].min()} - {combined_df['loop_size_bins'].max()}\n\n")
        
        f.write("Lifespan Bin Statistics:\n")
        f.write("-" * 30 + "\n")
        for stat in bin_stats:
            f.write(f"Range {stat['range']}: {stat['count']} loops, "
                   f"median={stat['median']:.1f}, mean={stat['mean']:.1f}\n")
        
        f.write("\nSample Statistics:\n")
        f.write("-" * 20 + "\n")
        for sample in combined_df['sample'].unique():
            sample_data = combined_df[combined_df['sample'] == sample]
            f.write(f"{sample}: {len(sample_data)} loops\n")
    
    print(f"Statistics saved as: {stats_file}")
    
    return fig, ax, bin_stats

def main():
    input_folder = "your_input_location" # write your own
    output_folder = "your_output_location" #write your own
    
    try:
        combined_df, file_info = read_all_loop_analysis_files(input_folder)
        
        fig, ax, stats = create_combined_loop_distribution_plot(combined_df, output_folder, n_bins=9)
        
        print("\nAnalysis complete!")
        print(f"Results saved in: {output_folder}")
        
    except Exception as e:
        print(f"Error: {e}")
        return None
    
    return combined_df

if __name__ == "__main__":
    combined_data = main()