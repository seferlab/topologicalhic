import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
from pathlib import Path

def read_loop_data_files(directory_path):
    """
    Read all loop size analysis files from the specified directory
    """
    # Get all .txt files in the directory
    file_pattern = os.path.join(directory_path, "*_loop_size_analysis.txt")
    files = glob.glob(file_pattern)
    
    print(f"Found {len(files)} files in {directory_path}")
    for f in files:
        print(f"  - {os.path.basename(f)}")
    
    all_data = []
    
    for file_path in files:
        try:
            # Extract cell type from filename
            filename = os.path.basename(file_path)
            if 'RUES' in filename.upper():
                cell_type = 'RUES2'
            elif 'WTC' in filename.upper():
                cell_type = 'WTC11'
            else:
                print(f"Skipping file (unknown cell type): {filename}")
                continue
            
            print(f"Processing: {filename} -> {cell_type}")
            
            # Try multiple methods to read the file
            df = None
            
            try:
                # Read the file line by line to find the actual data
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                
                print(f"  File has {len(lines)} lines")
                
                # Print first few lines for debugging
                print("  First 5 lines:")
                for i, line in enumerate(lines[:5]):
                    print(f"    {i}: {line.strip()}")
                
                # Look for the actual data header (not the comment line) if exists
                data_start = 0
                actual_header = None
                
                for i, line in enumerate(lines):
                    line = line.strip()
                    # Skip comment lines that start with # if exists 
                    if line.startswith('#'):
                        continue
                    # Look for the line that starts with "loop_id" (actual header)
                    if line.startswith('loop_id'):
                        actual_header = line
                        data_start = i
                        print(f"  Found actual header at line {i}: {line}")
                        break
                
                if actual_header is None:
                    print("  Could not find data header starting with 'loop_id'")
                    # Print first few non-comment lines for debugging
                    print("  First few non-comment lines:")
                    for i, line in enumerate(lines[:15]):
                        if not line.strip().startswith('#'):
                            print(f"    Line {i}: '{line.strip()}'")
                    continue
                
                # Read the data starting from the actual header
                df = pd.read_csv(file_path, skiprows=data_start, sep=r'\s+', engine='python')
                
                # Clean column names (remove any trailing commas or extra characters)
                original_columns = list(df.columns)
                df.columns = [col.rstrip(',') for col in df.columns]
                
                print(f"  Original columns: {original_columns}")
                print(f"  Cleaned columns: {list(df.columns)}")
                
            except Exception as e:
                print(f"  Error reading file: {e}")
                import traceback
                traceback.print_exc()
                continue
            
            if df is None or df.empty:
                print(f"  Failed to read data from {filename}")
                continue
            
            print(f"  Raw data shape: {df.shape}")
            
            # Check if we have the required columns
            required_cols = ['lifespan', 'loop_size_bins']
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                print(f"  Missing required columns: {missing_cols}")
                print(f"  Available columns: {list(df.columns)}")
                continue
            
            # Filter out loops with lifespan < 0.1 (artifacts as mentioned)
            original_len = len(df)
            df_filtered = df[df['lifespan'] >= 0.1].copy()
            
            if len(df_filtered) == 0:
                print(f"  No loops with lifespan >= 0.1 found")
                continue
            
            # Add cell type column
            df_filtered['cell_type'] = cell_type
            
            all_data.append(df_filtered)
            
            print(f"  - Loaded {len(df_filtered)} loops (filtered from {original_len}, lifespan >= 0.1)")
            print(f"  - Loop size range: {df_filtered['loop_size_bins'].min()} - {df_filtered['loop_size_bins'].max()}")
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    if not all_data:
        raise ValueError("No valid data files found! Check file format and directory path.")
    
    # Combine all data
    combined_data = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal loops loaded: {len(combined_data)}")
    print(f"Cell types found: {combined_data['cell_type'].value_counts()}")
    
    return combined_data

def create_boxplot(data, output_path):
    """
    Create a customized box plot with adjusted layout and styling
    """
    fig, ax = plt.subplots(figsize=(8, 12))
    
    # Prepare data
    wtc_data = data[data['cell_type'] == 'WTC11']['loop_size_bins']
    rues_data = data[data['cell_type'] == 'RUES2']['loop_size_bins']
    
    box_data = []
    labels = []
    
    if len(wtc_data) > 0:
        box_data.append(wtc_data[wtc_data > 0])
        labels.append('WTC11')
    
    if len(rues_data) > 0:
        box_data.append(rues_data[rues_data > 0])
        labels.append('RUES2')
    
    if not box_data:
        raise ValueError("No valid data found for plotting!")
    
    # Positions closer together
    positions = [1, 1.5]
    box_width = 0.3

    bp = ax.boxplot(box_data, labels=labels, patch_artist=True,
                    positions=positions,
                    widths=box_width,
                    showfliers=True,
                    flierprops=dict(marker='D', markersize=2, alpha=0.5, markerfacecolor='gray'),
                    medianprops=dict(color='black', linewidth=1.5),
                    boxprops=dict(linewidth=1.2),
                    whiskerprops=dict(linewidth=1.2),
                    capprops=dict(linewidth=1.2))

    # Colors
    colors = ['#7fb3a3', '#f0e68c']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')

    # Log scale
    ax.set_yscale('log')
    ax.set_ylim(1, 1000) # you can change here and make it more generic 
    ax.set_yticks([10, 100])
    ax.set_yticklabels(['10¹', '10²'], fontsize=22)

    # Labels
    ax.set_xlabel('Cell line', fontsize=24, fontweight='bold')

    # Tick size
    ax.tick_params(axis='x', labelsize=22)
    ax.tick_params(axis='y', labelsize=22)

  
    ax.grid(False)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.2)
        spine.set_edgecolor('black')


    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.show()
    print(f"Plot saved to: {output_path}")



def main():
    # Define paths
    input_directory = "/human_vr_figure2" #update here
    output_directory = "your_output_file" #update here
    output_filename = "loop_size_distributions_boxplot.png"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    # Full output path
    output_path = os.path.join(output_directory, output_filename)
    
    try:
        # Read all data files
        print("Reading loop size analysis files...")
        combined_data = read_loop_data_files(input_directory)
        
        # Print summary statistics
        print("\nSummary Statistics:")
        for cell_type in combined_data['cell_type'].unique():
            subset = combined_data[combined_data['cell_type'] == cell_type]
            print(f"{cell_type}:")
            print(f"  - Number of loops: {len(subset)}")
            print(f"  - Loop size bins (log) - Mean: {subset['log_loop_size_bins'].mean():.3f}, "
                  f"Median: {subset['log_loop_size_bins'].median():.3f}")
            print(f"  - Loop size bins range: {subset['loop_size_bins'].min()} - {subset['loop_size_bins'].max()}")
        
        # Create the boxplot
        print(f"\nCreating boxplot...")
        create_boxplot(combined_data, output_path)
        
        # Save the combined data as well for reference
        data_output_path = os.path.join(output_directory, "combined_loop_data.csv")
        combined_data.to_csv(data_output_path, index=False)
        print(f"Combined data saved to: {data_output_path}")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()