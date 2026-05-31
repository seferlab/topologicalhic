import numpy as np
import gudhi as gd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
import os
import glob

# Set font to DejaVu Sans globally (widely available)
plt.rcParams['font.family'] = 'DejaVu Sans'

# Cell type definitions - change here
class_list = [
    ("RUES", "RUES CM", "/VR/RUES_CM"),
    ("RUES", "RUES CP", "/VR/RUES_CP"),
    ("RUES", "RUES ESC", "/VR/RUES_ESC"),
    ("RUES", "RUES Fetal Heart", "/VR/RUES_Fetal_Heart"),
    ("RUES", "RUES MES", "/VR/RUES_MES"),
    ("WTC", "WTC CM", "/VR/WTC_CM"),
    ("WTC", "WTC CP", "/VR/WTC_CP"),
    ("WTC", "WTC MES", "/VR/WTC_MES"),
    ("WTC", "WTC PSC", "/VR/WTC_PSC"),
]

# read in persistence file: return persistence diagram as [[b,d],[b,d],...] for use in bottleneckdist
def readPersisFile(filename):
    persisdiag = []
    fullpersisinfo = []
    with open(filename,'rt') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            persisdiag.append((int(line[0]),(float(line[1]),float(line[2]))))
            pair = eval(line[3])
            fullpersisinfo.append([int(line[0]), float(line[1]), float(line[2]), pair])
    return persisdiag, fullpersisinfo

def findPersistenceFiles(directory):
    """Find all persistence diagram files in a directory, excluding specific patterns"""
    pattern = os.path.join(directory, "*_persisdiagram.txt")
    all_files = glob.glob(pattern)
    
    # Filter out files containing distperm, edgeperm, or randperm
    filtered_files = []
    for file in all_files:
        filename = os.path.basename(file)
        if not any(exclude in filename for exclude in ['distperm', 'edgeperm', 'randperm']):
            filtered_files.append(file)
    
    return filtered_files

def loadAllPersistenceData():
    """Load persistence data from all cell types and chromosomes"""
    all_data = {
        'h0_birth': [], 'h0_death': [], 'h0_lifespan': [],
        'h1_birth': [], 'h1_death': [], 'h1_lifespan': []
    }
    
    total_files_processed = 0
    
    for cell_group, cell_name, directory in class_list:
        print(f"Processing {cell_name}...")
        
        if not os.path.exists(directory):
            print(f"Warning: Directory {directory} does not exist. Skipping {cell_name}.")
            continue
        
        # Find all persistence files in this directory
        persis_files = findPersistenceFiles(directory)
        
        if not persis_files:
            print(f"Warning: No persistence files found in {directory}")
            continue
        
        files_processed = 0
        for file_path in persis_files:
            try:
                # Read persistence data
                persisdiag, _ = readPersisFile(file_path)
                
                # Process each persistence point
                for pt in persisdiag:
                    dim = pt[0]
                    birth = pt[1][0]
                    death = pt[1][1] if pt[1][1] != np.inf else 1.0  # Cap infinite deaths at 1.0
                    lifespan = death - birth
                    
                    if dim == 0:  # H0
                        all_data['h0_birth'].append(birth)
                        all_data['h0_death'].append(death)
                        all_data['h0_lifespan'].append(lifespan)
                    elif dim == 1:  # H1
                        all_data['h1_birth'].append(birth)
                        all_data['h1_death'].append(death)
                        all_data['h1_lifespan'].append(lifespan)
                
                files_processed += 1
                
            except Exception as e:
                print(f"Error processing file {file_path}: {e}")
        
        print(f"  Processed {files_processed} files for {cell_name}")
        total_files_processed += files_processed
    
    print(f"\nTotal files processed: {total_files_processed}")
    print(f"Total H0 features: {len(all_data['h0_birth'])}")
    print(f"Total H1 features: {len(all_data['h1_birth'])}")
    
    return all_data

def plotCombinedViolinPlot(all_data, filename):
    """Create a combined violin plot for all cell types and chromosomes"""
    
    # Check if we have data
    if not all_data['h0_birth'] and not all_data['h1_birth']:
        print("No data found to plot!")
        return
    
    # Create the violin plot
    fig, ax = plt.subplots(figsize=(12, 8))
    positions = [1, 2, 3, 4, 5, 6]
    
    # Prepare data for violin plot
    data = [
        all_data['h0_birth'], all_data['h1_birth'],
        all_data['h0_death'], all_data['h1_death'],
        all_data['h0_lifespan'], all_data['h1_lifespan']
    ]
    
    # Remove empty datasets and corresponding positions
    filtered_data = []
    filtered_positions = []
    for i, dataset in enumerate(data):
        if dataset:  # Only include non-empty datasets
            filtered_data.append(dataset)
            filtered_positions.append(positions[i])
    
    if not filtered_data:
        print("No valid data to plot!")
        return
    
    # Create violin plot
    violin_parts = ax.violinplot(filtered_data, filtered_positions, showmedians=True, widths=0.8)
    
    # Color the violin plots according to the figure (H0 orange, H1 blue)
    for i, pc in enumerate(violin_parts['bodies']):
        pos = filtered_positions[i]
        if pos in [1, 3, 5]:  # H0 positions (birth, death, lifespan for H0)
            pc.set_facecolor('#fc8d59')  # Orange for H0
            pc.set_alpha(0.7)
        else:  # H1 positions (birth, death, lifespan for H1)
            pc.set_facecolor('#91bfdb')  # Blue for H1
            pc.set_alpha(0.7)
    
    # Set the median line colors
    violin_parts['cmedians'].set_colors(['black'] * len(filtered_data))
    
    # Add labels - adjust based on available data
    if len(filtered_positions) == 6:
        plt.xticks([1.5, 3.5, 5.5], ['Birth time', 'Death time', 'Lifespan'], fontsize=20)
    else:
        # Create custom labels based on available positions
        label_positions = []
        labels = []
        if 1 in filtered_positions and 2 in filtered_positions:
            label_positions.append(1.5)
            labels.append('Birth time')
        if 3 in filtered_positions and 4 in filtered_positions:
            label_positions.append(3.5)
            labels.append('Death time')
        if 5 in filtered_positions and 6 in filtered_positions:
            label_positions.append(5.5)
            labels.append('Lifespan')
        
        if label_positions:
            plt.xticks(label_positions, labels, fontsize=20)
    
    plt.yticks(fontsize=20)
    plt.ylabel('Radius', fontsize=22)

    
    # Add legend with correct colors matching the figure
    import matplotlib.patches as mpatches
    h0_patch = mpatches.Patch(color='#fc8d59', label='H0')  # Orange for H0
    h1_patch = mpatches.Patch(color='#91bfdb', label='H1')  # Blue for H1
    plt.legend(handles=[h0_patch, h1_patch], title="Dimension", loc='upper right', fontsize=20, title_fontsize=22)
    
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight', dpi=600)
    plt.close()
    print(f'Saved combined violin plot to {filename}')

def main():
    output_dir = "your_output_location" #write your own
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print("Loading persistence data from all cell types and chromosomes...")
    all_data = loadAllPersistenceData()
    
    # Generate combined violin plot
    violin_file = os.path.join(output_dir, "combined_all_celltypes_violinplot.png")
    plotCombinedViolinPlot(all_data, violin_file)
    
    print(f"Analysis complete. Combined figure saved to {output_dir}")

if __name__ == "__main__":
    main()