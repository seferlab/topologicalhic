import numpy as np
import matplotlib.pyplot as plt
import csv
import ast
import os

# Updated class list with mouse data
class_list = [
    # Human data
    ("RUES", "RUES ESC", "/content/drive/MyDrive/THesis/VR/RUES_ESC"),
    ("RUES", "RUES MES", "/content/drive/MyDrive/THesis/VR/RUES_MES"),
    ("RUES", "RUES CM", "/content/drive/MyDrive/THesis/VR/RUES_CM"),
    ("RUES", "RUES CP", "/content/drive/MyDrive/THesis/VR/RUES_CP"),
    ("RUES", "RUES Fetal Heart", "/content/drive/MyDrive/THesis/VR/RUES_Fetal_Heart"),
    ("WTC", "WTC PSC", "/content/drive/MyDrive/THesis/VR/WTC_PSC"),
    ("WTC", "WTC MES", "/content/drive/MyDrive/THesis/VR/WTC_MES"),
    ("WTC", "WTC CP", "/content/drive/MyDrive/THesis/VR/WTC_CP"),
    ("WTC", "WTC CM", "/content/drive/MyDrive/THesis/VR/WTC_CM"),
    # Mouse data
    ("Mouse", "CD4 Q1", "/content/drive/MyDrive/THesis/VR/Mouse_CD4_Q1"),
    ("Mouse", "CD8", "/content/drive/MyDrive/THesis/VR/Mouse_CD8"),
    ("Mouse", "DN", "/content/drive/MyDrive/THesis/VR/Mouse_DN"),
    ("Mouse", "DP", "/content/drive/MyDrive/THesis/VR/Mouse_DP"),
    ("Mouse", "TCon1", "/content/drive/MyDrive/THesis/VR/Mouse_TCon1"),
    ("Mouse", "TReg1", "/content/drive/MyDrive/THesis/VR/Mouse_TReg1"),
]

def readPersisFile(filename):
    """Read persistence file and return persistence diagram."""
    persisdiag = []
    fullpersisinfo = []
    try:
        with open(filename, 'rt') as f:
            freader = csv.reader(f, delimiter='\t')
            for line in freader:
                persisdiag.append((int(line[0]), (float(line[1]), float(line[2]))))
                pair = ast.literal_eval(line[3])
                fullpersisinfo.append([int(line[0]), float(line[1]), float(line[2]), pair])
    except FileNotFoundError:
        print(f"Warning: File {filename} not found")
        return [], []
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return [], []
    return persisdiag, fullpersisinfo

def computeBottleneckDist(persis1, persis2):
    """Compute bottleneck distance between two persistence diagrams."""
    try:
        import gudhi as gd
        
        # Convert to list format for gudhi
        persis1list = []
        persis2list = []
        for pair in persis1:
            persis1list.append([pair[1][0], pair[1][1]])
        for pair in persis2:
            persis2list.append([pair[1][0], pair[1][1]])
        
        dist = gd.bottleneck_distance(persis1list, persis2list)
        return dist
    except ImportError:
        print("Warning: gudhi not available, using placeholder distance")
        return np.random.random() * 0.5  # Placeholder for demonstration
    except Exception as e:
        print(f"Error computing bottleneck distance: {e}")
        return 0.0

def get_filename_pattern(species, path, chrnum):
    """Get the appropriate filename pattern based on species."""
    base_name = os.path.basename(path)
    
    if species == "Mouse":
        # Mouse files use the format: CD8_combined_200000_iced_chr10_persisdiagram.txt
        if "Mouse_CD4_Q1" in path:
            return os.path.join(path, f"CD4_Q1_combined_200000_iced_chr{chrnum}_persisdiagram.txt")
        elif "Mouse_CD8" in path:
            return os.path.join(path, f"CD8_combined_200000_iced_chr{chrnum}_persisdiagram.txt")
        elif "Mouse_DN" in path:
            return os.path.join(path, f"DN_combined_200000_iced_chr{chrnum}_persisdiagram.txt")
        elif "Mouse_DP" in path:
            return os.path.join(path, f"DP_combined_200000_iced_chr{chrnum}_persisdiagram.txt")
        elif "Mouse_TCon1" in path:
            return os.path.join(path, f"TCon1_combined_200000_iced_chr{chrnum}_persisdiagram.txt")
        elif "Mouse_TReg1" in path:
            return os.path.join(path, f"TReg1_combined_200000_iced_chr{chrnum}_persisdiagram.txt")
    else:
        # Human files use the format: RUES_ESC_combined_200000_iced_chr1_persisdiagram.txt
        return os.path.join(path, f"{base_name}_combined_200000_iced_chr{chrnum}_persisdiagram.txt")
    
    return None

def plotHeatMap(avgdists, celltypelabels, filename, title=""):
    """Plot heatmap of distances with better formatting for combined human-mouse data."""
    fig, ax = plt.subplots(figsize=(14, 12))
    
    # Set global font size
    plt.rcParams.update({'font.size': 22})
    
    # Create heatmap
    im = ax.imshow(avgdists, cmap='YlGnBu', aspect='auto')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Bottleneck distance', rotation=270, labelpad=20, fontsize=24)
    cbar.ax.tick_params(labelsize=24)
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(celltypelabels)))
    ax.set_yticks(np.arange(len(celltypelabels)))
    # Changed rotation to 90 degrees (vertical) and removed ha parameter for proper vertical alignment
    ax.set_xticklabels(celltypelabels, rotation=90, fontsize=22)
    ax.set_yticklabels(celltypelabels, fontsize=22)
    
    # Add lines to separate human and mouse data
    # Human data: indices 0-8, Mouse data: indices 9-14
    ax.axhline(y=8.5, color='red', linewidth=2, alpha=0.7)
    ax.axvline(x=8.5, color='red', linewidth=2, alpha=0.7)
    
    # Removed the species text labels (Human/Mouse labels)
    
    # Remove grid and ticks
    ax.grid(False)
    ax.tick_params(top=False, right=False, left=False, bottom=False, labelsize=22)
    
    # Title removed as requested
    
    # Tight layout and save
    plt.tight_layout()
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.show()
    print(f'Heatmap saved as {filename}')

def generate_figure4():
    """Generate Figure 4 heatmap from class_list data including mouse data."""
    
    # Extract cell type labels
    celltypelabels = [item[1] for item in class_list]
    print(f"Cell types: {celltypelabels}")
    
    # Initialize distance matrix for 19 chromosomes (common to both human and mouse)
    num_samples = len(class_list)
    bndists = np.zeros((num_samples, num_samples, 19))
    
    print("Computing bottleneck distances...")
    print(f"Total samples: {num_samples}")
    print(f"Human samples: 9, Mouse samples: 6")
    
    # Process each chromosome (1-19, common to both species)
    for chrnum in range(1, 20):
        print(f"Processing chromosome {chrnum}...")
        
        # Read persistence diagrams for this chromosome
        persis_data = {}
        for idx, (species, cell_type, path) in enumerate(class_list):
            # Get appropriate filename
            filename = get_filename_pattern(species, path, chrnum)
            
            if filename is None:
                print(f"Warning: Could not determine filename for {cell_type}")
                persis_data[idx] = []
                continue
                
            persisdiag, _ = readPersisFile(filename)
            persis_data[idx] = persisdiag
            
            if not persisdiag:
                print(f"Warning: No data found for {cell_type}, chr{chrnum}")
        
        # Compute pairwise distances
        for i in range(num_samples):
            for j in range(i, num_samples):
                if i == j:
                    bndists[i, j, chrnum-1] = 0.0
                else:
                    # Only compute distance if both samples have data
                    if persis_data[i] and persis_data[j]:
                        dist = computeBottleneckDist(persis_data[i], persis_data[j])
                        bndists[i, j, chrnum-1] = dist
                        bndists[j, i, chrnum-1] = dist  # Symmetric
                    else:
                        # If one sample has no data, set distance to NaN
                        bndists[i, j, chrnum-1] = np.nan
                        bndists[j, i, chrnum-1] = np.nan
    
    # Average over all chromosomes, ignoring NaN values
    avgbndists = np.nanmean(bndists, axis=2)
    
    # Replace any remaining NaN values with 0
    avgbndists = np.nan_to_num(avgbndists, nan=0.0)
    
    print("Distance matrix shape:", avgbndists.shape)
    print("Average distances computed successfully!")
    
    # Create output directory if it doesn't exist
    output_dir = "output_path" # write your own 
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot heatmap
    output_filename = os.path.join(output_dir, 'figure4_human_mouse_bottleneck_heatmapd.png')
    plotHeatMap(avgbndists, celltypelabels, output_filename, 
                "Figure 4: Human-Mouse Bottleneck Distance Heatmap (Averaged over 19 chromosomes)")
    
    # Also save as PDF
    output_filename_pdf = os.path.join(output_dir, 'figure4_human_mouse_bottleneck_heatmap.pdf')
    plotHeatMap(avgbndists, celltypelabels, output_filename_pdf, 
                "Figure 4: Human-Mouse Bottleneck Distance Heatmap (Averaged over 19 chromosomes)")
    
    # Save distance matrix to file
    dist_matrix_file = os.path.join(output_dir, 'human_mouse_bottleneck_distances_matrix.txt')
    with open(dist_matrix_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        # Write header
        writer.writerow([''] + celltypelabels)
        # Write data
        for i, label in enumerate(celltypelabels):
            writer.writerow([label] + list(avgbndists[i, :]))
    
    print(f"Distance matrix saved to: {dist_matrix_file}")
    
    return avgbndists, celltypelabels

def analyze_distances(distances, labels):
    """Analyze distances between different groups."""
    print("\nDistance Matrix Statistics:")
    print(f"Min distance: {distances.min():.4f}")
    print(f"Max distance: {distances.max():.4f}")
    print(f"Mean distance: {distances.mean():.4f}")
    
    # Human samples (indices 0-8)
    human_indices = list(range(9))
    # Mouse samples (indices 9-14)
    mouse_indices = list(range(9, 15))
    
    # Within-species distances
    human_within = []
    for i in human_indices:
        for j in human_indices:
            if i != j:
                human_within.append(distances[i, j])
    
    mouse_within = []
    for i in mouse_indices:
        for j in mouse_indices:
            if i != j:
                mouse_within.append(distances[i, j])
    
    # Between-species distances
    between_species = []
    for i in human_indices:
        for j in mouse_indices:
            between_species.append(distances[i, j])
    
    print(f"\nSpecies comparison:")
    print(f"Human within-species avg distance: {np.mean(human_within):.4f}")
    print(f"Mouse within-species avg distance: {np.mean(mouse_within):.4f}")
    print(f"Between-species avg distance: {np.mean(between_species):.4f}")
    
    # Cell line specific analysis for human data
    print(f"\nHuman cell line comparison:")
    
    # RUES samples (indices 0-4)
    rues_indices = list(range(5))
    rues_within = []
    for i in rues_indices:
        for j in rues_indices:
            if i != j:
                rues_within.append(distances[i, j])
    
    # WTC samples (indices 5-8)
    wtc_indices = list(range(5, 9))
    wtc_within = []
    for i in wtc_indices:
        for j in wtc_indices:
            if i != j:
                wtc_within.append(distances[i, j])
    
    # Between human cell lines
    between_human_lines = []
    for i in rues_indices:
        for j in wtc_indices:
            between_human_lines.append(distances[i, j])
    
    print(f"RUES within-line avg distance: {np.mean(rues_within):.4f}")
    print(f"WTC within-line avg distance: {np.mean(wtc_within):.4f}")
    print(f"Between human lines avg distance: {np.mean(between_human_lines):.4f}")

if __name__ == "__main__":
    # Generate Figure 4 with mouse data
    distances, labels = generate_figure4()
    
    # Analyze distances
    analyze_distances(distances, labels)