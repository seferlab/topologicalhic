import numpy as np
import gudhi as gd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv

# Set font to DejaVu Sans globally (widely available)
plt.rcParams['font.family'] = 'DejaVu Sans'
 
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

def plotBarcode(persis, filename):
    colors=['#fc8d59', '#91bfdb']  # Orange for H0, Blue for H1
    maxy = len(persis)
    # Sort by homology dimension (H0 first, then H1)
    persissorted = sorted(persis, key=lambda x: x[0])
    
    plt.figure(figsize=(12, 8))  # Standardized figure size
    
    for idx, bar in enumerate(persissorted):
        if bar[1][1] == np.inf:  # For features that never die (persistent forever)
            plt.plot([0,1], [maxy-idx, maxy-idx], c=colors[bar[0]])
        else:
            plt.plot(bar[1], [maxy-idx, maxy-idx], c=colors[bar[0]])
            
    # Format and save the plot
    plt.axis([0,1,0,idx])
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.xlabel('Radius', fontsize=20)
    
    # Create legend with the correct color order to match figure
    import matplotlib.patches as mpatches
    h0_patch = mpatches.Patch(color='#fc8d59', label='H0')  # Orange for H0
    h1_patch = mpatches.Patch(color='#91bfdb', label='H1')  # Blue for H1
    plt.legend(handles=[h0_patch, h1_patch], frameon=False, fontsize=22)
    
    plt.savefig(filename, bbox_inches='tight', dpi=600)
    plt.close()
    print('saved barcode plot to', filename)

def plotPersisDiagram(persis, filename):
    colors=['#fc8d59', '#91bfdb']  # Orange for H0, Blue for H1
    
    plt.figure(figsize=(12, 8))  # Standardized figure size (changed from 8,8)
    
    for pt in persis:
        if pt[1][1] == np.inf:  # For features that never die
            plt.scatter(0, 1, c=colors[pt[0]], s=50)  # Increased point size from 5 to 50
        else:
            plt.scatter(pt[1][0], pt[1][1], c=colors[pt[0]], s=50)  # Increased point size from 5 to 50
            
    plt.plot([0,1], [0,1], c='k')  # Diagonal line
    plt.axis([-0.05, 1.05, -0.05, 1.05])
    
    # Set tick parameters
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
    # Create legend with the correct color order to match figure
    import matplotlib.patches as mpatches
    h0_patch = mpatches.Patch(color='#fc8d59', label='H0')  # Orange for H0
    h1_patch = mpatches.Patch(color='#91bfdb', label='H1')  # Blue for H1
    plt.legend(handles=[h0_patch, h1_patch], frameon=False, fontsize=22)
    
    plt.xlabel('birth radius', fontsize=20)
    plt.ylabel('death radius', fontsize=20)
    plt.savefig(filename, bbox_inches='tight', dpi=600)
    plt.close()
    print('saved persistence diagram to', filename)

def plotViolinPlot(persis, filename):
    # Extract birth, death, and lifespan data for both dimensions
    h0_birth = []
    h0_death = []
    h0_lifespan = []
    h1_birth = []
    h1_death = []
    h1_lifespan = []
    
    for pt in persis:
        dim = pt[0]
        birth = pt[1][0]
        death = pt[1][1] if pt[1][1] != np.inf else 1.0  # Cap infinite deaths at 1.0
        lifespan = death - birth
        
        if dim == 0:  # H0
            h0_birth.append(birth)
            h0_death.append(death)
            h0_lifespan.append(lifespan)
        else:  # H1
            h1_birth.append(birth)
            h1_death.append(death)
            h1_lifespan.append(lifespan)
    
    # Create the violin plot
    fig, ax = plt.subplots(figsize=(12, 8))  # Standardized figure size
    positions = [1, 2, 3, 4, 5, 6]
    
    data = [h0_birth, h1_birth, h0_death, h1_death, h0_lifespan, h1_lifespan]
    violin_parts = ax.violinplot(data, positions, showmedians=True)
    
    # Color the violin plots according to the figure (H0 orange, H1 blue)
    for i, pc in enumerate(violin_parts['bodies']):
        if i % 2 == 0:  # H0 (positions 0, 2, 4)
            pc.set_facecolor('#fc8d59')  # Orange for H0
            pc.set_alpha(0.7)
        else:  # H1 (positions 1, 3, 5)
            pc.set_facecolor('#91bfdb')  # Blue for H1
            pc.set_alpha(0.7)
    
    # Add labels
    plt.xticks([1.5, 3.5, 5.5], ['Birth time', 'Death time', 'Lifespan'], fontsize=20)
    plt.yticks(fontsize=18)
    plt.ylabel('Radius', fontsize=22)
    
    # Add legend with correct colors matching the figure
    import matplotlib.patches as mpatches
    h0_patch = mpatches.Patch(color='#fc8d59', label='H0')  # Orange for H0
    h1_patch = mpatches.Patch(color='#91bfdb', label='H1')  # Blue for H1
    plt.legend(handles=[h0_patch, h1_patch], title="Dimension", loc='upper right', fontsize=20)
    
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight', dpi=600)
    plt.close()
    print('saved violin plot to', filename)

def main():
    # Settings for WTC11 CP cell, chromosome 5
    input_file = "/VR/WTC_CP/WTC_CP_combined_200000_iced_chr5_persisdiagram.txt"  # Update this path
    output_dir = "your_output_location"  # Update this path
    
    # Create output directory if it doesn't exist
    import os
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Read persistence diagram data
    persisdiag, fullpersisinfo = readPersisFile(input_file)
    
    # Generate plots
    barcode_file = output_dir + "WTC11_CP_chr5_barcodeplot.png"
    persis_file = output_dir + "WTC11_CP_chr5_persistenceplot.png"
    violin_file = output_dir + "WTC11_CP_chr5_violinplot.png"
    
    plotBarcode(persisdiag, barcode_file)
    plotPersisDiagram(persisdiag, persis_file)
    plotViolinPlot(persisdiag, violin_file)
    
    print("Analysis complete. Figures saved to", output_dir)

if __name__ == "__main__":
    main()