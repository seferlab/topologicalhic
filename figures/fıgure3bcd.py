import csv
import glob
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Set global font sizes
plt.rcParams.update({
    'font.size': 22,        # Normal labels
    'axes.labelsize': 22,   # X and Y axis labels
    'xtick.labelsize': 16,  # X tick labels
    'ytick.labelsize': 16,  # Y tick labels
    'legend.fontsize': 22,  # Legend font size
})

# Define classes and their paths
classes = {
    "RUES_CM": "/VR/RUES_CM",
    "RUES_CP": "/VR/RUES_CP", 
    "RUES_ESC": "/VR/RUES_ESC",
    "RUES_Fetal_Heart": "/VR/RUES_Fetal_Heart",
    "RUES_MES": "/VR/RUES_MES"
}

#need to read in persistence diagrams
def readPersisFile(filename):
    persisdiag = []
    try:
        with open(filename,'rt') as f:
            freader = csv.reader(f,delimiter='\t')
            for line in freader:
                if len(line) < 3:
                    continue
                if int(line[0]) != 1: 
                    continue
                birth, death = float(line[1]), float(line[2])
                # Sadece geçerli persistence çiftlerini ekle (death > birth)
                if death > birth:
                    persisdiag.append((birth, death))
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
    return persisdiag

def find_files_in_classes(pattern, chrnum):
    """Find files matching pattern in all class directories for given chromosome"""
    found_files = []
    for class_name, class_path in classes.items():
        # Search for files in the class directory
        search_pattern = os.path.join(class_path, f"*{chrnum}{pattern}")
        files = glob.glob(search_pattern)
        for file in files:
            found_files.append((class_name, file))
    return found_files

def load_data_from_classes():
    """Load all data from class directories"""
    origpersis = {}
    distperms = {}
    randperms = {}
    edgeperms = {}
    
    for chrnum in range(13, 20):
        
        # Find original persistence diagrams
        orig_files = find_files_in_classes("_persisdiagram.txt", chrnum)
        orig_files_filtered = []
        for class_name, filepath in orig_files:
            # Skip files that contain perm patterns (these are permutation files)
            if 'distperm' in filepath or 'randperm' in filepath or 'edgeperm' in filepath:
                continue
            orig_files_filtered.append((class_name, filepath))
            
        for class_name, filepath in orig_files_filtered:
            persis_data = readPersisFile(filepath)
            if persis_data:  # Sadece veri varsa ekle
                if chrnum in origpersis:
                    origpersis[chrnum].append((class_name, persis_data))
                else:
                    origpersis[chrnum] = [(class_name, persis_data)]
        
        # Find distance permutation files
        dist_files = find_files_in_classes("_distperm*_persisdiagram.txt", chrnum)
        for class_name, filepath in dist_files:
            persis_data = readPersisFile(filepath)
            if persis_data:  # Sadece veri varsa ekle
                if chrnum in distperms:
                    distperms[chrnum].append((class_name, persis_data))
                else:
                    distperms[chrnum] = [(class_name, persis_data)]
        
        # Find random permutation files
        rand_files = find_files_in_classes("_randperm*_persisdiagram.txt", chrnum)
        for class_name, filepath in rand_files:
            persis_data = readPersisFile(filepath)
            if persis_data:  # Sadece veri varsa ekle
                if chrnum in randperms:
                    randperms[chrnum].append((class_name, persis_data))
                else:
                    randperms[chrnum] = [(class_name, persis_data)]
        
        # Find edge permutation files - daha geniş arama
        edge_files = find_files_in_classes("_edgeperm*_persisdiagram.txt", chrnum)
        # Alternatif arama desenleri dene
        if not edge_files:
            edge_files.extend(find_files_in_classes("*edgeperm*_persisdiagram.txt", chrnum))
        if not edge_files:
            edge_files.extend(find_files_in_classes("*edge*perm*_persisdiagram.txt", chrnum))
            
        for class_name, filepath in edge_files:
            persis_data = readPersisFile(filepath)
            if persis_data:  # Sadece veri varsa ekle
                if chrnum in edgeperms:
                    edgeperms[chrnum].append((class_name, persis_data))
                else:
                    edgeperms[chrnum] = [(class_name, persis_data)]
    
    return origpersis, distperms, randperms, edgeperms

#plot birth times
def plotBirthTimes(origpersis,randperms,edgeperms,distperms,filename):
    origtimes,randtimes,edgetimes,disttimes = [],[],[],[]
    
    for chrnum in range(13,20):
        if chrnum in origpersis:
            origtimes.extend([x[0] for celltype in origpersis[chrnum] for x in celltype[1]])
        if chrnum in randperms:
            randtimes.extend([x[0] for perm in randperms[chrnum] for x in perm[1]])
        if chrnum in edgeperms:
            edge_births = [x[0] for perm in edgeperms[chrnum] for x in perm[1]]
            edgetimes.extend(edge_births)
        if chrnum in distperms:
            disttimes.extend([x[0] for perm in distperms[chrnum] for x in perm[1]])
    
    if len(origtimes) == 0:
        print("No original data found, skipping birth times plot")
        return
        
    fig,ax = plt.subplots(figsize=(10, 8))

    # Plotting order değiştir - edge perm'i daha görünür yap
    if len(disttimes) > 0:
        sns.distplot(disttimes, label="Linear dependence", ax=ax, kde=True,kde_kws={'linewidth': 3},norm_hist=True,color='#ef5675',hist_kws={"edgecolor": "none","alpha":0.2})
    if len(randtimes) > 0:
        sns.distplot(randtimes, label="Random permutation", ax=ax, kde=True,kde_kws={'linewidth': 3}, norm_hist=True,color='#ffa600',hist_kws={"edgecolor": "none","alpha":0.2})
    if len(edgetimes) > 0:
        sns.distplot(edgetimes, label="Edge permutation", ax=ax, kde=True,kde_kws={'linewidth': 3},norm_hist=True,color='#7a5195',hist_kws={"edgecolor": "none","alpha":0.3})
    sns.distplot(origtimes, label="Hi-C", ax=ax, kde=True,kde_kws={'linewidth': 3}, norm_hist=True,color='#003f5c',hist_kws={"edgecolor": "none","alpha":0.2})
    
    ax.set_xlim([0,1])
    ax.set_yticklabels([])
    plt.xlabel('Loop birth time')
    ax.legend(frameon=False,loc='upper left')
    plt.savefig(filename, bbox_inches='tight', dpi=600)
    print('saved birth times histogram to',filename)

#plot lifespans
def plotLifeSpans(origpersis,randperms,edgeperms,distperms,filename):
    origtimes,randtimes,edgetimes,disttimes = [],[],[],[]
    
    for chrnum in range(13,20):
        if chrnum in origpersis:
            origtimes.extend([x[1]-x[0] for celltype in origpersis[chrnum] for x in celltype[1] if x[1] > x[0]])
        if chrnum in randperms:
            randtimes.extend([x[1]-x[0] for perm in randperms[chrnum] for x in perm[1] if x[1] > x[0]])
        if chrnum in edgeperms:
            edge_lifespans = [x[1]-x[0] for perm in edgeperms[chrnum] for x in perm[1] if x[1] > x[0]]
            edgetimes.extend(edge_lifespans)
        if chrnum in distperms:
            disttimes.extend([x[1]-x[0] for perm in distperms[chrnum] for x in perm[1] if x[1] > x[0]])
    
    if len(origtimes) == 0:
        print("No original data found, skipping life spans plot")
        return
        
    fig,ax = plt.subplots(figsize=(10, 8))

    # Plotting order değiştir - edge perm'i daha görünür yap
    if len(disttimes) > 0:
        sns.distplot(disttimes, label="Linear dependence", ax=ax, kde=True,kde_kws={'linewidth': 3},norm_hist=True,color='#ef5675',hist_kws={"edgecolor": "none","alpha":0.2})
    if len(randtimes) > 0:
        sns.distplot(randtimes, label="Random permutation", ax=ax, kde=True,kde_kws={'linewidth': 3},norm_hist=True,color = '#ffa600',hist_kws={"edgecolor": "none","alpha":0.2})
    if len(edgetimes) > 0:
        sns.distplot(edgetimes, label="Edge permutation", ax=ax, kde=True,kde_kws={'linewidth': 3},norm_hist=True,color='#7a5195',hist_kws={"edgecolor": "none","alpha":0.3})
    sns.distplot(origtimes, label="Hi-C", ax=ax, kde=True,kde_kws={'linewidth': 3},norm_hist=True,color='#003f5c',hist_kws={"edgecolor": "none","alpha":0.2})
    
    ax.set_xlim([0,0.6])
    ax.set_yticklabels([])
    plt.xlabel('Loop life span')
    ax.legend(frameon=False,loc='upper right')
    plt.savefig(filename, bbox_inches='tight', dpi=600)
    print('saved life spans histogram to',filename)

#plot numbers of loops (ratio w/ original)
def plotLoopNumbers(origpersis,randperms,edgeperms,distperms,filename):
    randtimes,edgetimes,disttimes = [],[],[]
    for chrnum in range(13,20):
        if chrnum not in origpersis:
            continue
            
        origloopnums = {}
        for celltypes in origpersis[chrnum]:
            origloopnums[celltypes[0]] = len(celltypes[1])
        
        if chrnum in randperms:
            for perm in randperms[chrnum]:
                celltype = perm[0]
                if celltype in origloopnums and origloopnums[celltype] > 0:
                    randtimes.append(float(len(perm[1]))/origloopnums[celltype])
        if chrnum in edgeperms:
            for perm in edgeperms[chrnum]:
                celltype = perm[0]
                if celltype in origloopnums and origloopnums[celltype] > 0:
                    edgetimes.append(float(len(perm[1]))/origloopnums[celltype])
        if chrnum in distperms:
            for perm in distperms[chrnum]:
                celltype = perm[0]
                if celltype in origloopnums and origloopnums[celltype] > 0:
                    disttimes.append(float(len(perm[1]))/origloopnums[celltype])
    
    if len(randtimes) == 0 and len(edgetimes) == 0 and len(disttimes) == 0:
        print("No permutation data found, skipping loop counts plot")
        return
        
    fig,ax = plt.subplots(figsize=(10, 8))

    if len(disttimes) > 0:
        sns.distplot(disttimes, label="Linear dependence", ax=ax, kde=False,norm_hist=True,color='#ef5675',hist_kws={"edgecolor": "none","alpha":0.75})
    if len(randtimes) > 0:
        sns.distplot(randtimes, label="Random permutation", ax=ax, kde=False,norm_hist=True,color='#ffa600',hist_kws={"edgecolor": "none","alpha":0.75})
    if len(edgetimes) > 0:
        sns.distplot(edgetimes, label="Edge permutation", ax=ax, kde=False,norm_hist=True,color='#7a5195',hist_kws={"edgecolor": "none","alpha":0.75})
    
    ax.set_yticklabels([])
    if len(edgetimes) > 0:
        plt.xticks(np.arange(0,np.ceil(max(edgetimes)),step=1))
    plt.xlabel('Ratio of loop count in model to corresponding Hi-C')
    ax.legend(frameon=False,loc='upper right')
    plt.savefig(filename, bbox_inches='tight', dpi=600)
    print('saved loop counts histogram to',filename)

# Main execution
if __name__ == "__main__":
    outputloc = "your_output_location" #write your own
    
    print("Loading data from class directories...")
    origpersis, distperms, randperms, edgeperms = load_data_from_classes()
    
    print("Creating high-resolution PNG plots...")
    plotBirthTimes(origpersis, randperms, edgeperms, distperms, outputloc + 'birthtimes.png')
    plotLifeSpans(origpersis, randperms, edgeperms, distperms, outputloc + 'lifespans.png') 
    plotLoopNumbers(origpersis, randperms, edgeperms, distperms, outputloc + 'loopcounts.png')
    
    print("Analysis complete!")