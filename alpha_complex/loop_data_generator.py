import numpy as np
import re
import os
import gudhi as gd
import time
import sys
import csv
import argparse
import networkx as nx
from sklearn.manifold import MDS
from sklearn.metrics import euclidean_distances
import math

def read_raw_HiC_data_no_split(file,reso):
    resolution=int(reso)
    print('resolution: ',resolution)
    all_data_list=[]
    with open(file,"r")as f:
        for line in f:
            line = line.strip()
            line = re.sub(r'\s+',' ', line)
            line = line.split(' ')
            all_data_list.append([float(line[i]) for i in range(len(line))])
    raw_data_matrix=np.array(all_data_list)
    temp_min=min(raw_data_matrix[:,0].min(),raw_data_matrix[:,1].min())
    temp_max=max(raw_data_matrix[:,0].max(),raw_data_matrix[:,1].max())
    dim=int((temp_max-temp_min)/resolution+1)
    data_frequency_matrix=np.zeros((dim,dim))
    for x,y,freq in all_data_list:
        data_frequency_matrix[int((x-temp_min)/resolution),int((y-temp_min)/resolution)]=freq
        data_frequency_matrix[int((y-temp_min)/resolution),int((x-temp_min)/resolution)]=freq
    return (resolution,data_frequency_matrix,temp_min,temp_max)

def matrix_normalize(TAD_matrix_all):
    distance_matrix_all=[]
    for TAD_matrix in TAD_matrix_all:
        TAD_matrix=np.log(TAD_matrix+1)
        max_num = TAD_matrix.max()
        TAD_matrix = TAD_matrix/(1.01*max_num)
        for i in range(TAD_matrix.shape[0]):
            TAD_matrix[i,i]=1.0
        distance_matrix=1-TAD_matrix
        
        # Ensure symmetry and non-negative diagonal
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        np.fill_diagonal(distance_matrix, 0.0)  # Distance to self should be 0
        
        distance_matrix_all.append(distance_matrix)
    return distance_matrix_all

def apply_mds_embedding(distance_matrix):
    """Apply MDS embedding to distance matrix"""
    print("Applying MDS embedding...")
    
    mds = MDS(
        n_components=2,
        dissimilarity='precomputed',
        random_state=42,
        metric=True,
        max_iter=3500,
        eps=1e-12
    )
    
    embedding = mds.fit_transform(distance_matrix)
    
    return embedding

def TDA_func(mds_embedding):
    """Apply TDA using pre-computed MDS embedding"""
    
    # Create Alpha Complex from MDS embedding
    alpha_complex = gd.AlphaComplex(points=mds_embedding)
    simplex_tree = alpha_complex.create_simplex_tree()
    print('done creating Alpha Complex simplex tree')
    
    # Get persistence
    diag = simplex_tree.persistence()
    
    pairs = simplex_tree.persistence_pairs()
    fullpersinfo = []
    for pair in pairs:
        btime = simplex_tree.filtration(pair[0])
        dtime = simplex_tree.filtration(pair[1])
        try:
            diag.index((0,(btime,dtime)))
            htype = 0
            fullpersinfo.append([htype, btime, dtime, pair])
        except:
            try:
                diag.index((1,(btime,dtime)))
                htype = 1
                fullpersinfo.append([htype, btime, dtime, pair])
            except:
                pass
    
    return fullpersinfo

def get_one_dim_persis(fullpersis):
    persis_1_dim = []
    for dim,birth,death,simp_pair in fullpersis:
        if(dim==1):
            persis_1_dim.append([birth,death,simp_pair[0]])
    persis_1_dim.sort(key=lambda x: [x[0]-x[1]])
    return persis_1_dim

def generate_1_dim_simp_list_from_embedding(mds_embedding):
    """Generate 1D simplex list from MDS embedding using Euclidean distances"""
    n_points = mds_embedding.shape[0]
    dim_1_simp_list = []
    
    for i in range(n_points-1):
        for j in range(i+1, n_points):
            dist = np.linalg.norm(mds_embedding[i] - mds_embedding[j])
            dim_1_simp_list.append([i, j, dist])
    
    dim_1_simp_list.sort(key=lambda x: [x[2], max(x[0], x[1]), min(x[0], x[1])])
    return dim_1_simp_list

def calculate_loop_sizes_and_lifespans(reso, minimal_bin, mds_embedding, persis_1_dim_list):
    """Calculate loop sizes and lifespans for each genomic loop using Alpha Complex"""
    loop_data = []
    
    # Generate simplex list from MDS embedding
    dim_1_simp_list = generate_1_dim_simp_list_from_embedding(mds_embedding)
    
    for i in range(len(persis_1_dim_list)):
        birth_time = persis_1_dim_list[i][0]
        death_time = persis_1_dim_list[i][1]
        lifespan = death_time - birth_time
        
        essen_edge = persis_1_dim_list[i][2]
        
        # Find the essential edge in the simplex list
        try:
            essen_index = next(idx for idx, edge in enumerate(dim_1_simp_list) 
                             if (edge[0] == min(essen_edge) and edge[1] == max(essen_edge)))
        except StopIteration:
            # If exact edge not found, use direct distance
            loop_size_bins = abs(essen_edge[1] - essen_edge[0]) + 1
            genomic_size = abs(int(minimal_bin + reso * essen_edge[1]) - int(minimal_bin + reso * essen_edge[0]))
            
            loop_data.append({
                'loop_id': i,
                'birth_time': birth_time,
                'death_time': death_time,
                'lifespan': lifespan,
                'loop_size_bins': loop_size_bins,
                'genomic_size': genomic_size,
                'genomic_start': int(minimal_bin + reso * min(essen_edge)),
                'genomic_end': int(minimal_bin + reso * max(essen_edge)),
                'path_nodes': [essen_edge[0], essen_edge[1]]
            })
            continue
        
        # Create graph without the essential edge to find alternative path
        G = nx.Graph()
        G.add_weighted_edges_from(dim_1_simp_list[:essen_index])
        
        try:
            # Find shortest path between the two nodes
            path = nx.shortest_path(G, source=essen_edge[0], target=essen_edge[1], weight="weight")
            loop_size_bins = len(path)  # Number of genomic bins in the loop
            
            # Convert to genomic coordinates
            loop_start = int(minimal_bin + reso * min(path))
            loop_end = int(minimal_bin + reso * max(path))
            genomic_size = loop_end - loop_start
            
            loop_data.append({
                'loop_id': i,
                'birth_time': birth_time,
                'death_time': death_time,
                'lifespan': lifespan,
                'loop_size_bins': loop_size_bins,
                'genomic_size': genomic_size,
                'genomic_start': loop_start,
                'genomic_end': loop_end,
                'path_nodes': path
            })
            
        except nx.NetworkXNoPath:
            # If no path exists, use direct distance
            loop_size_bins = abs(essen_edge[1] - essen_edge[0]) + 1
            genomic_size = abs(int(minimal_bin + reso * essen_edge[1]) - int(minimal_bin + reso * essen_edge[0]))
            
            loop_data.append({
                'loop_id': i,
                'birth_time': birth_time,
                'death_time': death_time,
                'lifespan': lifespan,
                'loop_size_bins': loop_size_bins,
                'genomic_size': genomic_size,
                'genomic_start': int(minimal_bin + reso * min(essen_edge)),
                'genomic_end': int(minimal_bin + reso * max(essen_edge)),
                'path_nodes': [essen_edge[0], essen_edge[1]]
            })
    
    return loop_data

def write_loop_analysis_file(loop_data, output_filename):
    """Write comprehensive loop analysis data for plotting"""
    with open(output_filename, 'w') as f:
        # Write header
        f.write("loop_id\tbirth_time\tdeath_time\tlifespan\tloop_size_bins\tgenomic_size\tgenomic_start\tgenomic_end\tlog_loop_size_bins\tlog_genomic_size\n")
        
        # Write data
        for loop in loop_data:
            log_loop_size_bins = math.log10(max(1, loop['loop_size_bins']))
            log_genomic_size = math.log10(max(1, loop['genomic_size']))
            
            f.write(f"{loop['loop_id']}\t{loop['birth_time']:.6f}\t{loop['death_time']:.6f}\t{loop['lifespan']:.6f}\t")
            f.write(f"{loop['loop_size_bins']}\t{loop['genomic_size']}\t{loop['genomic_start']}\t{loop['genomic_end']}\t")
            f.write(f"{log_loop_size_bins:.6f}\t{log_genomic_size:.6f}\n")
    
    print(f'Wrote loop analysis data to {output_filename}')
    
    # Write summary statistics
    summary_filename = output_filename.replace('.txt', '_summary.txt')
    with open(summary_filename, 'w') as f:
        f.write("# Loop Analysis Summary Statistics\n")
        f.write(f"Total number of loops: {len(loop_data)}\n")
        
        if loop_data:
            lifespans = [loop['lifespan'] for loop in loop_data]
            loop_sizes = [loop['loop_size_bins'] for loop in loop_data]
            genomic_sizes = [loop['genomic_size'] for loop in loop_data]
            
            f.write(f"Lifespan range: {min(lifespans):.6f} to {max(lifespans):.6f}\n")
            f.write(f"Loop size (bins) range: {min(loop_sizes)} to {max(loop_sizes)}\n")
            f.write(f"Genomic size range: {min(genomic_sizes)} to {max(genomic_sizes)}\n")
            f.write(f"Average lifespan: {np.mean(lifespans):.6f}\n")
            f.write(f"Average loop size (bins): {np.mean(loop_sizes):.2f}\n")
            f.write(f"Average genomic size: {np.mean(genomic_sizes):.2f}\n")
            
            # Lifespan bins for plotting
            f.write("\n# Suggested lifespan bins for plotting:\n")
            lifespan_bins = np.linspace(min(lifespans), max(lifespans), 10)
            for i, bin_edge in enumerate(lifespan_bins):
                f.write(f"Bin {i}: {bin_edge:.6f}\n")
    
    print(f'Wrote summary statistics to {summary_filename}')

def main(input_file, output_path, output_name, resol):
    t0 = time.time()
    
    # Read Hi-C data
    print("Reading Hi-C data...")
    resolut, freq_mat, min_bin, _ = read_raw_HiC_data_no_split(input_file, resol)
    
    # Generate distance matrix
    print("Generate distance matrix...")
    distance_matrix = matrix_normalize([freq_mat])[0]
    
    # Apply MDS embedding 
    mds_embedding = apply_mds_embedding(distance_matrix)
    
    # Calculate persistent homology using Alpha Complex
    print("Calculate persistent homology with Alpha Complex...")
    full_persis_info = TDA_func(mds_embedding)  
    
    # Extract 1D persistence and calculate loop properties
    print("Analyzing loop sizes and lifespans...")
    persis_dim_1_list = get_one_dim_persis(full_persis_info)
    
    # Calculate loop sizes and lifespans using Alpha Complex approach
    loop_data = calculate_loop_sizes_and_lifespans(resolut, min_bin, mds_embedding, persis_dim_1_list)
    
    # Write analysis file
    output_filename = output_path + output_name + '_alpha_loop_size_analysis.txt'
    write_loop_analysis_file(loop_data, output_filename)
    
    print(f'Analysis complete. Total time (s): {time.time()-t0}')
    print(f'Found {len(loop_data)} genomic loops')
    print(f'Data written to: {output_filename}')
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help="input file of HiC contact matrix")
    parser.add_argument('-o', type=str, help='the name of output files')
    parser.add_argument('-p', type=str, help="the path of output files")
    parser.add_argument('-r', type=str, help="resolution of HiC input file")
    
    args = parser.parse_args()
    main(args.i, args.p, args.o, args.r)
