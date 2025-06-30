import numpy as np
import re
import os
import gudhi as gd
import time
import sys
import csv
import argparse
import networkx as nx
import math

def read_persisdiag(input_file):
    persisdiag = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue  # skip malformed lines
            dim = int(parts[0])
            birth = float(parts[1])
            death = float(parts[2])
            if (death-birth < 0.0000000001): continue
            one_simp = eval(parts[3])  # e.g., [617]
            pair_simp = eval(parts[4])  # e.g., [617, 462]
            #print(one_simp)
            persisdiag.append((dim, birth, death, (one_simp, pair_simp)))
    return persisdiag

def read_distmat(dist_file):
    dist_data = []
    with open(dist_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            row = []
            for val in parts:
                if val == "1" or val == "1.0":
                    row.append(1.0)
                else:
                    try:
                        row.append(float(val))
                    except:
                        row.append(1.0)  # fallback for malformed values
            dist_data.append(row)

    # Convert to a symmetric square matrix
    size = len(dist_data)
    matrix = np.ones((size, size))
    for i in range(size):
        for j in range(len(dist_data[i])):
            matrix[i, j] = dist_data[i][j]
            matrix[j, i] = dist_data[i][j]  # symmetry
    return matrix

def matrix_normalize(TAD_matrix_all):
    distance_matrix_all=[]
    for TAD_matrix in TAD_matrix_all:
        TAD_matrix=np.log(TAD_matrix+1)
        max_num = TAD_matrix.max()
        TAD_matrix = TAD_matrix/(1.01*max_num)
        for i in range(TAD_matrix.shape[0]):
            TAD_matrix[i,i]=1.0
        distance_matrix=1-TAD_matrix
        distance_matrix_all.append(distance_matrix)
    return distance_matrix_all

def generate_1_dim_simp_list_from_dist_mat(dist_matrix):
    mat_dim = dist_matrix.shape[0]
    dim_1_simp_list = [[i,j,dist_matrix[i,j]] for i in range(0,mat_dim-1) for j in range(i+1,mat_dim)]
    dim_1_simp_list.sort(key=lambda x: [x[2],max(x[0],x[1]),min(x[0],x[1])])
    return dim_1_simp_list

def get_one_dim_persis(fullpersis):
    persis_1_dim = []
    for dim,birth,death,simp_pair in fullpersis:
        if(dim==1):
            persis_1_dim.append([birth,death,simp_pair[0]])
    persis_1_dim.sort(key=lambda x: [x[0]-x[1]])
    return persis_1_dim

def calculate_loop_sizes_and_lifespans(dist_mat, dim_1_simp_list, persis_1_dim_list):
    """Calculate loop sizes and lifespans for each genomic loop"""
    loop_data = []
    
    for i in range(len(persis_1_dim_list)):
        birth_time = persis_1_dim_list[i][0]
        death_time = persis_1_dim_list[i][1]
        lifespan = death_time - birth_time
        
        essen_edge = persis_1_dim_list[i][2]
        essen_index = dim_1_simp_list.index([min(essen_edge[0],essen_edge[1]),max(essen_edge[0],essen_edge[1]),dist_mat[essen_edge[0],essen_edge[1]]])
        
        # Create graph without the essential edge to find alternative path
        G = nx.Graph()
        G.add_weighted_edges_from(dim_1_simp_list[:essen_index])
        
        try:
            # Find shortest path between the two nodes
            path = nx.shortest_path(G, source=essen_edge[0], target=essen_edge[1], weight="weight")
            loop_size_bins = len(path)  # Number of genomic bins in the loop
            
            loop_data.append({
                'loop_id': i,
                'birth_time': birth_time,
                'death_time': death_time,
                'lifespan': lifespan,
                'loop_size_bins': loop_size_bins,
                'path_nodes': path
            })
            
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            # If no path exists, use direct distance
            loop_size_bins = abs(essen_edge[1] - essen_edge[0]) + 1

            loop_data.append({
                'loop_id': i,
                'birth_time': birth_time,
                'death_time': death_time,
                'lifespan': lifespan,
                'loop_size_bins': loop_size_bins,
                'path_nodes': [essen_edge[0], essen_edge[1]]
            })
    
    return loop_data

def write_loop_analysis_file(loop_data, output_filename):
    """Write loop analysis data (only birth/death/lifespan/loop size)"""
    with open(output_filename, 'w') as f:
        # Write header
        f.write("loop_id\tbirth_time\tdeath_time\tlifespan\tloop_size_bins\tlog_loop_size_bins\n")
        
        # Write data
        for loop in loop_data:
            log_loop_size_bins = math.log10(max(1, loop['loop_size_bins']))
            f.write(f"{loop['loop_id']}\t{loop['birth_time']:.6f}\t{loop['death_time']:.6f}\t{loop['lifespan']:.6f}\t")
            f.write(f"{loop['loop_size_bins']}\t{log_loop_size_bins:.6f}\n")
    
    print(f'Wrote loop analysis data to {output_filename}')
    
    # Write summary statistics
    summary_filename = output_filename.replace('.txt', '_summary.txt')
    with open(summary_filename, 'w') as f:
        f.write("# Loop Analysis Summary Statistics\n")
        f.write(f"Total number of loops: {len(loop_data)}\n")
        
        lifespans = [loop['lifespan'] for loop in loop_data]
        loop_sizes = [loop['loop_size_bins'] for loop in loop_data]
        
        f.write(f"Lifespan range: {min(lifespans):.6f} to {max(lifespans):.6f}\n")
        f.write(f"Loop size (bins) range: {min(loop_sizes)} to {max(loop_sizes)}\n")
        f.write(f"Average lifespan: {np.mean(lifespans):.6f}\n")
        f.write(f"Average loop size (bins): {np.mean(loop_sizes):.2f}\n")
        
        # Lifespan bins for plotting
        f.write("\n# Suggested lifespan bins for plotting:\n")
        lifespan_bins = np.linspace(min(lifespans), max(lifespans), 10)
        for i, bin_edge in enumerate(lifespan_bins):
            f.write(f"Bin {i}: {bin_edge:.6f}\n")
    
    print(f'Wrote summary statistics to {summary_filename}')


def main(persis_file, dist_file, output_path, output_name):
    t0 = time.time()
    
    # Read persistance diagram
    print("Reading Persistance diagram...")
    persisdiag = read_persisdiag(persis_file)
    
    # Read distance matrix
    print("Reading distance matrix...")
    distance_matrix = read_distmat(dist_file)
    
    # Extract 1D persistence and calculate loop properties
    print("Analyzing loop sizes and lifespans...")
    persis_dim_1_list = get_one_dim_persis(persisdiag)
    dim_1_simp = generate_1_dim_simp_list_from_dist_mat(distance_matrix)
    
    # Calculate loop sizes and lifespans
    loop_data = calculate_loop_sizes_and_lifespans(distance_matrix, dim_1_simp, persis_dim_1_list)
    
    # Write analysis file
    output_filename = output_path + output_name + '_loop_size_analysis.txt'
    write_loop_analysis_file(loop_data, output_filename)
    
    print(f'Analysis complete. Total time (s): {time.time()-t0}')
    print(f'Found {len(loop_data)} genomic loops')
    print(f'Data written to: {output_filename}')
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-ip', type=str, help="input file of persistance diagram")
    parser.add_argument('-id', type=str, help="input file of distance matrix")
    parser.add_argument('-o', type=str, help='the name of output files')
    parser.add_argument('-p', type=str, help="the path of output files")
    
    args = parser.parse_args()

    main(args.ip, args.id, args.p, args.o)
