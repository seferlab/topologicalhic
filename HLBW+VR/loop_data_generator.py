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
try:
    import matrix_optimization
except ImportError:
    print("Error: The 'matrix_optimization' module is not installed.")
    print("Please compile and install it using the provided setup.py script (e.g., 'pip install .')")
    sys.exit(1)

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
        distance_matrix_all.append(distance_matrix)
    return distance_matrix_all

def TDA_func(distance_matrix):
    print('distance matrix size =',distance_matrix.shape)
    rips_complex = gd.RipsComplex(distance_matrix=distance_matrix,max_edge_length=1.1)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
    print('done creating simplex tree')
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
            diag.index((1,(btime,dtime)))
            htype = 1
            fullpersinfo.append([htype, btime, dtime, pair])
    return fullpersinfo

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

def calculate_loop_sizes_and_lifespans(reso, minimal_bin, dist_mat, dim_1_simp_list, persis_1_dim_list):
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

def project_to_kernel(D, mu=1, maxits=100, tolconv=1.0e-6, toleigs=1.0e-5):
    """Stage I: Projects a distance matrix to a new distance matrix with a PSD kernel."""
    if np.max(D) == 0: return D # Avoid division by zero
    gamma = -mu / np.max(D)
    K = np.exp(gamma * D)
    n = K.shape[0]
    low_val = np.exp(-mu)
    high_val = 1
    diag_val = 1
    
    Y = (K + K.T) / 2
    U = np.zeros((n, n))
    iter_count = 0

    while iter_count < maxits:
        T = Y - U
        eigvals, eigvecs = np.linalg.eigh(T)
        pos_eigvals = eigvals > toleigs
        X = eigvecs[:, pos_eigvals] @ np.diag(eigvals[pos_eigvals]) @ eigvecs[:, pos_eigvals].T
        U = X - T
        
        if np.linalg.norm(Y - X, np.inf) / np.linalg.norm(Y, np.inf) <= tolconv:
            break
        
        Y = X
        np.fill_diagonal(Y, diag_val)
        Y = np.clip(Y, low_val, high_val)
        iter_count += 1
        
    C = np.log(np.maximum(Y, 1e-15)) / gamma # Add epsilon to avoid log(0)
    np.fill_diagonal(C, 0)
    D_new = (C + C.T) / 2
    
    return D_new

def hlwb_algorithm(D, n_projection=100):
    """The full HLWB algorithm to find the nearest metric distance matrix."""
    # Stage I: Get a high-quality initial approximation
    D_proj = project_to_kernel(D, maxits=10)
    
    # Optional heuristic improvement (requires the C++ implementation)
    X0 = matrix_optimization.heuristic_improve(D_proj.copy(), D.copy(), n_improve=1)

    # Stage II: Refine the solution using HLWB projections (the C++ implementation)
    X = matrix_optimization.hlwb_projection(X0.copy(), D.copy(), n_projection=n_projection)
    
    # Final cleanup
    X = (X + X.T) / 2
    np.fill_diagonal(X, 0)

    return X

def main(input_file, output_path, output_name, resol):
    t0 = time.time()
    
    # Read Hi-C data
    print("Reading Hi-C data...")
    resolut, freq_mat, min_bin, _ = read_raw_HiC_data_no_split(input_file, resol)
    
    # Generate distance matrix
    print("Generate distance matrix...")
    distance_matrix = matrix_normalize([freq_mat])[0]

    print("\n--- Correcting distance matrix using HLWB algorithm ---")
    distance_matrix = hlwb_algorithm(distance_matrix.copy())
    print("--- Correction complete ---")
    
    # Calculate persistent homology
    print("Calculate persistent homology...")
    full_persis_info = TDA_func(distance_matrix)
    
    # Extract 1D persistence and calculate loop properties
    print("Analyzing loop sizes and lifespans...")
    persis_dim_1_list = get_one_dim_persis(full_persis_info)
    dim_1_simp = generate_1_dim_simp_list_from_dist_mat(distance_matrix)
    
    # Calculate loop sizes and lifespans
    loop_data = calculate_loop_sizes_and_lifespans(resolut, min_bin, distance_matrix, dim_1_simp, persis_dim_1_list)
    
    # Write analysis file
    output_filename = output_path + output_name + '_loop_size_analysis.txt'
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
