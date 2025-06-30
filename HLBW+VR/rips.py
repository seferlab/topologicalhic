import numpy as np
import re
import os
import gudhi as gd
import time
import sys
import csv
import argparse
import networkx as nx
try:
    import matrix_optimization
except ImportError:
    print("Error: The 'matrix_optimization' module is not installed.")
    print("Please compile and install it using the provided setup.py script (e.g., 'pip install .')")
    sys.exit(1)

def read_raw_HiC_data(file):
    resolution=re.split('[_.]',os.path.basename(file).strip())[1]
    if(resolution[-2:]=='kb'):
        resolution=int(resolution[:-2])*1000
    elif(resolution[-2:]=='mb'):
        resolution=int(resolution[:-2])*1000000
    print(resolution)
    all_data_list=[]
    with open(file,"r")as f:
        for line in f:
            line = line.strip()
            line = re.sub(r'\s+',' ', line)
            line = line.split(' ')
            all_data_list.append([float(line[i]) for i in range(len(line))])
    raw_data_matrix=np.array(all_data_list)
    dim=int(max(raw_data_matrix[:,0].max(),raw_data_matrix[:,1].max())/resolution+1)
    data_frequency_matrix=np.zeros((dim,dim))
    for x,y,freq in all_data_list:
        data_frequency_matrix[int(x/resolution),int(y/resolution)]=freq
        data_frequency_matrix[int(y/resolution),int(x/resolution)]=freq
    return (resolution,data_frequency_matrix)

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
        # Avoid division by zero if max_num is 0
        max_num = TAD_matrix.max()
        if max_num > 0:
            TAD_matrix = TAD_matrix / (1.01 * max_num)
        for i in range(TAD_matrix.shape[0]):
            TAD_matrix[i,i]=1.0
        distance_matrix=1-TAD_matrix
        distance_matrix_all.append(distance_matrix)
    return distance_matrix_all

def TDA_func(distance_matrix, persfilename):
    print(f'Running TDA on matrix of size {distance_matrix.shape}...')
    rips_complex = gd.RipsComplex(distance_matrix=distance_matrix,max_edge_length=1.1)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
    print('...done creating simplex tree')
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
    writePersistencePairsToFile(fullpersinfo,persfilename)

def writePersistencePairsToFile(perspairs, filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f, delimiter='\t')
        for pers in perspairs:
            fwriter.writerow(pers)
    print(f'Wrote persistence pairs to {filename}')
    return

def randomlyPermuteDistMat(distance_matrix,flag=''):
    distmatsize = distance_matrix.shape
    n = distmatsize[0]
    permmat = np.zeros(distmatsize)
    if flag == 'edge':
        randperm = np.random.permutation(n)
        for rownum in range(n):
            permmat[randperm[rownum], :] = distance_matrix[rownum,randperm]
    elif flag == 'rand':
        permidx = np.triu_indices(n,1)
        alldistvals = distance_matrix[permidx]
        permvals = np.random.permutation(alldistvals)
        permmat[permidx] = permvals
        permidx_lower = np.tril_indices(n,-1)
        permmat[permidx_lower] = permmat.T[permidx_lower]
    elif flag == 'dist':
        for diagnum in range(n-1):
            diagvals = np.diag(distance_matrix,diagnum+1)
            avgval = np.mean(diagvals)
            std = np.std(diagvals)
            newdiag = np.random.normal(avgval,std,len(diagvals))
            newdiag = np.maximum(newdiag, 0)
            diagmat = np.diag(newdiag,diagnum+1)
            permmat += diagmat
            permmat += diagmat.T
    return permmat

def writeDistMatToFile(distmat,filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f,delimiter='\t')
        for row in distmat:
            fwriter.writerow(row)
    print(f'Wrote distance matrix to {filename}')
    return

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

def main(input_file, output_path, output_name, resol, run_perms, correct_matrix):
    t0 = time.time()
    
    # 1. Read and normalize data to get the initial distance matrix
    _, freq_mat, _, _ = read_raw_HiC_data_no_split(input_file, resol)
    print("Generating initial distance matrix...")
    distance_matrix = matrix_normalize([freq_mat])[0]

    if correct_matrix:
        print("\n--- Correcting distance matrix using HLWB algorithm ---")
        distance_matrix = hlwb_algorithm(distance_matrix.copy())
        print("--- Correction complete ---")
    
    writeDistMatToFile(distance_matrix, output_path + output_name + '_distmat.txt')
    
    print(f"\n--- Calculating Persistent Homology on {('Corrected' if correct_matrix else 'Original')} Matrix ---")
    persis_filename = output_path + output_name + '_persisdiagram.txt'
    TDA_func(distance_matrix, persis_filename)

    if run_perms:
        print(f"\n--- Running Permutation Analysis on {('Corrected' if correct_matrix else 'Original')} Matrix ---")
        for method in ["edge", "rand", "dist"]:
            print(f"\nPermutation method: '{method}'")
            permmat = randomlyPermuteDistMat(distance_matrix, method)
            perm_persis_filename = output_path + output_name + "_" + method + 'perm_persisdiagram.txt'
            TDA_func(permmat, perm_persis_filename)
    else:
        print("\n--- Skipping Permutation Analysis ---")

    print('\nTotal time (s):', time.time()-t0)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Hi-C Data TDA Analysis with Optional Metric Correction.")

    parser.add_argument('-i', type=str, required=True, help="Input file of HiC contact matrix")
    parser.add_argument('-o', type=str, required=True, help='The base name for output files')
    parser.add_argument('-p', type=str, required=True, help="The path for output files")
    parser.add_argument('-r', type=str, required=True, help="Resolution of HiC input file")
    parser.add_argument(
        '-c', '--correct_matrix',
        action='store_true',
        help="Enable metric correction using the HLWB algorithm. If not set, analysis runs on the original matrix."
    )
    parser.add_argument(
        '-perm', '--run_permutations',
        action='store_true',
        help="Enable the permutation analysis. If not set, this step is skipped."
    )

    args = parser.parse_args()
    
    if not os.path.exists(args.p):
        os.makedirs(args.p)
        print(f"Created output directory: {args.p}")

    main(args.i, args.p, args.o, args.r, args.run_permutations, args.correct_matrix)