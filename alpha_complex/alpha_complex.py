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
            #print(line)
            all_data_list.append([float(line[i]) for i in range(len(line))])
    raw_data_matrix=np.array(all_data_list)
    dim=int(max(raw_data_matrix[:,0].max(),raw_data_matrix[:,1].max())/resolution+1)
    #print(dim)
    data_frequency_matrix=np.zeros((dim,dim))
    for x,y,freq in all_data_list:
        data_frequency_matrix[int(x/resolution),int(y/resolution)]=freq
        data_frequency_matrix[int(y/resolution),int(x/resolution)]=freq
    #print(data_frequency_matrix)
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
            #print(line)
            all_data_list.append([float(line[i]) for i in range(len(line))])
    raw_data_matrix=np.array(all_data_list)
    temp_min=min(raw_data_matrix[:,0].min(),raw_data_matrix[:,1].min())
    temp_max=max(raw_data_matrix[:,0].max(),raw_data_matrix[:,1].max())
    dim=int((temp_max-temp_min)/resolution+1)
    #print('Hi-C matrix size =',dim)
    data_frequency_matrix=np.zeros((dim,dim))
    for x,y,freq in all_data_list:
        data_frequency_matrix[int((x-temp_min)/resolution),int((y-temp_min)/resolution)]=freq
        data_frequency_matrix[int((y-temp_min)/resolution),int((x-temp_min)/resolution)]=freq
    #print(data_frequency_matrix)
    return (resolution,data_frequency_matrix,temp_min,temp_max)

def split_TAD(freq_matrix,TAD_result_file,resol):
    TAD_matrix_list=[]
    with open(TAD_result_file,"r")as f:
        for line in f:
            line = line.strip()
            line = re.sub(r'\s+',' ', line)
            line = line.split(' ')
            #print([line[1],line[2]])
            index_x=int(int(line[1])/resol)
            index_y=int(int(line[2])/resol)+1
            TAD_matrix_list.append(freq_matrix[index_x:index_y,index_x:index_y])
            #print(freq_matrix[index_x:index_y,index_x:index_y])
            print(freq_matrix[index_x:index_y,index_x:index_y].shape)
    return TAD_matrix_list

def matrix_normalize(TAD_matrix_all):
    distance_matrix_all=[]
    for TAD_matrix in TAD_matrix_all:
        TAD_matrix=np.log(TAD_matrix+1)
        max_num = TAD_matrix.max()
        #print(max_num)
        TAD_matrix = TAD_matrix/(1.01*max_num)
        #print(TAD_matrix.shape[0])
        for i in range(TAD_matrix.shape[0]):
            TAD_matrix[i,i]=1.0
        #print(TAD_matrix)
        distance_matrix=1-TAD_matrix
        
        # Ensure symmetry and non-negative diagonal
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        np.fill_diagonal(distance_matrix, 0.0)  # Distance to self should be 0
        
        distance_matrix_all.append(distance_matrix)
        #print(distance_matrix)
    return distance_matrix_all

def calculate_lifespan_stats(diag):
    """Calculate min and max lifespan values from persistence diagram"""
    lifespans = []
    for dim, (birth, death) in diag:
        if death != float('inf'):  # Ignore infinite persistence
            lifespan = death - birth
            lifespans.append(lifespan)
    
    if lifespans:
        min_lifespan = min(lifespans)
        max_lifespan = max(lifespans)
        print(f"Min lifespan: {min_lifespan:.6f}")
        print(f"Max lifespan: {max_lifespan:.6f}")
        return min_lifespan, max_lifespan
    else:
        print("No finite lifespans found")
        return None, None

def calculate_stress(original_distances, embedded_distances):
    """Calculate stress value for MDS embedding"""
    stress = np.sum((original_distances - embedded_distances)**2) / np.sum(original_distances**2)
    print(f"Stress value: {stress:.6f}")
    return stress

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
    
    # Calculate stress using original distances
    embedded_distances = euclidean_distances(embedding)
    stress = calculate_stress(distance_matrix, embedded_distances)
    
    return embedding, stress

def TDA_func(distance_matrix, persfilename):
    print('distance matrix size =', distance_matrix.shape)
    
    # Apply MDS embedding first
    mds_embedding, stress = apply_mds_embedding(distance_matrix)
    
    # Create Alpha Complex from MDS embedding
    alpha_complex = gd.AlphaComplex(points=mds_embedding)
    simplex_tree = alpha_complex.create_simplex_tree()
    print('done creating Alpha Complex simplex tree')
    
    # Get persistence
    diag = simplex_tree.persistence()
    
    # Calculate lifespan statistics
    min_lifespan, max_lifespan = calculate_lifespan_stats(diag)
    
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
                # Handle case where pair is not found in either dimension
                pass
    
    writePersistencePairsToFile(fullpersinfo, persfilename)
    
    return (diag, fullpersinfo, mds_embedding, stress, min_lifespan, max_lifespan)

def writePersistencePairsToFile(perspairs, filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f, delimiter='\t')
        for pers in perspairs:
            fwriter.writerow(pers)
    print('wrote persistence pairs to',filename)
    return

def randomlyPermuteDistMat(distance_matrix,flag=''):
    distmatsize = distance_matrix.shape
    n = distmatsize[0]
    permmat = np.zeros(distmatsize)
    if flag == 'edge':
        # randomly permute by row
        randperm = np.random.permutation(n)
        for rownum in range(n):
            permmat[randperm[rownum],:] = distance_matrix[rownum,randperm]
        # Ensure symmetry
        permmat = (permmat + permmat.T) / 2
    elif flag == 'rand':
        # randomly permute all dist values (in upper section, to preserve symmetry)
        permidx = np.triu_indices(n,1)
        alldistvals = distance_matrix[permidx]
        permvals = np.random.permutation(alldistvals)
        permmat[permidx] = permvals
        permidx_lower = np.tril_indices(n,-1)
        permmat[permidx_lower] = permmat.T[permidx_lower]
    elif flag == 'dist':
        # matrix is purely distance dependent (same averages along non-main diagonals as original) + noise
        for diagnum in range(n-1):
            diagvals = np.diag(distance_matrix,diagnum+1)
            avgval = np.mean(diagvals)
            std = np.std(diagvals)
            newdiag = np.random.normal(avgval,std,len(diagvals))
            diagmat = np.diag(newdiag,diagnum+1)
            permmat += diagmat
            permmat += diagmat.T
    return permmat

def writeDistMatToFile(distmat,filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f,delimiter='\t')
        for row in distmat:
            fwriter.writerow(row)
    print('wrote distance matrix to',filename)
    return

def process_permutations(distance_matrix, output_path, output_name):
    # Generate and process random permutation
    print("Processing random permutation...")
    rand_perm = randomlyPermuteDistMat(distance_matrix, flag='rand')
    TDA_func(rand_perm, 
             os.path.join(output_path, output_name + '_randperm_persisdiagram.txt'))
    writeDistMatToFile(rand_perm, os.path.join(output_path, output_name + '_randperm_distmat.txt'))
    
    # Generate and process edge permutation
    print("Processing edge permutation...")
    edge_perm = randomlyPermuteDistMat(distance_matrix, flag='edge')
    TDA_func(edge_perm,
             os.path.join(output_path, output_name + '_edgeperm_persisdiagram.txt'))
    writeDistMatToFile(edge_perm, os.path.join(output_path, output_name + '_edgeperm_distmat.txt'))
    
    # Generate and process distance permutation
    print("Processing distance permutation...")
    dist_perm = randomlyPermuteDistMat(distance_matrix, flag='dist')
    TDA_func(dist_perm,
             os.path.join(output_path, output_name + '_distperm_persisdiagram.txt'))
    writeDistMatToFile(dist_perm, os.path.join(output_path, output_name + '_distperm_distmat.txt'))

def main(input_file,output_path,output_name,resol):
    t0 = time.time()
    resolut,freq_mat,min_bin,_=read_raw_HiC_data_no_split(input_file,resol)
    print("Generate distance matrix...")
    distance_matrix=matrix_normalize([freq_mat])[0]
    writeDistMatToFile(distance_matrix, os.path.join(output_path, output_name + '_distmat.txt'))
    
    print("Calculate persistent homology with Alpha Complex...")
    persisdiag,full_persis_info,mds_embedding,stress,min_lifespan,max_lifespan=TDA_func(distance_matrix, os.path.join(output_path, output_name + '_persisdiagram.txt'))
    
    print("Processing permutations...")
    process_permutations(distance_matrix, output_path, output_name)
    
    print('total time (s):',time.time()-t0)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str,help= "input file of HiC contact matrix")
    parser.add_argument('-o',type=str, help='the name of output files')
    parser.add_argument('-p',type=str,help="the path of output files")
    parser.add_argument('-r',type=str,help="resolution of HiC input file")

    args=parser.parse_args()
    main(args.i,args.p,args.o,args.r)