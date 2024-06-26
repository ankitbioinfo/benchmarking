import warnings
warnings.filterwarnings("ignore") 
import os, time
import pandas as pd
import numpy as np
import scanpy as sc 
import random
import scipy.sparse as sparse
from scipy.sparse import csr_matrix, issparse

from banksy.initialize_banksy import initialize_banksy
from banksy.run_banksy import run_banksy_multiparam
from banksy_utils.color_lists import spagcn_color
from memory_profiler import memory_usage

from banksy.main import median_dist_to_nearest_neighbour
from banksy.embed_banksy import generate_banksy_matrix
from banksy.main import concatenate_all
from banksy_utils.umap_pca import pca_umap
from banksy.cluster_methods import run_Leiden_partition
from banksy.plot_banksy import plot_results

start = time.perf_counter_ns()
random_seed = 1234
cluster_algorithm = 'leiden'
np.random.seed(random_seed)
random.seed(random_seed)

from memory_profiler import memory_usage

# put banksy dir from banksy_py into same path as this script 

time_start= time.time()
memory_before = memory_usage()[0]  # Memory before execution


def create_directory(outputFolder):
    "This function creates empty directory."
    answer=os.path.isdir(outputFolder)
    if answer==True:
        pass
    else:
        os.mkdir(outputFolder)


from banksy_utils.load_data import load_adata
file_path = os.path.join("inputQuery", "./")
adata_filename = "spatial_intestine.h5ad"

gcm_filename = ""
locations_filename = ""
load_adata_directly = True
save_fig = True
#output_folder = os.path.join(os.getcwd(), 'data', 'starmap', 'tmp_png', f'{cluster_algorithm}', f'seed{random_seed}')
output_folder = './tmp_intestine_banksy/'
create_directory(output_folder)
# Colour map
c_map = 'tab20'

coord_keys = ('x', 'y', 'spatial')
num_clusters = 20
resolutions = [.9] # clustering resolution for Leiden clustering
pca_dims = [20] # number of dimensions to keep after PCA
lambda_list = [.8] # lambda
k_geom = 15 # 15 spatial neighbours
max_m = 1 # use AGF
nbr_weight_decay = "scaled_gaussian" # can also be "reciprocal", "uniform" or "ranked"

sample = 'starmap'


raw_y, raw_x, adata = load_adata(file_path,
                                 load_adata_directly,
                                 adata_filename,
                                 gcm_filename,
                                 locations_filename,
                                 coord_keys)
#adata.var_names_make_unique()


print('1',adata)
index=[]
for i in range(len(adata.obs_names)):
    flag=1
    if adata.obs['nico_ct'][i]=='NM':
        flag=0
    if adata.obs['nico_ct'][i]=='Cycling/GC B cell':
        flag=1
    if adata.obs['nico_ct'][i]=='pDC':
        flag=1
    if flag==1:
        index.append(i)

adata=adata[index]
print('2',adata)



banksy_dict = initialize_banksy(
    adata,
    coord_keys,
    k_geom,
    nbr_weight_decay=nbr_weight_decay,
    max_m=max_m,
    plt_edge_hist=True,
    plt_nbr_weights=True,
    plt_agf_angles=True,
    plt_theta=True,
)

banksy_dict, banksy_matrix = generate_banksy_matrix(adata,
                                                    banksy_dict,
                                                    lambda_list,
                                                    max_m)


banksy_dict["nonspatial"] = {
    # Here we simply append the nonspatial matrix (adata.X) to obtain the nonspatial clustering results
    0.0: {"adata": concatenate_all([adata.X], 0, adata=adata), }
}

print(banksy_dict['nonspatial'][0.0]['adata'])


pca_umap(banksy_dict,
         pca_dims = pca_dims,
         add_umap = True,
         plt_remaining_var = False,
         )


seed=123456
results_df, max_num_labels = run_Leiden_partition(
    banksy_dict,
    resolutions,
    num_nn = 50,
    num_iterations = -1,
    partition_seed = seed,
    match_labels = True,
)


c_map =  'tab20' # specify color map
weights_graph =  banksy_dict['scaled_gaussian']['weights'][0]

plot_results(
    results_df,
    weights_graph,
    c_map,
    match_labels = True,
    coord_keys = coord_keys,
    max_num_labels  =  max_num_labels, 
    save_path = os.path.join(output_folder, 'tmp_png'),
    save_fig = False, # save the spatial map of all clusters
    save_seperate_fig = False, # save the figure of all clusters plotted seperately
)


'''
A=banksy_dict['nonspatial'][0.0]['adata'].obs['labels_nonspatial_pc20_nc0.00_r0.90']
#B=banksy_dict['scaled_gaussian'][0.2]['adata'].obs['labels_scaled_gaussian_pc20_nc0.20_r0.90']
C=banksy_dict['scaled_gaussian'][0.8]['adata'].obs['labels_scaled_gaussian_pc20_nc0.80_r0.90']
A.to_csv(output_folder+'nonspatial_L0.csv')
#B.to_csv(output_folder+'SG_L02.csv')
C.to_csv(output_folder+'SG_L08.csv')
'''





#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)

