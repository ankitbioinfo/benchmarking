import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os, time 
import sys
from memory_profiler import memory_usage
from scipy import sparse

import torch
import STAGATE_pyG
#import STAGATE


from sklearn.metrics.cluster import adjusted_rand_score


time_start= time.time()
#memory_start = tracemalloc.start()
#memory_start = psutil.Process().memory_info().rss / (1024 * 1024)
memory_before = memory_usage()[0]

section_id = '151676'

#input_dir = os.path.join('Data', section_id)
#adata = sc.read_visium(path=input_dir, count_file=section_id+'_filtered_feature_bc_matrix.h5')
#adata.var_names_make_unique()
adata=sc.read_h5ad('inputQuery/spatial_intestine.h5ad')

print(adata.X)

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

#Normalization
#sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# read the annotation
#Ann_df = pd.read_csv(os.path.join('Data', section_id, section_id+'_truth.txt'), sep='\t', header=None, index_col=0)
#Ann_df.columns = ['Ground Truth']

#adata.obs['Ground Truth'] = Ann_df.loc[adata.obs_names, 'Ground Truth']
adata.obs['Ground Truth'] = adata.obs['nico_ct']

#adata.X = sparse.csr_matrix(adata.X)


STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=83)
STAGATE_pyG.Stats_Spatial_Net(adata)

adata = STAGATE_pyG.train_STAGATE(adata, device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu'))


print(adata)

sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)

#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)


sc.tl.louvain(adata,resolution=0.5,key_added="louvain0.5")
sc.tl.louvain(adata,resolution=0.8,key_added="louvain0.8")
sc.tl.louvain(adata,resolution=1,key_added="louvain1")
sc.tl.louvain(adata,resolution=0.2,key_added="louvain0.2")

adata.write_h5ad('stagate_intestine.h5ad')


#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)



