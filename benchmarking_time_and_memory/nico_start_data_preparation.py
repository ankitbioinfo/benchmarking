import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
#import pickle
#from scipy.spatial import cKDTree
#from SCTransform import SCTransform


#from nico import Annotations as sann
import Annotations as sann

from memory_profiler import memory_usage 
import warnings
import time
warnings.filterwarnings('ignore')
#export PYTHONWARNINGS='ignore:Multiprocessing-backed parallel loops:UserWarning'
os.environ["PYTHONWARNINGS"] = "ignore::UserWarning"



time_start= time.time()
#memory_start = tracemalloc.start()
#memory_start = psutil.Process().memory_info().rss / (1024 * 1024)
memory_before = memory_usage()[0]



# This function find the common gene index between two data.  
def find_index(sp_genename,sc_genename):
    index_sc=[]
    index_sp=[]
    d={}
    for j in range(len(sc_genename)):
        name=sc_genename[j]
        d[name]=j

    for i in range(len(sp_genename)):
        name=sp_genename[i]
        try:
            d[name]
            flag=1
        except KeyError:
            flag=0
        if flag==1:
            index_sc.append(d[name])
            index_sp.append(i)
    return index_sp,index_sc


# This is input data path for the scRNA-seq and spatial data 
# If the data is not in h5ad or csv format then please adjust following the standard scanpy routine. 

scdatapath='./inputRef/'
spdatapath='./inputQuery/'

ad_spatial_ori=sc.read(spdatapath+'gene_by_cell.csv').transpose()
ad_seq_ori=sc.read_h5ad(scdatapath+'input_ref.h5ad') 

# This is the coordinate file of the cell centroids from the spatial transcriptomics experiment.  
coordinate = pd.read_csv(spdatapath+'tissue_positions_list.csv')
coordinate=coordinate.to_numpy()
ad_spatial_ori.obsm['spatial']=coordinate[:,1:].astype(float)

print('check equality',np.array_equal(coordinate[:,0],ad_spatial_ori.obs_names))


# data size of the variables 
print(ad_spatial_ori)
print(ad_seq_ori)
print(coordinate.shape)


# Filter the cells and genes 
sc.pp.filter_cells(ad_spatial_ori, min_counts=5)
sc.pp.filter_cells(ad_seq_ori, min_counts=5)

sc.pp.filter_genes(ad_spatial_ori, min_cells=1)
sc.pp.filter_genes(ad_seq_ori, min_cells=1)

print(ad_spatial_ori)
print(ad_seq_ori)


Original_counts=ad_seq_ori.copy()
Original_counts.raw=Original_counts.copy()

# Standard scanpy analysis 

sc.pp.normalize_total(Original_counts)
sc.pp.log1p(Original_counts)

sc.tl.pca(Original_counts)
sc.pp.neighbors(Original_counts)
sc.tl.umap(Original_counts)
#sc.pl.umap(Original_counts)


Original_counts.write_h5ad(scdatapath+'Original_counts.h5ad')

sp_genename=ad_spatial_ori.var_names.to_numpy()
sc_genename=ad_seq_ori.var_names.to_numpy()

index_sp,index_sc=find_index(sp_genename,sc_genename)
#ad_seq_common=ad_seq_ori[:,index_sc].copy()
ad_seq_common=ad_seq_ori.copy()
ad_spatial_common=ad_spatial_ori[:,index_sp].copy()

print(len(sp_genename[index_sp]))
print(len(sc_genename[index_sc]))


# Alternative 1 
# The sctransform normalization function used from scanpy 


ad_seq_common.raw=ad_seq_common.copy()
ad_spatial_common.raw=ad_spatial_common.copy()
# perform scTranform normalization common gene space for spatial data and scRNAseq data  
sc.experimental.pp.normalize_pearson_residuals(ad_seq_common,inplace=True) #ad_seq_common.X[ad_seq_common.X<0]=0
ad_seq_common=ad_seq_common[:,index_sc]

ad_seq_common.write_h5ad(scdatapath+'sct_singleCell.h5ad')
sc.experimental.pp.normalize_pearson_residuals(ad_spatial_common,inplace=True) #ad_spatial_common.X[ad_spatial_common.X<0]=0
#print(ad_spatial_common.X.toarray()


# standard scanpy analysis 
sc.pp.pca(ad_spatial_common)
sc.pp.neighbors(ad_spatial_common,n_pcs=30)
sc.tl.umap(ad_spatial_common)

sc.tl.leiden(ad_spatial_common, resolution=0.4,key_added="leiden0.4")
sc.tl.leiden(ad_spatial_common, resolution=0.5,key_added="leiden0.5")

ad_spatial_common.write_h5ad(spdatapath+'sct_spatial.h5ad')



#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)

