import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import time 

from memory_profiler import memory_usage

time_start= time.time()
memory_before = memory_usage()[0]  # Memory before execution


def create_directory(outputFolder):
    "This function creates empty directory."
    answer=os.path.isdir(outputFolder)
    if answer==True:
        pass
    else:
        os.mkdir(outputFolder)


adata = sc.read_h5ad("./inputQuery/spatial_intestine.h5ad")
coord= adata.obsm['spatial']
#adata.obs["x_array"]=adata.obs["x2"]
#adata.obs["y_array"]=adata.obs["x3"]
adata.obs["x_pixel"]=coord[:,0]
adata.obs["y_pixel"]=coord[:,1]
#Select captured samples
#adata=adata[adata.obs["x1"]==1]
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")

'''
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

#maindata=adata.copy()
'''

x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["y_pixel"].tolist()

for i in range(len(x_pixel)):
    x=x_pixel[i]
    y=y_pixel[i]


#Calculate adjacent matrix
s=1
b=49
#adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, beta=b, alpha=s, histology=True) # image=img,
#If histlogy image is not available, SpaGCN can calculate the adjacent matrix using the fnction below
adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
np.savetxt('adj.csv', adj, delimiter=',')


#adata=sc.read("../data/151673/sample_data.h5ad")
adj=np.loadtxt('adj.csv', delimiter=',')
adata.var_names_make_unique()
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)
#Normalize and take log for UMI
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)


#Set hyper-parameters
#p: Percentage of total expression contributed by neighborhoods.
#l: Parameter to control p.

p=0.5
#Find the l value given p
l=spg.search_l(p, adj, start=0.001, end=1000, tol=0.01, max_run=100)


r_seed=t_seed=n_seed=100
'''
#If the number of clusters known, we can use the spg.search_res() fnction to search for suitable resolution(optional)
#For this toy data, we set the number of clusters=7 since this tissue has 7 layers
n_clusters=7
#Set seed

#Seaech for suitable resolution
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05,
 max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
'''

res=1.0
seed=100

clf=spg.SpaGCN()
clf.set_l(l)
#Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
#Run
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
y_pred, prob=clf.predict()
adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')
#Do cluster refinement(optional)
#shape="hexagon" for Visium data, "square" for ST data.

'''
adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
adata.obs["refined_pred"]=refined_pred
adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
#Save results
'''
adata.write_h5ad("spagcn_output.h5ad")




#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)

