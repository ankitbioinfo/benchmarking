# Import necessary modules
import time
import anndata as ad
import squidpy as sq
import cellcharter as cc
import pandas as pd
import scanpy as sc
import scvi
import numpy as np
import matplotlib.pyplot as plt
import pickle
from memory_profiler import memory_usage

time_start= time.time()
memory_before = memory_usage()[0]  # Memory before execution


spdatapath='./inputQuery/'
adata=sc.read_h5ad(spdatapath+'common_counts_sp.h5ad')

#adata=sc.read(spdatapath+'gene_by_cell.csv').transpose()
coordinate = pd.read_csv(spdatapath+'tissue_positions_list.csv')
coordinate=coordinate.to_numpy()


sample=[]
for i in range(len(adata.obs_names)):
    sample.append('intestine')

print(coordinate.shape,adata)

df=pd.read_csv(spdatapath+'MNN_based_annotations/3_nico_annotation_ct_name.csv')
sc_ctype_name=df.to_numpy()

df=pd.read_csv(spdatapath+'MNN_based_annotations/3_nico_annotation_cluster.csv')
nico_cluster=df.to_numpy()
cellname=nico_cluster[:,0]

print("equal",np.array_equal(cellname,adata.obs_names))


d={}
for i in range(len(sc_ctype_name)):
    d[i]=sc_ctype_name[i][1]
ctname=[]
for i in range(len(cellname)):
        ctname.append( d[ nico_cluster[i,1]])

adata.obs['nico_ct']=ctname

sc_ctype_name=coordinate[:,0]
d={}
for i in range(len(sc_ctype_name)):
    d[sc_ctype_name[i]]=i
index=[]
for i in range(len(cellname)):
    index.append( d[ cellname[i]])
coordinate1=coordinate[index]
adata.obsm['spatial']=coordinate1[:,[1,2]].astype(float) 
print("equal",np.array_equal(cellname,coordinate1[:,0]),coordinate1.shape)

'''
ad_spatial_ori=sc.read_h5ad('small_slide1120.h5ad')
ad_spatial_ori.obs['nico_ct']=ad_spatial_ori.obs['cluster']

adata=ad_spatial_ori
'''


adata.obs['sample']=np.array(sample)
adata.obs['sample']=pd.Categorical(adata.obs['sample'])

print('1',adata)

index=[]
for i in range(len(ctname)):
    flag=1
    if adata.obs['nico_ct'][i]=='NM':
        flag=0
    if adata.obs['nico_ct'][i]=='HsPCs':
        flag=0
    if adata.obs['nico_ct'][i]=='Basophils':
        flag=0
    if adata.obs['nico_ct'][i]=='NK cells':
        flag=0
    if flag==1:
        index.append(i)


adata=adata[index]
print('2',adata)



adata.raw = adata.copy()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

scvi.settings.seed = 12345
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key='sample'
)
model = scvi.model.SCVI(adata)

model.train(early_stopping=True, enable_progress_bar=True)
adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
sq.gr.spatial_neighbors(adata, library_key='sample', coord_type='generic', delaunay=True)
cc.gr.remove_long_links(adata)

cc.gr.aggregate_neighbors(adata, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter',
                          sample_key='sample')

autok = cc.tl.ClusterAutoK(
    n_clusters=(2,25),
    max_runs=50,
    model_params=dict(
        random_state=14533
        # If running on GPU
        #trainer_params=dict(accelerator='gpu', devices=1)
    )
)
autok.fit(adata, use_rep='X_cellcharter')
cc.pl.autok_stability(autok,save='stabilitynew.pdf')
adata.obs['cluster_cellcharter'] = autok.predict(adata, use_rep='X_cellcharter')





cc.gr.enrichment(adata, group_key='nico_ct', label_key='cluster_cellcharter')
cc.pl.enrichment(adata, group_key='nico_ct', label_key='cluster_cellcharter', color_threshold=0.58,
                 size_threshold=0.58,save="new_niche.pdf")





fout='save_data_withoutNM.p'
myfile = open(fout, 'wb')
pickle.dump(adata,myfile)
myfile.close()




#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)

