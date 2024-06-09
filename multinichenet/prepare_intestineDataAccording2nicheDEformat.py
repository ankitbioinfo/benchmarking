import time
import anndata as ad
#import squidpy as sq
#import cellcharter as cc
import pandas as pd
import scanpy as sc
#import scvi
import numpy as np
import matplotlib.pyplot as plt
#from memory_profiler import memory_usage
from sklearn.preprocessing import OneHotEncoder


def get_hot_array(spatial_data,cell_type_key):
    one_hot_enc = OneHotEncoder().fit(np.asarray(list(set(spatial_data.obs[cell_type_key]))).reshape([-1,1]))
    cell_type_one_hot = one_hot_enc.transform(np.asarray(spatial_data.obs[cell_type_key]).reshape([-1,1])).reshape([spatial_data.obs[cell_type_key].shape[0], -1]).todense()
    myct=one_hot_enc.categories_[0]
    return cell_type_one_hot,myct



def remove_extra_character_from_name(name):
    """
    This function removes special characters from the cell type names to avoid throwing an error while saving the figures.
    """
    name=name.replace('/','_')
    name=name.replace(' ','_')
    name=name.replace('"','')
    name=name.replace("'",'')
    name=name.replace(')','')
    name=name.replace('(','')
    name=name.replace('+','p')
    name=name.replace('-','n')
    name=name.replace('.','')
    return name

def format_ct_in_good_string(array):
    value=[]
    for i in range(len(array)):
        value.append(remove_extra_character_from_name(array[i]))

    return value


spdatapath='./inputQuery/'
ad_spatial_ori=sc.read(spdatapath+'gene_by_cell.csv').transpose()
coordinate = pd.read_csv(spdatapath+'tissue_positions_list.csv')
coordinate=coordinate.to_numpy()

#because it is 2d tissue slide so only upload the X and Y coordiante
ad_spatial_ori.obsm['spatial']=coordinate[:,[1,2]].astype(float)
adata=ad_spatial_ori
sample=[]
for i in range(len(ad_spatial_ori.obs_names)):
    sample.append('intestine')


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
adata.obs['sample']=np.array(sample)
adata.obs['sample']=pd.Categorical(adata.obs['sample'])

print('1',adata)

index=[]
for i in range(len(ctname)):
    flag=1
    if adata.obs['nico_ct'][i]=='NM':
        flag=0
    if adata.obs['nico_ct'][i]=='Cycling/GC B cell':
        flag=0
    if adata.obs['nico_ct'][i]=='pDC':
        flag=0
    if flag==1:
        index.append(i)

adata=adata[index]
print('2',adata)

outdir='./intestine_for_nicheDE/'

mat=adata.X
print(mat)
print(mat.shape)
df=pd.DataFrame(mat,index=adata.obs_names,columns=adata.var_names)
df.to_csv(outdir+'counts.csv')

df=pd.DataFrame(adata.obsm['spatial'],index=adata.obs_names,columns=['imagerow',	'imagecol'])
df.to_csv(outdir+'coord.csv')

#unique_spct= sorted(list(set(adata.obs['nico_ct'])))

one_hot_encoding,unique_spct=get_hot_array(adata,'nico_ct')
print(one_hot_encoding.shape)
print(one_hot_encoding[0:3])

unique_spct_good=format_ct_in_good_string(unique_spct)
df=pd.DataFrame(one_hot_encoding,index=adata.obs_names,columns=unique_spct_good)
df.to_csv(outdir+'deconv.csv')



scdata=sc.read_h5ad('inputRef/Original_counts.h5ad')
mat=scdata.raw.X.todense()
cluster=scdata.obs['cluster']
print(mat,mat.shape,cluster.shape)
unique_scct= sorted(list(set(cluster)))

n=len(unique_scct)
m=len(unique_spct)
ngene=scdata.var_names

avg_scExp=np.zeros((m,len(ngene)))
for i in range(m):
    for j in range(n):
        if unique_scct[j]==unique_spct[i]:
            index=np.where(unique_scct[j]==cluster)
            a=mat[index[0]]
            b=np.mean(a,axis=0)
            print(a.shape,b.shape)
            avg_scExp[i]=b
    #print(ct[i],a.shape,b.shape,len(index[0]))
    #print(ct[i],b[0:10])

print(avg_scExp.shape,len(ngene),len(unique_spct))
df=pd.DataFrame(avg_scExp,index=unique_spct_good,columns=ngene)
df.to_csv(outdir+'library_avgExp.csv')




'''
cd <- read.csv('author_data/counts.csv', row.names = 1)
pos <- read.csv('author_data/coord.csv', row.names = 1)
celltype <- read.csv('author_data/deconv.csv', row.names = 1)
library <- read.csv('author_data/library.csv',row.names=1)
'''
