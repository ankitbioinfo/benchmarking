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
vec=adata.obs['nico_ct']
sample = []

'''
for i in range(len(vec)):
	if vec[i]=='cDC/monocyte':
		vec[i]='cDCmonocyte'
	if vec[i]=='Stem/TA':
		vec[i]='StemTA'
	if vec[i]=='neurons/enteroendocrine':
		vec[i]='neuronsEnteroendocrine'
	if vec[i]=='Rest B':
		vec[i]='RestB'
	if vec[i]=='Blood vasc.':
		vec[i]='BloodVasc'
'''


for i in range(len(vec)):
    r=np.random.rand(1)
    if r[0]>0.5:
        sample.append('g1')
    else:
        sample.append('g2')

sample=np.array(sample)
print(sample[0:5])

vec = np.reshape(vec,(len(vec),1))
sample = np.reshape(sample,(len(sample),1))
print("samp",sample.shape)
mat1=np.hstack((adata.obsm['spatial'],vec,sample))


df=pd.DataFrame(mat1,index=adata.obs_names,columns=['x','y','annot','sample'])
df.to_csv(outdir+'meta.csv')


