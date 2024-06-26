import uniport as up
import scanpy as sc
import pandas as pd
import numpy as np 
import warnings 
np.warnings = warnings 


import time
from memory_profiler import memory_usage 


time_start= time.time()
#memory_start = tracemalloc.start()
#memory_start = psutil.Process().memory_info().rss / (1024 * 1024)
memory_before = memory_usage()[0]


scdatapath='./inputRef/'
spdatapath='./inputQuery/'

spa_data = sc.read_h5ad(spdatapath+'spatial_intestine.h5ad')
seq_data = sc.read_h5ad(scdatapath+'input_ref.h5ad')


sc.pp.filter_cells(spa_data, min_counts=5)
sc.pp.filter_cells(seq_data, min_counts=5)
print(spa_data)
print(seq_data)


seq_data.obs['cell_type']=seq_data.obs['cluster']
seq_data.obs['domain_id']=0
seq_data.obs['domain_id']=seq_data.obs['domain_id'].astype('category')
seq_data.obs['source']='RNA'
#seq_data.obs['labels']=sc_cluster[:,1]


spa_data.obs['cell_type']=spa_data.obs['nico_ct']
spa_data.obs['domain_id']=1
spa_data.obs['domain_id']=spa_data.obs['domain_id'].astype('category')
spa_data.obs['source']='MERFISH'

adata_cm=spa_data.concatenate(seq_data,join='inner',batch_key='domain_id')
print(adata_cm)

sc.pp.normalize_total(adata_cm)
sc.pp.log1p(adata_cm)
sc.pp.highly_variable_genes(adata_cm, n_top_genes=2000,batch_key='domain_id',inplace=False,subset=True)
up.batch_scale(adata_cm)


sc.pp.normalize_total(seq_data)
sc.pp.log1p(seq_data)
sc.pp.highly_variable_genes(seq_data, n_top_genes=2000,batch_key='domain_id',inplace=False,subset=True)
up.batch_scale(seq_data)

sc.pp.normalize_total(spa_data)
sc.pp.log1p(spa_data)
sc.pp.highly_variable_genes(spa_data, n_top_genes=2000,batch_key='domain_id',inplace=False,subset=True)
up.batch_scale(spa_data)


#spa_data.write('MERFISH/MERFISH_processed.h5ad', compression='gzip')
#seq_data.write('MERFISH/RNA_processed.h5ad', compression='gzip')
#adata_cm.write('MERFISH/MERFISH_and_RNA.h5ad', compression='gzip')

adata_cm_copy = adata_cm.copy()
sc.pp.pca(adata_cm_copy)
sc.pp.neighbors(adata_cm_copy)
sc.tl.umap(adata_cm_copy, min_dist=0.1)
#sc.pl.umap(adata_cm_copy, color=['source', 'cell_type'])

adata = up.Run(adatas=[spa_data, seq_data], adata_cm=adata_cm)

RNA_ref=adata[adata.obs['source']=='RNA']
MERFISH_query=adata[adata.obs['source']=='MERFISH']
MERFISH_predict_label=up.metrics.label_transfer(ref=RNA_ref,query=MERFISH_query,label='cell_type')


'''
df=pd.DataFrame(MERFISH_predict_label,index=sp_cellname)
df.to_csv('annotationsPrediction_uniPort_Embryo_CTname.csv',index=True) #,header=['barcode','uniPortAnnotation']

ctname=np.unique(MERFISH_predict_label)
d={}
fw=open('uniPort_celltypeName.dat','w')
for i in range(len(ctname)):
    fw.write(str(i)+','+ctname[i]+'\n')
    d[ctname[i]]=i
fw.close()

fw=open('uniPort_CLID.csv','w')
for i in range(len(sp_cellname)):
    clid=d[MERFISH_predict_label[i]]
    fw.write(sp_cellname[i]+','+str(clid)+'\n')
fw.close()
'''    


#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)


