
import scanpy as sc
import pandas as pd
import numpy as np

df1=pd.read_csv('gene_names.csv',) #index_col=0
d1=df1.to_numpy()

df2=pd.read_csv('cell_metadata.csv') #,index_col=0
d2=df2.to_numpy()

#df3=pd.read_csv('coords.csv')
#coord=df3.to_numpy()

df4= pd.read_csv('umap_coordinates.csv')
umap=df4.to_numpy()

annot = d2[:,1]


#0  NK
#1  T
#2  Cd14+
#3  FCGR3A+
#4  T
#5  (FCER1A) DC
#6  BEC
#7  KC (VSIG4 und wenig CD5l CLEC4f aber auch TREM2)
#8  LSEC
#9  B
#10 FCGR3A+
#11 B
#12 CEC
#13 cycling myeloid (MKI67+)
#14 BEC
#15 HEP
#16 unclear - doublet
#17 FIBS (PDGRFB+)
#18 DC (SCT PTCRA pDC daher auch als DC deklariert)
#19 Doublet

celltype=["NK", "T", "CD14+", "FCGR3A+", "T", "DC_FCER1A", "BEC", "KC", "LSEC", "B", "FCGR3A+", "B", "CEC", "Cycling Myeloid", "BEC", "HEP","NM", "FIBS", "DC","Doublet"]

print("unique",np.unique(annot))
for i in range(len(annot)):
    #print(annot[i])
    annot[i]=celltype[annot[i]]


print(d1.shape,d2.shape,umap.shape)
print(d1[0:3])
print(d2[0:3])

ad = sc.read_mtx('counts_matrix.mtx').transpose()
ad.obs_names = d2[:,0]
ad.var_names = d1[:,0]

ad.obs['cell_type']=annot
ad.obsm['X_umap'] = umap[:,[1,2]].astype(float)

index=[]
annot =ad.obs['cell_type']
for i in range(len(annot)):
    flag=1
    if annot[i]=="NM":
        flag=0
    if annot[i]=="Doublet":
        flag=0
    if flag==1:
        index.append(i)

print(ad)
adata=ad[index,:].copy()
print(adata)
ad.write_h5ad('reference.h5ad')
