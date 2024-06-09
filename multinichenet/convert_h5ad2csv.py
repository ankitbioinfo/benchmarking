import scanpy as sc
import pandas as pd
import numpy as np

ad=sc.read_h5ad('./inputQuery/spatial_aba.h5ad')
print(ad)
mat=ad.X.toarray()
df=pd.DataFrame(mat,index=ad.obs_names,columns=ad.var_names)
df.to_csv('countsCbyG.csv')

#df=pd.DataFrame(mat.T,columns=ad.obs_names,index=ad.var_names)
#df.to_csv('countsGbyC.csv')

vec=ad.obs['nico_ct']
vec = np.reshape(vec,(len(vec),1))

mat1=np.hstack((ad.obsm['spatial'],vec))
print(mat.shape,mat1.shape)

df=pd.DataFrame(mat1,index=ad.obs_names,columns=['x','y','annot'])
df.to_csv('meta.csv')

#print(ad.X)
