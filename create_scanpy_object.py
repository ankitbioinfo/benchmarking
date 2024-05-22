

import scanpy as sc
import pandas as pd
import numpy as np

'''
https://github.com/dmcable/spacexr/tree/master/documentation

Doublet mode results: The results of ‘doublet_mode’ are stored in @results$results_df.
More specifically, the results_df object contains one column per pixel (barcodes as rownames).
Important columns are: * spot_class, a factor variable representing RCTD’s classification in doublet mode:
“singlet” (1 cell type on pixel),
“doublet_certain” (2 cell types on pixel),
“doublet_uncertain” (2 cell types on pixel, but only confident of 1),
“reject” (no prediction given for pixel). Typically, reject pixels should not be used, and
doublet_uncertain pixels can only be used for applications that do not require knowledge of all cell types on a pixel.
* Next, the first_type column gives the first cell type predicted on the bead (for all spot_class conditions except “reject”).
* The second_type column gives the second cell type predicted on the bead for doublet spot_class conditions
 (not a confident prediction for “doublet_uncertain”).
'''

df1=pd.read_csv('gname.csv',index_col=0)
d1=df1.to_numpy()

df2=pd.read_csv('cname.csv',index_col=0)
d2=df2.to_numpy()

df3=pd.read_csv('coords.csv')
coord=df3.to_numpy()

df4= pd.read_csv('rctd_doublet.csv')
annot=df4.to_numpy()


print(d1.shape,d2.shape,coord.shape,annot.shape)


ad = sc.read_mtx('matrix.mtx').transpose()
ad.obs_names = d2[:,0]
ad.var_names = d1[:,0]
print(ad, np.array_equal(d2[:,0],coord[:,0]), np.array_equal(d2[:,0],annot[:,0]),       )

ad.obsm['spatial']=coord[:,1:].astype(float)
ad.obs['rctd_class']= annot[:,1]
ad.obs['rctd_first_type']= annot[:,2]

#index= np.where(ad.obs['rctd_class']!='reject')

index= np.where((ad.obs['rctd_class']=='singlet')|(ad.obs['rctd_class']=='doublet_uncertain'))


adata = ad[index].copy()

a=set(ad.obs['rctd_class'])
b=set(adata.obs['rctd_class'])
print(a)
print(b)

print(adata)

ct = set(adata.obs['rctd_first_type'])
print(ct)

adata.write_h5ad('cerebellum.h5ad')
