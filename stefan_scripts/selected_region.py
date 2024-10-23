import scanpy as sc
import numpy as np

ad = sc.read_h5ad('nico_celltype_annotation_full.h5ad')
print(ad)

coord=ad.obsm['spatial']

print(coord[0:4])

index=np.where( (4500<=coord[:,1])&(coord[:,1]<=5500) &(1500<=coord[:,0]) & (coord[:,0]<=2500))

#x = c(4500, 5500), y = c(1500, 2500))

print(index)

adata = ad[index[0],:].copy()
adata.write('nico_analysis/nico_celltype_annotation.h5ad')


ad = sc.read_h5ad('sct_spatial_full.h5ad')
adata = ad[index[0],:].copy()
adata.write('inputQuery/sct_spatial.h5ad')
