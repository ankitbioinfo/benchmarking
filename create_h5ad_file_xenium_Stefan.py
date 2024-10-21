

import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_10x_h5('cell_feature_matrix.h5')
print(adata)
cellname=adata.obs_names

df=pd.read_csv('cells.csv')
print(df.columns)

data = df[['x_centroid', 'y_centroid']]
cell_id= df[['cell_id']].to_numpy()
print(cellname.shape,cell_id.shape)

print('equal',np.array_equal(cell_id[:,0],cellname))
print(cellname[0:5],cellname[-1])
print(cell_id[0:5],cell_id[-1])

adata.obsm['spatial']=data.to_numpy()

adata.write_h5ad('query.h5ad')
