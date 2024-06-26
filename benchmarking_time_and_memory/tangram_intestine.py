import scanpy as sc 
import tangram as tg
import pandas as pd 
import numpy as np 


import time
from memory_profiler import memory_usage 


time_start= time.time()
#memory_start = tracemalloc.start()
#memory_start = psutil.Process().memory_info().rss / (1024 * 1024)
memory_before = memory_usage()[0]


scdatapath='./inputRef/'
spdatapath='./inputQuery/'

ad_sp = sc.read_h5ad(spdatapath+'spatial_intestine.h5ad')
ad_sc = sc.read_h5ad(scdatapath+'input_ref.h5ad')

print(ad_sc)

tg.pp_adatas(ad_sc, ad_sp, genes=None)

ad_map = tg.map_cells_to_space(ad_sc, ad_sp,mode='clusters',cluster_label='cluster')
#ad_map.write_h5ad("ad_map.h5ad")

ad_ge = tg.project_genes(ad_map, ad_sc,cluster_label='cluster')
#ad_ge.write_h5ad("Tangram_sc_genes_spatial_cell_cluster_level.h5ad")
tg.plot_training_scores(ad_map, bins=10, alpha=.5)




#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)

