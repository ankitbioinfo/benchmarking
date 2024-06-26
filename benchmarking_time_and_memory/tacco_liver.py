

import os
import sys
import json
import matplotlib

import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import seaborn as sns
import sklearn.metrics
import tacco as tc
import time 
import tracemalloc 
from memory_profiler import memory_usage 

time_start= time.time()
#memory_start = tracemalloc.start()
#memory_start = psutil.Process().memory_info().rss / (1024 * 1024)
memory_before = memory_usage()[0]


# The notebook expects to be executed either in the workflow directory or in the repository root folder...
#sys.path.insert(1, os.path.abspath('workflow' if os.path.exists('workflow/common_code.py') else '..'))
#import common_code

path1='./28.12.2023/inputQuery/'
path2='./28.12.2023/inputRef/'


adseq = sc.read_h5ad(path2+'input_ref.h5ad')
adspa = sc.read_h5ad(path1+'common_counts_sp.h5ad')

#adseq.X.toarray()

print(adseq)
print(adspa)


# Filter the cells and genes 
sc.pp.filter_cells(adspa, min_counts=5)
sc.pp.filter_cells(adseq, min_counts=5)

#sc.pp.filter_genes(adspa, min_cells=1)
sc.pp.filter_genes(adseq, min_cells=3)

print(adspa)
print(adseq)

np.unique(adseq.obs['cluster'])

adseq.obs['cell'] = adseq.obs.index


tc.tl.annotate(adspa,adseq,'cluster',result_key='Tacco_singleCell_OT',multi_center=None,method='OT')
#tc.tl.annotate(adspa,adseq,'cluster',result_key='Tacco_singleCell_WOT',multi_center=None,method='WOT')
tc.tl.annotate(adspa,adseq,'cluster',result_key='Tacco_singleCell_withoutmulticenter')


adspa.write_h5ad('tacco_output.h5ad')



#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)










