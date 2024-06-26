#from nico import Annotations as sann
#from nico import Interactions as sint
#from nico import Covariations as scov

import Annotations as sann
import Interactions as sint
import Covariations as sppath
import scanpy as sc
import numpy as np

#import scanpy as sc
#import gseapy
#import xlsxwriter

#import numpy as np
import time
from memory_profiler import memory_usage 
#import os

#font_dir = '/Library/Fonts'

import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']
plt.rcParams['pdf.fonttype'] = 42  # Embed fonts in PDF files
#plt.rcParams['font.serif'] = [font_dir + 'Helvetica.ttc']
#plt.rcParams['font.sans-serif'] = [font_dir + 'Helvetica.ttc']

#import warnings
#warnings.filterwarnings("ignore")
#plt.rc('font', family='sans-serif')
#plt.rc('font','DejaVu Sans')

#https://github.com/satijalab/sctransform/issues/45




time_start= time.time()
#memory_start = tracemalloc.start()
#memory_start = psutil.Process().memory_info().rss / (1024 * 1024)
memory_before = memory_usage()[0]





ref_datapath='./inputRef/'
query_datapath='./inputQuery/'



#output_annotation_dir='./MNN_based_annotations/'
output_annotation_dir=None
#output_nico_dir='./spatial_ct_ct_interactions/'
output_nico_dir=None


anchors_and_neighbors_info=sann.find_anchor_cells_between_ref_and_query(refpath=ref_datapath,quepath=query_datapath,
output_annotation_dir=output_annotation_dir)

output_info=sann.nico_based_annotation(anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.5',
across_spatial_clusters_dispersion_cutoff=0.15,
resolved_tie_issue_with_weighted_nearest_neighbor='No')

'''
sann.visualize_umap_and_cell_coordinates_with_all_celltypes(quepath=query_datapath,
saveas='pdf',output_annotation_dir=output_annotation_dir)

# For visualizing every cell type individually, leave list choose_celltypes empty.

sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(quepath=query_datapath,
choose_celltypes=[],saveas='pdf',output_annotation_dir=output_annotation_dir)

sann.save_annotations_in_spatial_object(output_info,output_h5ad_name='spatial_intestine_with_nico_annotations')
'''


#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)


