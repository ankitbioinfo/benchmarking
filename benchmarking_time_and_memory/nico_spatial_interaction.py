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


time_start= time.time()
#memory_start = tracemalloc.start()
#memory_start = psutil.Process().memory_info().rss / (1024 * 1024)
memory_before = memory_usage()[0]



ref_datapath='./inputRef/'
query_datapath='./inputQuery/'
output_annotation_dir=None
output_nico_dir=None
ct_not_good=['Basophils','HsPCs']



#clusterFilename=output_annotation_folder+'3_nico_annotation_cluster.csv'
#celltypeFilename=output_annotation_folder+'3_nico_annotation_ct_name.csv'
clusterFilename=None
celltypeFilename=None


niche_pred_output=sint.spatial_neighborhood_analysis(Radius=0,clusterFilename=clusterFilename,
quepath=query_datapath,
celltypeFilename=celltypeFilename,
output_annotation_dir=output_annotation_dir,
output_niche_prediction_dir=output_nico_dir,removed_CTs_before_finding_CT_CT_interactions=ct_not_good)


sint.find_interacting_cell_types(niche_pred_output,choose_celltypes=[],
coeff_cutoff=30,saveas='png',transparent_mode='False',figsize=(4.0,2.0))
#sint.plot_roc_results(niche_pred_output)
sint.plot_confusion_matrix(niche_pred_output,saveas='png')
sint.plot_coefficient_matrix(niche_pred_output,saveas='png')
#st.plot_predicted_probabilities(niche_pred_output)

#Find the niche interactions for other cutoff
#sint.plot_niche_interactions_with_edge_weight(niche_pred_output,saveas='png',niche_cutoff=0.18)
#sint.plot_niche_interactions_without_edge_weight(niche_pred_output,saveas='png',niche_cutoff=0.18)

sint.plot_evaluation_scores(niche_pred_output,saveas='png', figsize=(4,3))


#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)


