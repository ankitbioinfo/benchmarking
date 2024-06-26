import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
import pandas as pd 
from cell2location.utils.filtering import filter_genes

from cell2location.models import RegressionModel

import time
from memory_profiler import memory_usage 


time_start= time.time()
#memory_start = tracemalloc.start()
#memory_start = psutil.Process().memory_info().rss / (1024 * 1024)
memory_before = memory_usage()[0]



results_folder = './c2l_results/c2l/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'


scdatapath='./inputRef/'
spdatapath='./inputQuery/'

spa_data = sc.read_h5ad(spdatapath+'spatial_intestine.h5ad')
seq_data = sc.read_h5ad(scdatapath+'input_ref.h5ad')

print("sp",spa_data.X)
print("seq",seq_data.X)


#default parameter 
cell_count_cutoff=5, 
cell_percentage_cutoff2=0.03, 
nonz_mean_cutoff=1.12


spa_data.obs['domain_id']=0
spa_data.obs['domain_id']=spa_data.obs['domain_id'].astype('category')


#cluSC=convert_celltype_name_into_array(sc_cluster,sc_CTname)
seq_data.obs['cell_type']=seq_data.obs['cluster']
seq_data.obs['domain_id']=0
seq_data.obs['domain_id']=seq_data.obs['domain_id'].astype('category')
seq_data.obs['Method']='RNA'
#seq_data.obs['labels']=sc_cluster[:,1



selected = filter_genes(seq_data, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# filter the object
adata_ref = seq_data[:, selected].copy()


print('\n\nI am here')



# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='domain_id',
                        # cell type, covariate used for constructing signatures
                        labels_key='cell_type',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['Method']
                       )


# create the regression model
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=250)


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)


# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']



# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(spa_data.var_names, inf_aver.index)
adata_vis = spa_data[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="domain_id")



N_cells_per_location=1
detection_alpha=20
detection_alpha=200



# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=1,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()



mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
#plt.legend(labels=['full data training']);



# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
print(adata_file)



#snapshot = tracemalloc.take_snapshot()
memory_after = memory_usage()[0]
time_end=time.time()

#print("total memory by malloc",display_top(snapshot))
print("total memory by profiler",memory_after-memory_before)
print("total execution time",time_end-time_start)


