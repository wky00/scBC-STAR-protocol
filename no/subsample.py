import scanpy as sc
import numpy as np

adata = sc.read('no-pathology_adata.h5ad')
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, layer='counts')

sc.pp.filter_genes(adata, min_counts=1)
sc.pp.log1p(adata, layer="counts")
stratify_key = "broad.cell.type" 
stratify_values = adata.obs[stratify_key].unique()
stratified_sample = []
for value in stratify_values:
    subset = adata[adata.obs[stratify_key] == value]
    n_cells_per_stratum = max(1, len(subset) // 10)
    if len(subset) >= n_cells_per_stratum: 
        subset_obs_index = list(subset.obs.index)
        np.random.shuffle(subset_obs_index)
        sampled_indices = subset_obs_index[:n_cells_per_stratum]
        stratified_sample.extend(sampled_indices)
adata = adata[stratified_sample]

sc.pp.highly_variable_genes(adata, n_top_genes=2000, layer="counts") 

adata = adata[:, adata.var['highly_variable']]

adata.write("no-pathology_adata.h5ad")
