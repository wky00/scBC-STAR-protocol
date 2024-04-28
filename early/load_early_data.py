import anndata as ad
import numpy as np
early_data = ad.read_h5ad("early-pathology_adata.h5ad")
from scBC.model import scBC
my_model = scBC(adata=early_data)
my_model.train_VI()
my_model.get_reconst_data(n_samples=10)
# my_model.get_edge()
edge=np.load("early_edge_arr.npy")
my_model.edge=edge
my_model.Biclustering(L=6)
with open('early_model_output.txt', 'w') as file:
    for i, entry in enumerate(my_model.S, start=1):
        file.write(f"bc{i}\n")
        strong_class_idx = np.where(np.array(my_model.strong_class()) == i-1)[0].tolist()
        file.write(" ".join(map(str, strong_class_idx)) + "\n")
        gene_names = early_data.var_names
        gene_names_clustered = [gene_names[idx] for idx in entry['measurements']]
        file.write(" ".join(gene_names_clustered) + "\n")
