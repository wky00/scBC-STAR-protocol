import anndata as ad
import numpy as np
no_data = ad.read_h5ad("no-pathology_adata.h5ad")
from scBC.model import scBC
my_model = scBC(adata=no_data)
my_model.train_VI()
my_model.get_reconst_data(n_samples=10)
# my_model.get_edge()
edge=np.load("no_edge_arr.npy")
my_model.edge=edge
my_model.Biclustering(L=6)
#print(my_model.strong_class())
with open('no_model_output.txt', 'w') as file:
    for i, entry in enumerate(my_model.S, start=1):
        file.write(f"bc{i}\n")
        strong_class_idx = np.where(np.array(my_model.strong_class()) == i-1)[0].tolist()
        file.write(" ".join(map(str, strong_class_idx)) + "\n")
        gene_names = no_data.var_names
        gene_names_clustered = [gene_names[idx] for idx in entry['measurements']]
        file.write(" ".join(gene_names_clustered) + "\n")
