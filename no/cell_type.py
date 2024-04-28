import anndata
adata=anndata.read_h5ad("no-pathology_adata.h5ad")
with open("no_cell_type.txt", 'w') as ct:
    cell_type=adata.obs["broad.cell.type"].tolist()
    ct.write(" ".join(map(str, cell_type)) + "\n")
