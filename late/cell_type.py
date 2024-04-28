import anndata
adata=anndata.read_h5ad("late-pathology_adata.h5ad")
with open("late_cell_type.txt", 'w') as ct:
    cell_type=adata.obs["broad.cell.type"].tolist()
    ct.write(" ".join(map(str, cell_type)) + "\n")
