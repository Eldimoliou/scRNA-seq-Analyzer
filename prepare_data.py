import scanpy as sc

adata = sc.datasets.pbmc3k()
adata.write('./data/h5ad/dataset1.h5ad')
print("Data saved!")
