# data/ap.py
import scanpy as sc

def adata_preprocessor(file_path, n_genes_min=1000, n_genes_max=10000):
    adata = sc.read_h5ad(file_path)
    sc.pp.filter_cells(adata, min_genes=n_genes_min)
    sc.pp.filter_cells(adata, max_genes=n_genes_max)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    return adata
