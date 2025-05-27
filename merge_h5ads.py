import os
import scanpy as sc
import anndata
import hdf5plugin

# Έκδοση αρχείου
VERSION = 1

# Φάκελοι
FILTERED_DIR = "./data/h5ad_filt"
CONCAT_DIR = "./data/h5ad_concat"
MERGED_NAME = "alzheimer_data_concat"
H5AD_CONCAT = f"{CONCAT_DIR}/{MERGED_NAME}_v{VERSION}.h5ad"

# Δημιουργία φακέλων αν δεν υπάρχουν
os.makedirs(FILTERED_DIR, exist_ok=True)
os.makedirs(CONCAT_DIR, exist_ok=True)

def get_h5ads(directory=FILTERED_DIR):
    h5ad_name_lst = []
    for file in os.listdir(directory):
        if file.endswith(".h5ad"):
            h5ad_name_lst.append(os.path.join(directory, file))
    return h5ad_name_lst

def merge_h5ads():
    h5ad_paths = get_h5ads()
    if len(h5ad_paths) == 0:
        raise ValueError(f"Δεν βρέθηκαν αρχεία .h5ad στον φάκελο {FILTERED_DIR}")

    # Φορτώνουμε όλα τα αρχεία
    adatas = [sc.read_h5ad(path) for path in h5ad_paths]

    # Συγχωνεύουμε (concatenate)
    adata_merged = anndata.AnnData.concatenate(*adatas, batch_key='batch', join="inner")

    # Έλεγχος και φιλτράρισμα .var για gene_ids αν υπάρχει
    if 'gene_ids' in adata_merged.var.columns:
        adata_merged.var = adata_merged.var[['gene_ids']]
    else:
        print("Προσοχή: Δεν βρέθηκε η στήλη 'gene_ids' στο var. Προχωράμε χωρίς φιλτράρισμα.")

    # Αποθηκεύουμε το αποτέλεσμα
    adata_merged.write_h5ad(
        H5AD_CONCAT,
        compression=hdf5plugin.FILTERS["zstd"],
        compression_opts=hdf5plugin.Zstd(clevel=5).filter_options
    )
    print(f"Merge complete. Αποθηκεύτηκε στο: {H5AD_CONCAT}")

if __name__ == "__main__":
    merge_h5ads()

