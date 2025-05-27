import streamlit as st


import streamlit as st
from PIL import Image
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import ap
import hdf5plugin
import scanorama

import os


import streamlit as st





def data_preprocessor(adata):
    st.write("📦 Κανονικοποίηση, Log1p, Scaling και PCA...")
    sc.pp.normalize_total(adata, target_sum=1e4)  # Κανονικοποίηση
    sc.pp.log1p(adata)  # Logarithm μετασχηματισμός
    adata.raw = adata  # Αποθηκεύει τα raw δεδομένα
    sc.pp.scale(adata, max_value=10)  # Standardize τα χαρακτηριστικά
    sc.pp.pca(adata)  # PCA
    return adata


st.set_page_config(page_title="scRNA-seq Analyzer", layout="wide")

# Sidebar
st.sidebar.title("⚙️ Επιλογές")
page = st.sidebar.radio("Μενού:", [
    "🏠 Αρχική",
    "📁 Φόρτωση Δεδομένων",
    "🧪 Προεπεξεργασία",
    "🔗 Ενοποίηση (Scanorama)",
    "📊 DEG Ανάλυση",
    "🌋 Volcano Plot",
    "🧬 Gene Expression Plots",
    "👥 Πληροφορίες Ομάδας"
])


if page == "🏠 Αρχική":
    st.title("🔬 Εφαρμογή Ανάλυσης scRNA-seq Δεδομένων")
    st.markdown("""
    Καλώς ήρθατε στην Streamlit εφαρμογή για ανάλυση δεδομένων μοριακής βιολογίας τύπου single-cell RNA-seq. 

    Χρησιμοποιήστε το μενού αριστερά για να περιηγηθείτε στο pipeline.
    """)

elif page == "📁 Φόρτωση Δεδομένων":


    st.subheader("🧬 Δεδομένα:")
    load_mode = st.radio("Τύπος αρχείου:", ["Αρχείο .h5ad"])

    if load_mode == "Αρχείο .h5ad":
        uploaded_file = st.file_uploader("📤 Επιλέξτε αρχείο .h5ad", type=["h5ad"])

        if uploaded_file is not None:
            try:
                adata = ad.read_h5ad(uploaded_file)
                st.success("✅ Το αρχείο ανέβηκε και διαβάστηκε επιτυχώς!")
                st.write(f"📐 Σχήμα: {adata.shape[0]} κύτταρα × {adata.shape[1]} γονίδια")
                st.dataframe(adata.obs.head())
            except Exception as e:
                st.error(f"❌ Σφάλμα κατά την ανάγνωση του αρχείου: {e}")







elif page == "🧪 Προεπεξεργασία":
    st.header("🧪 Προεπεξεργασία Δεδομένων")

    directory = "./data/h5ad"
    processed_dir = "./data/h5ad_filt"
    os.makedirs(processed_dir, exist_ok=True)

    st.info("Εδώ μπορείτε να φιλτράρετε όλα τα αρχεία .h5ad στον φάκελο ./data/h5ad")

    # Λίστα αρχείων h5ad
    h5ad_files = [f for f in os.listdir(directory) if f.endswith(".h5ad")]
    if not h5ad_files:
        st.warning(f"Δεν βρέθηκαν αρχεία .h5ad στον φάκελο: {directory}")
    else:
        st.write(f"Βρέθηκαν τα αρχεία: {h5ad_files}")

        if st.button("Εκτέλεση φίλτρων σε όλα τα αρχεία"):
            import hdf5plugin 
            import ap 

            for file in h5ad_files:
                file_path = os.path.join(directory, file)
                st.write(f"Επεξεργασία αρχείου: {file_path}")

                try:
                    adata_filtered = ap.adata_preprocessor(file_path, n_genes_min=1000, n_genes_max=10000)
                    output_path = os.path.join(processed_dir, f"filtered_{file}")
                    adata_filtered.write_h5ad(output_path, compression=hdf5plugin.FILTERS["zstd"])
                    st.success(f"✅ Αποθηκεύτηκε το φιλτραρισμένο αρχείο: {output_path}")
                except Exception as e:
                    st.error(f"❌ Σφάλμα στην επεξεργασία του αρχείου {file}: {e}")



elif page == "🔗 Ενοποίηση (Scanorama)":
    st.header("🔗 Ενοποίηση Δεδομένων με Scanorama")
    st.info("Ανεβάστε τουλάχιστον 2 αρχεία .h5ad για ενοποίηση με Scanorama.")

    uploaded_files = st.file_uploader(
        "Ανεβάστε πολλαπλά αρχεία (.h5ad):",
        type=["h5ad"],
        accept_multiple_files=True,
    )

    if not uploaded_files:
        st.warning("⚠️ Παρακαλώ ανεβάστε τουλάχιστον 2 αρχεία για ενοποίηση.")
    else:
        if len(uploaded_files) < 2:
            st.error("❌ Χρειάζονται τουλάχιστον 2 αρχεία για ενοποίηση.")
        else:
            if st.button("🚀 Εκτέλεση Scanorama Ενοποίησης"):
                try:
                    adatas = []
                    for file in uploaded_files:
                        adata = sc.read_h5ad(file)
                        st.write(f"📄 {file.name}: Διαστάσεις δεδομένων {adata.shape}")

                        # Εδώ προσθέτεις τη στήλη condition
                        if "control" in file.name.lower():
                            adata.obs["condition"] = "control"
                        elif "disease" in file.name.lower():
                            adata.obs["condition"] = "disease"
                        else:
                            adata.obs["condition"] = "unknown"




                        # Ελαφριά προεπεξεργασία για Scanorama
                        sc.pp.normalize_total(adata, target_sum=1e4)
                        sc.pp.log1p(adata)
                        sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor='seurat_v3')
                        sc.pp.pca(adata, n_comps=50)

                        adatas.append(adata)

                    shared_genes = set.intersection(*(set(ad.var_names) for ad in adatas))
                    st.write(f"🔬 Κοινά γονίδια: {len(shared_genes)}")
                    if len(shared_genes) == 0:
                        st.error("❌ Δεν υπάρχουν κοινά γονίδια μεταξύ των αρχείων!")
                        st.stop()

                    adatas = [adata[:, list(shared_genes)] for adata in adatas]

                    st.write("🧠 Εκτέλεση ενοποίησης με Scanorama...")
                    scanorama.integrate_scanpy(adatas, dimred=50)

                    scanorama_integrated = [adata.obsm["X_scanorama"] for adata in adatas]
                    integrated_matrix = np.concatenate(scanorama_integrated)

                    combined_adata = adatas[0].concatenate(adatas[1:], join='inner')



                    st.write("Batch IDs στο combined_adata.obs:")
                    st.write(combined_adata.obs['batch'].unique())




                    
                    condition_map = {'0': 'control', '1': 'disease'}
                    combined_adata.obs['condition'] = combined_adata.obs['batch'].map(condition_map).astype(str)

                    st.write("Unique values στο condition:")
                    st.write(combined_adata.obs['condition'].unique())



                    combined_adata.obsm["Scanorama"] = integrated_matrix


                    st.write("📌 Υπολογισμός UMAP στα ενοποιημένα δεδομένα...")
                    sc.pp.neighbors(combined_adata, use_rep="Scanorama")
                    sc.tl.umap(combined_adata)

                    import matplotlib.pyplot as plt

                    # Δημιουργία UMAP χωρίς εμφάνιση
                    sc.pl.umap(combined_adata, color=["batch"], title="Batch Visualization", show=False)

                    # Εμφάνιση στο Streamlit
                    st.pyplot(plt.gcf())

                    
                    plt.savefig("./data/h5ad_integrat/umap_batch_visualization.png", dpi=300)

                    # Εκκαθάριση για επόμενη χρήση
                    plt.clf()

                    os.makedirs("./data/h5ad_integrat", exist_ok=True)
                    output_file = "./data/h5ad_integrat/integrated_scanorama.h5ad"
                    combined_adata.write_h5ad(output_file)
                    st.success(f"✅ Ενοποιημένα δεδομένα αποθηκεύτηκαν στο: `{output_file}`")

                except Exception as e:
                    st.error(f"❌ Σφάλμα: {e}")







elif page == "📊 DEG Ανάλυση":
    st.header("📊 Ανάλυση Διαφορικής Έκφρασης")
    st.info("Επιλέξτε ομάδες κυττάρων για σύγκριση.")

    ADATA_PATH ="./data/h5ad_integrat/integrated_scanorama.h5ad"


    if not os.path.exists(ADATA_PATH):
        st.error(f"❌ Το αρχείο {ADATA_PATH} δεν βρέθηκε. Παρακαλώ συγχώνευσε πρώτα τα αρχεία.")
    else:
        try:
            adata = sc.read_h5ad(ADATA_PATH)

            # Αν δεν έχει γίνει log1p, το κάνουμε
            if "log1p" not in adata.uns:
                sc.pp.log1p(adata)

            # Ελέγχουμε αν υπάρχει η στήλη condition για groupby
            if "condition" not in adata.obs.columns:
                st.error("❌ Η στήλη 'condition' δεν βρέθηκε στο .obs. Χρειάζεται για το groupby.")
            else:
                sc.tl.rank_genes_groups(
                    adata,
                    groupby="condition",
                    groups=["disease"],
                    reference="control",
                    method="wilcoxon"
                )

                result = adata.uns["rank_genes_groups"]

                degs = pd.DataFrame({
                    "genes": result["names"]["disease"],
                    "pvals": result["pvals"]["disease"],
                    "pvals_adj": result["pvals_adj"]["disease"],
                    "logfoldchanges": result["logfoldchanges"]["disease"],
                })

                st.subheader("🧾 Αποτελέσματα DEG:")
                st.dataframe(degs)

                degs_filtered = degs[
                    (degs["pvals"] <= 0.05) &
                    (degs["pvals"] != 0.0) &
                    (degs["logfoldchanges"].abs() > 0.5)
                ].reset_index(drop=True)

                st.subheader("✅ Φιλτραρισμένα DEG (p-val <=0.05 & |logFC| > 0.5):")
                st.dataframe(degs_filtered)

                out_dir = "./data/deg_data"
                os.makedirs(out_dir, exist_ok=True)
                out_path = os.path.join(out_dir, "alzheimer_data_degs.csv")
                degs_filtered.to_csv(out_path, index=False)
                st.success(f"✅ Φιλτραρισμένα DEG αποθηκεύτηκαν στο: {out_path}")

        except Exception as e:
            st.error(f"❌ Σφάλμα κατά την ανάλυση DEG: {e}")







elif page == "🌋 Volcano Plot":
    st.header("🌋 Volcano Plot")
    st.info("Οπτικοποίηση DEG αποτελεσμάτων.")

    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    import pandas as pd
    import numpy as np

    DEG_CSV_PATH = "./data/deg_data/alzheimer_data_degs.csv"

    if not os.path.exists(DEG_CSV_PATH):
        st.error(f"❌ Το αρχείο DEG δεν βρέθηκε: {DEG_CSV_PATH}")
    else:
        df = pd.read_csv(DEG_CSV_PATH)
        if df.empty:
            st.warning("⚠️ Το αρχείο DEG είναι κενό.")
        else:
            # Υπολογισμός -log10(p-value)
            df["neg_log10_pval"] = -np.log10(df["pvals"])

            # Κατηγοριοποίηση διαφοροποιημένης έκφρασης
            df["diffexpressed"] = "NS"
            df.loc[(df["logfoldchanges"] > 1) & (df["pvals"] < 0.001), "diffexpressed"] = "UP"
            df.loc[(df["logfoldchanges"] < -1) & (df["pvals"] < 0.001), "diffexpressed"] = "DOWN"

            plt.figure(figsize=(10, 6))
            sns.scatterplot(
                data=df,
                x="logfoldchanges",
                y="neg_log10_pval",
                hue="diffexpressed",
                palette={"UP": "#bb0c00", "DOWN": "#00AFBB", "NS": "grey"},
                alpha=0.7,
                edgecolor=None,
                legend="full",
            )

            # Όρια 
            plt.axhline(y=-np.log10(0.05), color='gray', linestyle='dashed')
            plt.axvline(x=-1, color='gray', linestyle='dashed')
            plt.axvline(x=1, color='gray', linestyle='dashed')

            plt.xlim(df["logfoldchanges"].min() - 1, df["logfoldchanges"].max() + 1)
            plt.ylim(0, df["neg_log10_pval"].max() + 5)
            plt.xlabel("log2 Fold Change", fontsize=14)
            plt.ylabel("-log10 p-value", fontsize=14)
            plt.title("Volcano Plot of DEGs (Disease vs Control)", fontsize=16)
            plt.legend(title="Expression", loc="upper right")

            st.pyplot(plt.gcf())








elif page == "🧬 Gene Expression Plots":
    st.header("🧬 Gene Expression Plots")
    st.info("Plot έκφρασης γονιδίων για επιλεγμένα genes.")

    import os
    import matplotlib.pyplot as plt
    import scanpy as sc

    H5AD_PATH = "./data/h5ad_integrat/integrated_scanorama.h5ad"

    if not os.path.exists(H5AD_PATH):
        st.error(f"❌ Το αρχείο δεν βρέθηκε: {H5AD_PATH}")
    else:
        adata = sc.read_h5ad(H5AD_PATH)

        # Έλεγχος για UMAP
        if "X_umap" not in adata.obsm:
            st.write("Δεν βρέθηκε embedding UMAP. Υπολογίζω...")
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
        else:
            st.write("Υπάρχει ήδη UMAP embedding.")

        # Εμφάνιση πρώτων 20 γονιδίων
        st.write("📋 Πρώτα 20 γονίδια:")
        st.write(list(adata.var_names)[:20])

        # Επιλογή γονιδίων με dropdown
        default_genes = ['MALAT1', 'KCNIP4', 'DPP10']
        selected_genes = st.multiselect(
            "🔍 Επιλέξτε γονίδια προς απεικόνιση:",
            options=sorted(adata.var_names),
            default=[g for g in default_genes if g in adata.var_names]
        )

        if len(selected_genes) == 0:
            st.warning("⚠️ Παρακαλώ επιλέξτε τουλάχιστον ένα γονίδιο.")
        else:
            rows = (len(selected_genes) + 2) // 3
            columns = 3

            def plot_gene_expr(adata, gene_lst, rows, columns):
                fig, axs = plt.subplots(rows, columns, figsize=(columns * 4, rows * 4))
                axs = axs.flatten()
                for i, gene in enumerate(gene_lst):
                    if i < len(axs):
                        if gene in adata.var_names:
                            sc.pl.umap(
                                adata,
                                color=gene,
                                add_outline=True,
                                legend_loc='on data',
                                legend_fontsize=12,
                                legend_fontoutline=2,
                                frameon=True,
                                title=f'{gene}',
                                palette='Set1',
                                ax=axs[i],
                                show=False
                            )
                        else:
                            axs[i].set_title(f"{gene} (not found)")
                            axs[i].axis("off")
                for j in range(i + 1, len(axs)):
                    fig.delaxes(axs[j])
                plt.tight_layout()
                return fig

            fig = plot_gene_expr(adata, selected_genes, rows, columns)
            st.pyplot(fig)






elif page == "👥 Πληροφορίες Ομάδας":
    st.header("👥 Πληροφορίες Ομάδας")
    st.markdown("""
    **Ονόματα Μελών:**
    - Ελένη Στύλου – UI , Latex , Φόρτωση Δεδομένων, Gene Expression Plots
    - Ελευθερία Δημολιού – Scanorama, DEG Ανάλυση, Volcano Plot, Docker
    - Παρασκευή Σιαμπανάι– Latex , Προεπεξεργασία










    **GitHub Repo:** https://github.com/Eldimoliou/scRNA-seq-Analyzer
    """)
