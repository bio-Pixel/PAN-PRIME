# inferCNV analysis is performed separately for each individual patient

import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import pandas as pd

sc.settings.set_figure_params(dpi=80, figsize=(5, 5))

# load AnnData object
adata = sc.read("INPUT_DATA.h5ad")

# load gene location table (hg38)
# df should be a DataFrame indexed by gene names and
# containing genomic coordinates required by infercnvpy
df = pd.DataFrame(...)  # gene location table (hidden)

# align genes between AnnData and gene-location table
common_genes = adata.var_names.intersection(df.index)
adata = adata[:, common_genes].copy()
df = df.loc[common_genes]

# attach genomic coordinates to adata.var
adata.var = df

# define reference cell populations
exclude_list = ["KRT19+", "Acinar", "Ductal", "Tuft", "ADM"]
unique_elements = adata.obs["Ann_Level1"].unique().tolist()
reference_cat = [x for x in unique_elements if x not in exclude_list]

# infer copy-number variation
cnv.tl.infercnv(
    adata,
    reference_key="Ann_Level1",
    reference_cat=reference_cat,
    window_size=250,
)

# downstream embedding and CNV scoring
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)

# save results
adata.write_h5ad("OUTPUT_CNV.h5ad")

# visualize chromosome-level CNV heatmap
cnv.pl.chromosome_heatmap(adata, groupby="Ann_Level1", show=False)
plt.savefig("CNV_heatmap.png", dpi=200, bbox_inches="tight")
plt.close()
