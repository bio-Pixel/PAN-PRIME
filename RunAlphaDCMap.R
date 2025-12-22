## Diffusion map analysis based on CSS embeddings for pancreatic Alpha cells
## All computations are performed on Alpha-cell–restricted data for each donor/patient independently

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(simspec)
  library(destiny)
})

# input/output placeholders (paths hidden for sharing)
in_rds  <- "ALPHA_CELL_OBJECT.rds"
out_css <- "ALPHA_CELL_CSS_OBJECT.rds"
out_dm  <- "ALPHA_CELL_CSS_DCmap.rds"

# load Alpha cell Seurat object
alpha <- readRDS(in_rds)

# standard Seurat preprocessing
alpha <- NormalizeData(alpha, normalization.method = "LogNormalize", scale.factor = 10000)
alpha <- FindVariableFeatures(alpha, selection.method = "vst", nfeatures = 5000)
alpha <- ScaleData(alpha, features = VariableFeatures(alpha))
alpha <- RunPCA(alpha, features = VariableFeatures(alpha))

# cluster similarity spectrum (CSS) embedding
alpha <- cluster_sim_spectrum(
  object = alpha,
  label_tag = "id",
  spectrum_type = "corr_ztransform",
  corr_method = "spearman"
)

# diffusion map and diffusion pseudotime (DPT) computed on CSS embeddings
css <- Embeddings(alpha, "css")
dmm <- DiffusionMap(css)
dpt <- DPT(dmm)

alpha_dcmap <- data.frame(
  eigenvectors(dmm)[, 1:3],
  time = dpt[["dpt"]],
  cell_id = rownames(css),
  louvain = alpha@meta.data[["louvain"]],
  row.names = rownames(css)
)

# save Alpha-cell diffusion map results
saveRDS(alpha_dcmap, out_dm)

