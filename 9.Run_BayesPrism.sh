# Run BayesPrism deconvolution

for i in ICGC_AU ICGC_CA TCGA CPTAC GSE224564
do
  script="${i}_metacell.R"

  cat > "${script}" << EOF
## BayesPrism deconvolution for bulk PDAC dataset: ${i}
## This script loads a bulk RNA-seq count matrix and a single-cell reference (Seurat),
## harmonizes gene features, and estimates cell-type fractions with BayesPrism.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(BayesPrism)
})

# ---- input/output placeholders (paths hidden) ----
dataset_dir <- "DATASET_DIR/${i}"                      # bulk dataset directory
bulk_pattern <- "_matrix.rds"                           # pattern for bulk count matrix file
sc_ref_rds <- "SC_REFERENCE_SEURAT_OBJECT.rds"          # Seurat object used as sc reference
meta_rds <- "SC_CELL_METADATA.rds"                      # metadata infomation
out_rds <- "OUTPUT_DIR/${i}_theta.rds"                  # output file (theta + theta.cv)

# ---- load bulk count matrix ----
setwd(dataset_dir)
matrix_file <- list.files(pattern = bulk_pattern)
bk.dat <- readRDS(matrix_file)
features <- rownames(bk.dat)
bk.dat <- bk.dat[features, , drop = FALSE]

# ---- load single-cell reference and attach refined annotations ----
seu_obj <- readRDS(sc_ref_rds)
meta <- readRDS(meta_rds)

# align metadata to the Seurat object by cell barcodes
meta <- meta[meta[["cell"]] %in% colnames(seu_obj), , drop = FALSE]
rownames(meta) <- meta[["cell"]]
meta <- meta[colnames(seu_obj), , drop = FALSE]

# ---- subset to overlapping genes (bulk features) ----
suba <- subset(seu_obj, features = features)
Idents(suba) <- "Ann_Level2"
print(levels(suba))

# ---- extract single-cell count matrix and harmonize genes with bulk ----
sc_count <- GetAssayData(object = suba, assay = "RNA", slot = "counts")
features <- intersect(rownames(bk.dat), rownames(sc_count))

bk.dat <- bk.dat[features, , drop = FALSE]
sc.dat.filtered.pc.sig <- t(as.matrix(sc_count))[, features, drop = FALSE]

# cell type/state labels used by BayesPrism
cell.type.labels <- as.character(Idents(suba))
cell.state.labels <- cell.type.labels

# ---- construct prism object ----
myPrism <- new.prism(
  reference = sc.dat.filtered.pc.sig,
  mixture = t(bk.dat),
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

# ---- run BayesPrism ----
bp.res <- run.prism(prism = myPrism, n.cores = 10)

# ---- extract estimated cell-type fractions (theta) and CV ----
theta <- get.fraction(bp = bp.res, which.theta = "final", state.or.type = "type")
theta.cv <- bp.res@posterior.theta_f@theta.cv

# ---- save results ----
re <- list(theta = theta, theta.cv = theta.cv)
saveRDS(re, out_rds)
EOF
done
