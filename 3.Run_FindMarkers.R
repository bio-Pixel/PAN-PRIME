## Differential expression analysis for cluster-level markers
## The analysis is performed on a pre-defined Seurat object (cell type already subset upstream)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# load Seurat object (cell type–restricted upstream)
aa <- readRDS("INPUT_SEURAT_OBJECT.rds")

# identify cluster markers
markers <- FindAllMarkers(
  object = aa,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  only.pos = TRUE
)

# rank markers by effect size
markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]

# filter by detection rate in the target cluster
markers <- markers[markers$pct.1 > 0.25, ]

return(markers)
