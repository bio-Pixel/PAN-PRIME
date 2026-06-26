source("src/consensus_graph_modules.R")

g <- readRDS("spatial_module_graph.rds")

res <- consensus_graph_modules(
  g,
  resolutions = seq(0.2, 2, by = 0.1),
  nrun = 100,
  algorithm = "leiden"
)

library(igraph)
library(ggraph)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

C <- res$consensus_matrix
cols <- readRDS("Network_layout/SpaM_cols.rds")

lab <- factor(
  paste0("SpaM", res$consensus_membership),
  levels = paste0("SpaM",1:14)
)

hc <- hclust(as.dist(1 - C), method = "ward.D2")

ord <- order(lab, match(seq_along(lab), hc$order))

C2 <- C[ord, ord]
lab2 <- lab[ord]

ha <- HeatmapAnnotation(
  cluster = lab2,
  col = list(cluster = cols)
)

pdf(
  "consensusLeiden/consensusLeiden_matrix.pdf",
  width = 7,
  height = 7
)

Heatmap(
  C2,
  name = "Consensus",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  rect_gp = gpar(col = NA),
  col = colorRamp2(c(0, 1), c("white", "blue")),
  use_raster = TRUE,
  raster_quality = 50,
  width = unit(8, "cm"),
  height = unit(8, "cm"),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 1)
)

dev.off()


lay <- create_layout(g, layout = "fr")
lay$SpaM = paste0("SpaM",sm[lay$Sample_cluster])

p1 <- ggraph(lay) +
  geom_edge_link( aes(width = weight), alpha = 0.2,color = "grey80") +
  geom_node_point( aes(color = SpaM, size = degree(g)) ,
                   alpha = 0.8) +
  scale_edge_width(range = c(0.2, 3)) +
  scale_size(range = c(1, 5)) +
  scale_color_manual(values = cols)+
  theme_void()+
  theme(aspect.ratio = 0.8)
ggsave(paste0("consensusLeiden/SpaM_network.pdf"), width = 8, height = 6)

