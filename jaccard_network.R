
tme_colors = colorspace::lighten(c(TME1 = "#cc6666", TME2 = "#6699cc", TME3 = "#008280", 
                                   TME4 = "#bb0021", TME5 = "#631879"), amount = 0.5)
names(tme_colors) = paste0("TME",1:5)
col_tumor <- colorspace::lighten(c(
  K1 = "#007B7F", K2 = "#DF75AE", K3 = "#007ABA", K4 = "#A8C7E9", K5 = "#00B7CA",
  K6 = "#F37121", K7 = "#8CCA7C", K8 = "#8B58A4"
), amount = 0.5)
names(col_tumor) = paste0("K",1:8)

plot_jaccard_network_custom <- function(jaccard_mat, threshold = 0.4) {
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(dplyr)
  library(ggplot2)
  library(colorspace)
 
  # Convert matrix to lower-triangle long format
  lower_tri <- lower.tri(jaccard_mat, diag = FALSE)
  jac_long <- data.frame(
    from = rownames(jaccard_mat)[row(jaccard_mat)[lower_tri]],
    to   = colnames(jaccard_mat)[col(jaccard_mat)[lower_tri]],
    jac  = jaccard_mat[lower_tri]
  ) %>% filter(jac > threshold)
  
  # Create node table
  nodes <- unique(c(jac_long$from, jac_long$to))
  TME <- gsub(".*_(TME\\d+)$", "\\1", nodes)
  Tumour <- gsub("^(K\\d+)_.*$", "\\1", nodes)
  node_df <- data.frame(name = nodes, TME = TME, Tumour = Tumour)
  
  # Graph object
  g <- graph_from_data_frame(jac_long, vertices = node_df, directed = FALSE)

}


# Plot
ggraph(plot_jaccard_network_custom(jaccard_mat), layout = 'kk') +
  geom_edge_link(aes(width = jac), color = "gray80", alpha = 0.5) +
  geom_node_point(aes(color = TME), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE ) +
  scale_color_manual(values = tme_colors) +
  scale_edge_width(range = c(0.2, 2)) +
  theme_void() +
  theme(legend.position = "right",
        aspect.ratio = 1.2)
