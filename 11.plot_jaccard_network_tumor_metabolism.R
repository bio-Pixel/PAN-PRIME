# ---- Palettes (lightened for softer look) ----
tme_colors <- colorspace::lighten(
  c(TME1 = "#cc6666", TME2 = "#6699cc", TME3 = "#008280",
    TME4 = "#bb0021", TME5 = "#631879"),
  amount = 0.5
)
names(tme_colors) <- paste0("TME", 1:5)

tumor_colors <- colorspace::lighten(
  c(K1 = "#007B7F", K2 = "#DF75AE", K3 = "#007ABA", K4 = "#A8C7E9",
    K5 = "#00B7CA", K6 = "#F37121", K7 = "#8CCA7C", K8 = "#8B58A4"),
  amount = 0.5
)
names(tumor_colors) <- paste0("K", 1:8)

#' Build an undirected graph from a Jaccard matrix
#'
#' @param jaccard_mat Numeric symmetric matrix with row/col names like "Kx_TMEy".
#' @param threshold Numeric in [0,1]. Keep edges with Jaccard > threshold.
#' @param min_degree Integer. Drop nodes with degree < min_degree after thresholding.
#'
#' @return tbl_graph object for downstream ggraph plotting.
build_jaccard_graph <- function(jaccard_mat, threshold = 0.4, min_degree = 1) {
  stopifnot(is.matrix(jaccard_mat), is.numeric(jaccard_mat))
  if (is.null(rownames(jaccard_mat)) || is.null(colnames(jaccard_mat))) {
    stop("jaccard_mat must have rownames and colnames.")
  }
  if (!all(rownames(jaccard_mat) == colnames(jaccard_mat))) {
    stop("Row names and column names must match and be in the same order.")
  }
  if (!isTRUE(all.equal(jaccard_mat, t(jaccard_mat), tolerance = 1e-8))) {
    warning("jaccard_mat is not perfectly symmetric; forcing symmetry by (M + t(M))/2.")
    jaccard_mat <- (jaccard_mat + t(jaccard_mat)) / 2
  }

  # Lower-triangle long form (exclude diagonal)
  lower_tri <- lower.tri(jaccard_mat, diag = FALSE)
  jac_long <- data.frame(
    from = rownames(jaccard_mat)[row(jaccard_mat)[lower_tri]],
    to   = colnames(jaccard_mat)[col(jaccard_mat)[lower_tri]],
    jac  = jaccard_mat[lower_tri],
    stringsAsFactors = FALSE
  )
  jac_long <- dplyr::filter(jac_long, jac > threshold)

  # Node attributes: extract Tumor Kx and TMEx from labels like "K3_TME4"
  nodes <- unique(c(jac_long$from, jac_long$to))
  tumor <- sub("^(K\\d+)_.*$", "\\1", nodes)
  tme   <- sub(".*_(TME\\d+)$", "\\1", nodes)
  node_df <- data.frame(name = nodes, Tumor = tumor, TME = tme, stringsAsFactors = FALSE)

  # Build graph and drop isolated nodes (degree < min_degree) if requested
  g <- igraph::graph_from_data_frame(jac_long, vertices = node_df, directed = FALSE)
  tg <- tidygraph::as_tbl_graph(g)

  if (min_degree > 0) {
    tg <- tg |>
      tidygraph::mutate(deg = tidygraph::centrality_degree(mode = "all")) |>
      tidygraph::filter(deg >= min_degree) |>
      tidygraph::select(-deg)
  }
  tg
}

#' Plot a Jaccard network with ggraph
#'
#' @param tg tbl_graph from build_jaccard_graph().
#' @param layout Character ggraph layout, e.g., "kk", "fr", "stress".
#' @param color_by One of c("TME", "Tumor"). Controls node color mapping.
#' @param edge_width_range Numeric length-2 vector for edge width range.
#' @param node_size Numeric node size.
#' @param label Logical, whether to show node labels.
#'
#' @return A ggplot object.
plot_jaccard_network <- function(
  tg,
  layout = "kk",
  color_by = c("TME", "Tumor"),
  edge_width_range = c(0.2, 2),
  node_size = 6,
  label = TRUE
) {
  color_by <- match.arg(color_by)

  # Choose palette based on color_by
  pal <- switch(color_by,
                TME   = tme_colors,
                Tumor = tumor_colors)

  # Prepare a safe manual scale: any unseen category becomes gray
  scale_node_color <- ggplot2::scale_color_manual(
    values = pal,
    na.value = "grey75",
    guide = ggplot2::guide_legend(override.aes = list(size = 4))
  )

  p <- ggraph::ggraph(tg, layout = layout) +
    ggraph::geom_edge_link(ggplot2::aes(width = .data$jac),
                           alpha = 0.5, color = "grey80") +
    ggraph::geom_node_point(ggplot2::aes(color = .data[[color_by]]),
                            size = node_size) +
    ggplot2::scale_edge_width(range = edge_width_range) +
    scale_node_color +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right", aspect.ratio = 1.2)

  if (isTRUE(label)) {
    p <- p + ggraph::geom_node_text(ggplot2::aes(label = .data$name),
                                    repel = TRUE, size = 3)
  }
  p
}

# ---- Example usage ----
# tg <- build_jaccard_graph(jaccard_mat, threshold = 0.4, min_degree = 1)
# plot_jaccard_network(tg, layout = "kk", color_by = "TME")
# plot_jaccard_network(tg, layout = "fr", color_by = "Tumor")
