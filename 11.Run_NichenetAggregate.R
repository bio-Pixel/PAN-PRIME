run_nichenet_aggregate <- function(
  sender_rds,
  receiver_rds,
  sender_label,
  receiver_label,
  groupby_col = "Ann_Level2",
  donor_id_col = "DonorID",
  subtype_rds,
  condition_col = "subtype",
  condition_oi = "M5",
  condition_ref = "other",
  expression_pct = 0.25,
  lr_network_rds,
  ligand_target_matrix_rds,
  weighted_networks_rds
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(nichenetr)
  })

  # load and subset sender
  sender <- readRDS(sender_rds)
  sender <- subset(sender, subset = .data[[groupby_col]] %in% sender_label)

  # load and subset receiver
  receiver <- readRDS(receiver_rds)
  receiver <- subset(receiver, subset = .data[[groupby_col]] %in% receiver_label)

  # merge sender and receiver
  seuratObj <- merge(sender, receiver)
  Idents(seuratObj) <- groupby_col

  # attach condition labels (e.g., TME subtype per donor)
  subtype_map <- readRDS(subtype_rds)
  seuratObj@meta.data[[condition_col]] <-
    subtype_map[as.character(seuratObj@meta.data[[donor_id_col]])]

  # load NicheNet priors
  lr_network <- readRDS(lr_network_rds) %>% distinct(from, to)
  ligand_target_matrix <- readRDS(ligand_target_matrix_rds)
  weighted_networks <- readRDS(weighted_networks_rds)

  # run NicheNet aggregate analysis
  nichenet_output <- nichenet_seuratobj_aggregate(
    seurat_obj = seuratObj,
    sender = sender_label,
    receiver = receiver_label,
    condition_colname = condition_col,
    condition_oi = condition_oi,
    condition_reference = condition_ref,
    expression_pct = expression_pct,
    ligand_target_matrix = ligand_target_matrix,
    lr_network = lr_network,
    weighted_networks = weighted_networks
  )

  list(
    seurat_obj = seuratObj,
    nichenet_output = nichenet_output,
    conditions = unique(seuratObj@meta.data[[condition_col]])
  )
}
