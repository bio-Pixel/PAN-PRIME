suppressPackageStartupMessages({
  library(scProportionTest)
  library(Seurat)
})

ptest <- function(seu, n_permutations = 1000, cluster_identity = "Ann_Level2", sample_identity = "TME") {
  tme_levels <- unique(as.character(seu[[sample_identity]][, 1]))
  tme_levels <- tme_levels[!is.na(tme_levels)]

  res <- lapply(tme_levels, function(x) {
    seu2 <- seu
    seu2[[sample_identity]] <- as.character(seu2[[sample_identity]][, 1])
    seu2[[sample_identity]][seu2[[sample_identity]][, 1] != x, 1] <- "other"

    prop_test <- sc_utils(seu2)

    prop_test <- permutation_test(
      prop_test,
      n_permutations = n_permutations,
      cluster_identity = cluster_identity,
      sample_1 = "other",
      sample_2 = x,
      sample_identity = sample_identity
    )

    list(
      target = x,
      prop_test = prop_test,
      plot = permutation_plot(prop_test)
    )
  })

  names(res) <- tme_levels
  res
}
