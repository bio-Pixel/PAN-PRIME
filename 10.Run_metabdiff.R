## Compare KEGG pathway activity shifts between TME-associated cell subsets
## Effect size (standardized difference) and significance are evaluated per pathway
## using permutation-defined case/control cell groups

suppressPackageStartupMessages({
  library(decoupleR)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(Seurat)
  library(parallel)
})

# load cell-level metadata (cell IDs, annotations, TME labels)
meta <- readRDS("METAINFO_ALLCELL.rds")

# load table defining differential TME–cell-type combinations
tme <- read.delim("TME_DIFF_ABUNDANCE.txt")

# build cell pools by different annotation levels
cpool  <- split(meta[["cell"]], meta[["Ann_Level1"]])   # level-1 cell types
cpool2 <- split(meta[["cell"]], meta[["Ann_Level2"]])   # level-2 cell subtypes
cpoolt <- split(meta[["cell"]], meta[["TME"]])          # TME groups

message("1")
# load KEGG pathway activity scores inferred by decoupleR (ULM)
dat <- readRDS("KEGG_DECOUPLER_SCORES.rds")

# split scores by pathway and convert to named numeric vectors (cell-level)
dat <- split(dat, dat[["source"]])
dat <- lapply(dat, function(x) {
  v <- x[["score"]]
  names(v) <- x[["condition"]]  # cell IDs
  v
})

message("2")
# initialize parallel backend
cl <- makeCluster(30)

message("3")
# iterate over each TME–cell-type combination
re <- parApply(cl, tme, 1, function(x, cpool, cpool2, cpoolt, dat) {

  # extract identifiers for current comparison
  tme_id  <- x[["TME"]]
  cell_id <- x[["Ann_Level2"]]
  l1_id   <- x[["Ann_Level1"]]

  # define case cells: cells belonging to both the target TME and cell subtype
  case <- intersect(cpoolt[[tme_id]], cpool2[[cell_id]])

  # define control cells: remaining cells of the same level-1 lineage
  control <- setdiff(cpool[[l1_id]], case)

  # compute pathway-wise effect sizes and statistical significance
  do.call(rbind, lapply(dat, function(y) {

    y_case <- y[case]
    y_control <- y[control]

    # summary statistics
    mean_1 <- mean(y_case, na.rm = TRUE)
    mean_2 <- mean(y_control, na.rm = TRUE)
    sd_1 <- sd(y_case, na.rm = TRUE)
    sd_2 <- sd(y_control, na.rm = TRUE)

    # standardized difference (Cohen's d–like metric)
    pooled_sd <- sqrt((sd_1^2 + sd_2^2) / 2)
    delta_std <- ifelse(
      is.na(pooled_sd) || pooled_sd == 0,
      0,
      (mean_1 - mean_2) / pooled_sd
    )

    # non-parametric significance testing (robust to non-normality)
    if (length(na.omit(y_control)) < 2 || length(na.omit(y_case)) == 0) {
      p <- 1
    } else {
      p <- tryCatch(
        wilcox.test(y_case, y_control)[["p.value"]],
        error = function(e) 1
      )
    }

    data.frame(delta_std = delta_std, p = p)
  }))
}, cpool, cpool2, cpoolt, dat)

# shut down parallel backend
stopCluster(cl)
