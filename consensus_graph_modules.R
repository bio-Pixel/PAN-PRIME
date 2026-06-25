library(igraph)
library(mclust)

consensus_graph_modules <- function(
    g,
    resolutions = seq(0.2, 2, by = 0.1),
    nrun = 100,
    algorithm = c("leiden", "louvain"),
    weight_col = "weight",
    seed = 123,
    min_size_cutoff = 3,
    modularity_tol = 0.995,
    ari_cutoff = 0.80,
    prefer_resolution = NULL
){
  
  algorithm <- match.arg(algorithm)
  set.seed(seed)
  
  w <- edge_attr(g, weight_col)
  if(is.null(w)) stop("No edge weight found: ", weight_col)
  
  n <- vcount(g)
  node_names <- V(g)$name
  if(is.null(node_names)) node_names <- as.character(seq_len(n))
  
  A <- as_adj(g, attr = weight_col, sparse = FALSE)
  
  calc_metrics <- function(memb){
    
    cls <- unique(memb)
    
    within_w <- sum(sapply(cls, function(k){
      ids <- which(memb == k)
      if(length(ids) <= 1) return(0)
      sum(A[ids, ids]) / 2
    }))
    
    total_w <- sum(w)
    coverage <- within_w / total_w
    
    conductance <- mean(sapply(cls, function(k){
      ids <- which(memb == k)
      other <- setdiff(seq_along(memb), ids)
      
      cut_w <- sum(A[ids, other])
      vol_in <- sum(A[ids, ])
      vol_out <- sum(A[other, ])
      
      cut_w / min(vol_in, vol_out)
    }), na.rm = TRUE)
    
    data.frame(
      modularity = modularity(g, memb, weights = w),
      n_cluster = length(unique(memb)),
      min_size = min(table(memb)),
      max_size = max(table(memb)),
      coverage = coverage,
      conductance = conductance
    )
  }
  
  run_once <- function(res){
    if(algorithm == "leiden"){
      membership(cluster_leiden(
        g,
        weights = w,
        resolution_parameter = res,
        objective_function = "modularity"
      ))
    } else {
      membership(cluster_louvain(
        g,
        weights = w,
        resolution = res
      ))
    }
  }
  
  consensus_matrix <- function(mem.list){
    C <- matrix(0, n, n)
    for(m in mem.list){
      C <- C + outer(m, m, "==")
    }
    C <- C / length(mem.list)
    diag(C) <- 1
    rownames(C) <- colnames(C) <- node_names
    C
  }
  
  pairwise_ari_mean <- function(mem.list){
    cmb <- combn(seq_along(mem.list), 2)
    mean(apply(cmb, 2, function(x){
      adjustedRandIndex(mem.list[[x[1]]], mem.list[[x[2]]])
    }))
  }
  
  get_consensus_membership <- function(C, k){
    hc <- hclust(as.dist(1 - C), method = "average")
    cutree(hc, k = k)
  }
  
  calc_pac <- function(C, lower = 0.1, upper = 0.9){
    x <- C[upper.tri(C)]
    mean(x > lower & x < upper)
  }
  
  # 1. resolution gradient + consensus
  res.list <- lapply(resolutions, function(res){
    
    mem.list <- lapply(seq_len(nrun), function(i){
      run_once(res)
    })
    
    mods <- sapply(mem.list, function(m){
      modularity(g, m, weights = w)
    })
    
    nclus <- sapply(mem.list, function(m){
      length(unique(m))
    })
    
    C <- consensus_matrix(mem.list)
    k <- round(median(nclus))
    cons.mem <- get_consensus_membership(C, k = k)
    
    met <- calc_metrics(cons.mem)
    
    stat <- data.frame(
      resolution = res,
      mean_run_modularity = mean(mods),
      sd_run_modularity = sd(mods),
      mean_run_ARI = pairwise_ari_mean(mem.list),
      PAC = calc_pac(C),
      consensus_modularity = met$modularity,
      n_cluster = met$n_cluster,
      min_size = met$min_size,
      max_size = met$max_size,
      coverage = met$coverage,
      conductance = met$conductance
    )
    
    list(
      resolution = res,
      stat = stat,
      consensus_matrix = C,
      consensus_membership = cons.mem,
      mem.list = mem.list
    )
  })
  
  stat.df <- do.call(rbind, lapply(res.list, function(x) x$stat))
  
  # 2. no-resolution reference algorithms
  ref.list <- list()
  
  ref.list$fast_greedy <- tryCatch({
    membership(cluster_fast_greedy(g, weights = w))
  }, error = function(e) NULL)
  
  ref.list$walktrap <- tryCatch({
    membership(cluster_walktrap(g, weights = w))
  }, error = function(e) NULL)
  
  ref.list$infomap <- tryCatch({
    membership(cluster_infomap(g, e.weights = w))
  }, error = function(e) NULL)
  
  ref.list <- ref.list[!sapply(ref.list, is.null)]
  
  if(length(ref.list) > 0){
    stat.df$ARI_to_reference <- sapply(seq_along(res.list), function(i){
      m <- res.list[[i]]$consensus_membership
      mean(sapply(ref.list, function(ref){
        adjustedRandIndex(m, ref)
      }))
    })
  } else {
    stat.df$ARI_to_reference <- NA_real_
  }
  
  # 3. choose stable plateau, not highest score
  max_mod <- max(stat.df$consensus_modularity, na.rm = TRUE)
  
  cand <- stat.df[
    stat.df$consensus_modularity >= max_mod * modularity_tol &
      stat.df$mean_run_ARI >= ari_cutoff &
      stat.df$min_size >= min_size_cutoff,
  ]
  
  if(nrow(cand) == 0){
    cand <- stat.df[
      stat.df$consensus_modularity >= max_mod * modularity_tol &
        stat.df$min_size >= min_size_cutoff,
    ]
  }
 
  if(!is.null(prefer_resolution)){
    cand$rank_value <- abs(cand$resolution - prefer_resolution)
    cand <- cand[order(
      cand$rank_value,
      cand$n_cluster,
      cand$conductance,
      -cand$consensus_modularity
    ), ]
  } else {
    cand <- cand[order(
      cand$n_cluster,
      cand$conductance,
      -cand$consensus_modularity,
      cand$resolution
    ), ]
  }
  
  best.row <- cand[1, ]
  best.idx <- which(stat.df$resolution == best.row$resolution)[1]
  best <- res.list[[best.idx]]
  
  V(g)$consensus_cluster <- best$consensus_membership
  
  list(
    best_resolution = best.row$resolution,
    best_stat = best.row,
    stat = stat.df,
    candidate_resolutions = cand,
    consensus_membership = best$consensus_membership,
    consensus_matrix = best$consensus_matrix,
    reference_memberships = ref.list,
    graph = g,
    all_results = res.list
  )
}
