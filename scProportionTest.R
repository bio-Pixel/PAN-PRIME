library("scProportionTest")
ptest <- function(seu){
  df = lapply(as.character(unique(seu$TME)), function(x){
    seu@meta.data$TME = as.character(seu@meta.data$TME)
    seu@meta.data$TME[-which(seu$TME == x)] = "other"
    
    prop_test <- sc_utils(seu)
    
    prop_test <- permutation_test(n_permutations = 1000,
      prop_test, cluster_identity = "Ann_Level2",
      sample_1 = "other", sample_2 = x,
      sample_identity = "TME"
    )
    p = permutation_plot(prop_test)
  })
}
