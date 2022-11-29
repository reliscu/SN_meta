plot_knn_model_boxplot <- function(model_res_paths){
  
  model_frac <- lapply(1:length(model_res_paths), function(i){
    model_res <- read.csv(model_res_paths[i])
    frac_per_k <- unlist(lapply(grep("k", colnames(model_res)), function(k){
      sum(model_res[,k]>.05/nrow(model_res))/nrow(model_res)
    }))
    df <- data.frame(
      Dataset=sapply(strsplit(sapply(
        strsplit(model_res_paths[i], "/"), function(x) x[length(x)]
      ), "_gene_knn", ), "[", 1),
      No.Neighbors=colnames(model_res)[grep("k", colnames(model_res))], 
      Frac=frac_per_k
    )
    return(df)
  })
  
  model_df <- do.call(rbind, model_frac)
  model_df$No.Neighbors <- gsub("k_", "", model_df$No.Neighbors)
  model_df$No.Neighbors <- factor(model_df$No.Neighbors, levels=sort(as.numeric(unique(model_df$No.Neighbors))))
  
  ggplot(model_df, aes(x=No.Neighbors, y=Frac, group=No.Neighbors)) +
    geom_boxplot()
  
}

# plot_nn_v_chance <- function(){
#   
#   nn_v_chance_frac <- lapply(1:length(nn_v_chance_paths), function(i){
#     nn_res <- read.csv(nn_v_chance_paths[i])
#     sum(nn_res$Empirical_Pval>0)
#     df <- data.frame(
#       Dataset=sapply(strsplit(sapply(
#         strsplit(nn_v_chance_paths[i], "/"), function(x) x[length(x)]
#       ), "_gene_knn", ), "[", 1),
#       No.Neighbors=colnames(model_res)[grep("k", colnames(model_res))], 
#       No.Genes=sum(nn_res$Empirical_Pval>0)
#     )
#     return(df)
#   })
#   
# }