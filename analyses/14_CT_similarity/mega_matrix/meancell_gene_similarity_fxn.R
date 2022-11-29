library(dplyr)
library(data.table)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

meancell_gene_similarity <- function(summary_vecs, sim_type){
  
  if(sim_type=="pearson"){
    
    sim_mat <- cor(t(summary_vecs[,-c(1)]))
    colnames(sim_mat) <- rownames(sim_mat) <- summary_vecs[,c(1)]
    
  }
  
  sim_df <- reshape2::melt(sim_mat)
  sim_df <- sim_df[sim_df$Var1!=sim_df$Var2,]
  sim_df <- sim_df[order(sim_df$value, decreasing=T),]
  
  sim_mat <- sim_mat[upper.tri(sim_mat)]
  hist(sim_mat)
  
}