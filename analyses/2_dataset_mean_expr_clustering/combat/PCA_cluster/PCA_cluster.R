setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/2_dataset_mean_expr_clustering/combat/PCA_cluster")

source("../PCA_cluster/PCA_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv")

data_type <- "author_data"
expr_type <- "normalized_counts"

dataset_expr <- fread(list.files(path=paste0("../data/", data_type, "/", expr_type), full.names=T), data.table=F)

pc_x="PC1"; pc_y="PC2"

top_n <- NULL

# for(top_n in c("NULL", 500)){
# 
#   if(top_n=="NULL"){
#     top_n <- NULL
#   }
  
plot_PCs(dataset_expr, datinfo, data_type, expr_type, top_n, pc_x="PC1", pc_y="PC2")
plot_PCs(dataset_expr, datinfo, data_type, expr_type, top_n, pc_x="PC2", pc_y="PC3")
plot_PCs(dataset_expr, datinfo, data_type, expr_type, top_n, pc_x="PC3", pc_y="PC4")

for(which_pc in c("PC1", "PC2", "PC3")){
  PC_vs_covariates(dataset_expr, datinfo, data_type, expr_type, top_n, which_pc)
}
  
# }

