setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/2_dataset_mean_expr_clustering/PCA_cluster")

source("PCA_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 

data_type <- "author_data"
expr_type <- "normalized_counts"

dataset_expr <- fread(list.files(path=paste0("../data/", data_type, "/", expr_type), full.names=T), data.table=F)

combat <- F
top_n <- NULL
pc_x="PC1"; pc_y="PC2"

plot_PCs(dataset_expr, datinfo, data_type, expr_type, combat, gene_list, top_n, pc_x="PC1", pc_y="PC2")
plot_PCs(dataset_expr, datinfo, data_type, expr_type, combat, gene_list, top_n, pc_x="PC2", pc_y="PC3")

for(which_pc in c("PC1", "PC2", "PC3")){
  PC_vs_covariates(dataset_expr, datinfo, data_type, expr_type, combat, gene_list, top_n, which_pc)
}
