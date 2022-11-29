setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/6_cell_class_mean_expr_clustering/combat/PCA_cluster")

source("../../PCA_cluster/PCA_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("")

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_cell_type_annotations/neuronal_subtype_annotations.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
top_n <- NULL
pc_x="PC1"; pc_y="PC2"
comat <- F

expr_paths <- list.files(path=paste0("../data/", data_type, "/", expr_type), full.names=T)

for(i in 1:length(expr_paths)){
  
  cell_class <- strsplit(sapply(strsplit(expr_paths[i], "/"), function(x) x[length(x)]), "_")[[1]][2]
  class_expr <- fread(expr_paths[i], data.table=F)
  
  plot_PCs(class_expr, datinfo, cell_class, data_type, expr_type, combat, pc_x="PC1", pc_y="PC2", gene_list, top_n, neu_subtypes)
  plot_PCs(class_expr, datinfo, cell_class, data_type, expr_type, combat, pc_x="PC2", pc_y="PC3", gene_list, top_n, neu_subtypes)
  
} 
