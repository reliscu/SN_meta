setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/6_cell_class_mean_expr_clustering/PCA_cluster")

source("PCA_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("")

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_cell_type_annotations/neuronal_subtype_annotations.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
top_n <- NULL
pc_x="PC1"; pc_y="PC2"
comat <- F

fid_paths <- list.files(path=paste0("../data/", data_type, "/", expr_type), full.names=T)
expr_paths <- ""

for(i in 1:length(fid_paths)){
  
  cell_class <- strsplit(sapply(strsplit(fid_paths[i], "/"), function(x) x[length(x)]), "_")[[1]][1]
  
  cat("\n")
  print(fid_paths[i])
  print(expr_paths[i])
  
  class_fid <- fread(fid_paths[i], data.table=F)
  class_expr <- fread(expr_paths[i], data.table=F)
  
  plot_PCs(class_fid, datinfo, cell_class, data_type, expr_type, pc_x="PC1", pc_y="PC2", top_n, class_expr, prop_scaled, neu_subtypes)
  plot_PCs(class_fid, datinfo, cell_class, data_type, expr_type, pc_x="PC2", pc_y="PC3", top_n, class_expr, prop_scaled, neu_subtypes)
  
}

