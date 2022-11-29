setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/6_cell_class_mean_expr_clustering/combat/hierarchical_cluster")

source("../hierarchical_cluster/hiearchical_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 

data_type <- "author_data"
expr_type <- "normalized_counts"
which_genes <- "intersection"
prop_metric <- "rho"
sim_type <- "pearson"
clust_method <- "average"

expr_paths <- list.files(path=paste0("../data/", data_type, "/", expr_type), full.names=T)

combat <- T

for(i in 1:length(expr_paths)){

  cell_class <- strsplit(sapply(strsplit(expr_paths[i], "/"), function(x) x[length(x)]), "_")[[1]][1]

  print(cell_class)

  class_expr <- fread(expr_paths[i], data.table=F)
  
  for(clust_method in c("complete", "ward.D2", "average")){
    plot_dendro(class_expr, datinfo, cell_class, data_type, expr_type, combat, sim_type, prop_metric, clust_method, top_n)
  }

}

