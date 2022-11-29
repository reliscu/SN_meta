setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/2_dataset_mean_expr_clustering/hierarchical_cluster")

source("hiearchical_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 

data_type <- "author_data"
expr_type <- "normalized_counts"
which_genes <- "intersection"
prop_metric <- "rho"
sim_type <- "pearson"
clust_method <- "average"

dataset_expr <- fread(list.files(path=paste0("../data/", data_type, "/", expr_type), full.names=T), data.table=F)

for(clust_method in c("complete", "ward.D2", "average")){
  
  plot_dendro(dataset_expr, datinfo, data_type, expr_type, combat=F, sim_type, prop_metric, clust_method, top_n)

}




