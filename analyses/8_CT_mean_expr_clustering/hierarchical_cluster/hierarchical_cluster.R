setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/8_CT_mean_expr_clustering/hierarchical_cluster")

source("hiearchical_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
which_genes <- "intersection"
prop_metric <- "rho"
n_clusters <- 9

ct_expr <- fread(list.files(path=paste0("../data/", data_type, "/", expr_type), full.names=T), data.table=F)

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

top_n <- NULL
sim_type <- "pearson"
clust_method <- "ward.D2"

# for(top_n in c("NULL", 500)){
#   
#   if(top_n=="NULL"){
#     top_n <- NULL
#   }
#   
  for(clust_method in c("complete", "ward.D2", "average")){
    #for(sim_type in c("pearson", "prop")){
      plot_dendro(ct_expr, datinfo, data_type, expr_type, which_genes, gene_list, sim_type, prop_metric, clust_method, top_n, n_clusters)
   # }
  }
#   
# }

