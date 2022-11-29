setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/9_CT_fidelity_clustering/hierarchical_cluster")

source("hiearchical_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
which_genes <- "intersection"
prop_metric <- "rho"
n_clusters <- 9

prop_scaled <- T

fid_pattern <- "fidelity_author"
if(prop_scaled){
  fid_pattern <- "prop"
}

ct_fid <- fread(list.files(path=paste0("../data/", data_type, "/", expr_type), pattern=fid_pattern, full.names=T), data.table=F)
ct_expr <- fread(list.files(path=paste0("../../8_CT_mean_expr_clustering/data/", data_type, "/", expr_type), full.names=T), data.table=F)

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

top_n <- NULL
sim_type <- "pearson"
clust_method <- "average"

for(clust_method in c("complete", "ward.D2", "average")){
  plot_dendro(ct_fid, datinfo, data_type, expr_type, which_genes, gene_list, sim_type, prop_metric, clust_method, top_n, ct_expr, n_clusters, prop_scaled)
}
