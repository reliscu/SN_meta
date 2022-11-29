setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/2_dataset_mean_expr_clustering/combat/hierarchical_cluster")

source("../hierarchical_cluster/hiearchical_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv")

data_type <- "author_data"
expr_type <- "normalized_counts"
which_genes <- "intersection"
prop_metric <- "rho"

dataset_expr <- fread(list.files(path=paste0("../data/", data_type, "/", expr_type), full.names=T), data.table=F)

## Restrict to protein coding genes:

pc_genes <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/protein_coding_19723_union_genes.csv", data.table=F)

dataset_expr <- dataset_expr[is.element(dataset_expr$SYMBOL, pc_genes$SYMBOL),]

sim_type <- "pearson"
clust_method <- "ward.D2"
top_n <- NULL

for(clust_method in c("complete", "ward.D2", "average")){

 # for(sim_type in c("pearson", "prop")){

    plot_dendro(dataset_expr, datinfo, data_type, expr_type, which_genes, sim_type, prop_metric, clust_method, top_n)

 # }

}
