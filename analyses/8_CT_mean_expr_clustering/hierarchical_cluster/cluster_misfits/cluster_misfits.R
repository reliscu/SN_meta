setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/8_CT_mean_expr_clustering/hierarchical_cluster/cluster_misfits")

source("cluster_misfits_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")
cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")

data_type <- "author_data"
expr_type <- "normalized_counts"
which_genes <- "intersection"
prop_metric <- "rho"

n_clusters <- 9

ct_expr <- fread(list.files(path=paste0("../../data/", data_type, "/", expr_type), pattern="union", full.names=T), data.table=F)

top_n <- NULL
sim_type <- "pearson"
clust_method <- "ward.D2"

model_cluster_misfits(
  ct_expr,
  datinfo,
  data_type,
  expr_type,
  which_genes,
  sim_type,
  prop_metric,
  clust_method,
  top_n,
  n_clusters
)

plot_cluster_misfits(
  ct_expr,
  datinfo,
  data_type,
  expr_type,
  which_genes,
  sim_type,
  prop_metric,
  clust_method,
  top_n,
  n_clusters
)

plot_highly_connected_class_members(
  ct_expr,
  datinfo,
  data_type,
  expr_type,
  cell_classes,
  which_genes,
  sim_type,
  prop_metric,
  top_n
)