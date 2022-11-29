setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/8_CT_mean_expr_clustering/hierarchical_cluster/neuron_subtypes")

source("neuron_subtypes_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_cell_type_annotations/neuronal_subtype_annotations.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
which_genes <- "intersection"
prop_metric <- "rho"
n_clusters <- 9

ct_expr <- fread(list.files(path=paste0("../../data/", data_type, "/", expr_type), full.names=T), data.table=F)

top_n <- NULL
clust_method <- "ward.D2"
sim_type <- "pearson"

plot_neuron_subtypes(ct_expr, datinfo, data_type, expr_type, which_genes, gene_list, sim_type, prop_metric, clust_method, top_n, n_clusters, neu_subtypes)
