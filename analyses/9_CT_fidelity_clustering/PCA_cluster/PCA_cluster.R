setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/9_CT_fidelity_clustering/PCA_cluster")

source("PCA_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")
 
gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_cell_type_annotations/neuronal_subtype_annotations.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"

prop_scaled <- T

fid_pattern <- "fidelity_author"
if(prop_scaled){
  fid_pattern <- "prop"
}

ct_fid <- fread(list.files(path=paste0("../data/", data_type, "/", expr_type), pattern=fid_pattern, full.names=T), data.table=F)
ct_expr <- fread(list.files(path=paste0("../../8_CT_mean_expr_clustering/data/", data_type, "/", expr_type), full.names=T), data.table=F)

top_n <- NULL
pc_x="PC1"; pc_y="PC2"

plot_PCs(ct_fid, datinfo, data_type, expr_type, pc_x="PC1", pc_y="PC2", gene_list, top_n, ct_expr, prop_scaled, neu_subtypes)
plot_PCs(ct_fid, datinfo, data_type, expr_type, pc_x="PC2", pc_y="PC3", gene_list, top_n, ct_expr, prop_scaled, neu_subtypes)
