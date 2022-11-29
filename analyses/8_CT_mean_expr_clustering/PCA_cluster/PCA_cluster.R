setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/8_CT_mean_expr_clustering/PCA_cluster")

source("PCA_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") 

data_type <- "author_data"
expr_type <- "normalized_counts"

datinfo <- datinfo %>% na_if("") 

ct_expr <- fread(list.files(path=paste0("../data/", data_type, "/", expr_type), full.names=T), data.table=F)

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_cell_type_annotations/neuronal_subtype_annotations.csv") %>% na_if("")

top_n <- NULL

# pc_x="PC1"; pc_y="PC2"

plot_PCs(ct_expr, datinfo, data_type, expr_type, pc_x="PC1", pc_y="PC2", gene_list, top_n, neu_subtypes)
plot_PCs(ct_expr, datinfo, data_type, expr_type, pc_x="PC2", pc_y="PC3", gene_list, top_n, neu_subtypes)
