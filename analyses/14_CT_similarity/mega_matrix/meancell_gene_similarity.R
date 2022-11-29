setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/14_CT_similarity/mega_matrix")

source("meancell_gene_similarity_fxn.R")

data_type <- "author_data"
expr_type <- "normalized_counts"
summary_type <- "mean"
sim_type <- "pearson"

cellinfo <- fread("/home/rebecca/SCSN_meta_analysis/prelim_new/16_mega_matrix/data/author_data/normalized_counts/cellinfo_author_data_normalized_counts_421993_nuclei_30_datasets_8596_intersection_genes.csv", data.table=F)

summary_vecs <- fread(Sys.glob(path=file.path("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/19_cell_modeling/mega_matrix/data", data_type, expr_type, paste0("cell_class_", summary_type, "*unscaled*"))), data.table=F)
