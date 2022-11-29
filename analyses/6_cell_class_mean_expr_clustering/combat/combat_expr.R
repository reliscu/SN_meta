setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/6_cell_class_mean_expr_clustering/combat")

source("/home/rebecca/SCSN_meta_analysis/code/combat_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv")

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

data_type <- "author_data"
expr_type <- "normalized_counts"

expr_paths <- list.files(path=paste0("../data/", data_type, "/", expr_type), pattern="union", full.names=T)

## Restrict to cortical datasets:

datinfo <- datinfo %>% dplyr::filter(grepl("C$", Region_Code))

combat_class_fxn(datinfo, expr_paths, data_type, expr_type, gene_list)