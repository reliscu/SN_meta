setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/2_dataset_mean_expr_clustering/combat")

source("/home/rebecca/SCSN_meta_analysis/code/combat_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv")

data_type <- "author_data"
expr_type <- "normalized_counts"

dataset_expr <- fread(list.files(path=paste0("../data/", data_type, "/", expr_type), pattern="union", full.names=T), data.table=F)

## Restrict to cortical datasets:

datinfo <- datinfo %>% dplyr::filter(grepl("C$", Region_Code))

## Restrict to intersection coding genes:

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

dataset_expr <- dataset_expr[is.element(dataset_expr$SYMBOL, gene_list),]

expr <- dataset_expr
