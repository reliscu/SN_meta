setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/CT_PC1_variance_explained")

source("CT_PC1_variance_explained_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

pc_genes <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/protein_coding_19723_union_genes.csv", data.table=F)

data_type <- "author_data"
expr_type <- "normalized_counts"

## Restrict to cortical datasets:

datinfo <- datinfo %>% dplyr::filter(grepl("C$", Region_Code))

expr_list_paths <- list.files("../data/author_data/QC_counts", pattern="downsampled", full.names=T)

datasets <- sapply(strsplit(sapply(strsplit(expr_list_paths, "/"), "[", 5), "_UMIs"), "[", 1)

datinfo <- datinfo[match(datasets, datinfo$Dataset),]

#CT_PC1_VE(datinfo, expr_list_paths, data_type, expr_type)