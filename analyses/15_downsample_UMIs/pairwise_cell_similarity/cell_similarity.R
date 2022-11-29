setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/cell_similarity_distribution")

source("cell_similarity_distribution_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

expr_list_paths <- list.files(path="../data/author_data/QC_counts", full.names=T)

data_type <- "author_data"
expr_type <- "normalized_counts"
cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")

datasets <- sapply(strsplit(sapply(strsplit(expr_list_paths, "/"), "[", 5), "_UMIs"), "[", 1)

datinfo <- datinfo[match(datasets, datinfo$Dataset),]

sim_type <- "pearson"
