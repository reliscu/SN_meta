setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/bulk_vs_SN_fidelity")

source("cell_class_fidelity_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
cell_classes <- c("ASC", "END", "NEU", "MIC", "OG", "OPC", "PER", "VSMC")

expr_list_paths <- list.files(path="../data/author_data/QC_counts", full.names=T)

datasets <- sapply(strsplit(sapply(strsplit(expr_list_paths, "/"), "[", 5), "_UMIs"), "[", 1)

datinfo <- datinfo[match(datasets, datinfo$Dataset),]

prop_scaled <- F

# cell_class_fidelity(
#   datinfo, 
#   expr_list_paths,
#   data_type, 
#   expr_type,
#   prop_scaled
# )
