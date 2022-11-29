setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/cell_similarity_distribution")

source("cell_similarity_distribution_hist_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
sim_type <- "cosine"

cell_sim_paths <- list.files(path=file.path("data", data_type, expr_type), pattern=toupper(sim_type), full.names=T)


