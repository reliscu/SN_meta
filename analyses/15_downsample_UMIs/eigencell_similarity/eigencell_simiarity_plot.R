setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/eigencell_similarity")

source("eigencell_simiarity_plot_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
sim_type <- "pearson"

sim_list <- readRDS(list.files(path=file.path("data", data_type, expr_type), pattern="stats", full.names=T))

df <- do.call(rbind, sim_list)


