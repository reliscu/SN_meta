setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/14_CT_similarity/within_datasets")

library(data.table)

source("CT_similarity_plot_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/15_CT_summary_vecs/datinfo_SCSN_meta_analysis_summary_vecs.csv") %>% na_if("") 

data_type <- "author_data"
expr_type <- "normalized_counts"
summary_type <- "mean"
n_genes <- 8596

cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "EXC", "INH")

sim_paths <- Sys.glob(path=file.path("data", data_type, expr_type, "*.RDS"))
sim_list <- lapply(sim_paths, readRDS)
