setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/14_CT_similarity/between_datasets")

source("CT_similarity_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/15_CT_summary_vecs/datinfo_SCSN_meta_analysis_summary_vecs.csv") %>% na_if("") %>% dplyr::filter(!is.na(Author_Barcode_Annotations_Class))

data_type <- "author_data"
expr_type <- "normalized_counts"
n_resamples <- 10
summary_type <- "mean"
sim_type <- "pearson"

real_sim <- readRDS(Sys.glob(paste0("data/", data_type, "/", expr_type, "/*CT_", summary_type, "*")))

rand_sim <- readRDS(Sys.glob(paste0("data/", data_type, "/", expr_type, "/*random_", summary_type, "*")))

