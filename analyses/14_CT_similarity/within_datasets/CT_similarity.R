setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/14_CT_similarity/within_datasets")

source("CT_similarity_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") %>% dplyr::filter(!is.na(CT_Mean_Vecs))

data_type <- "author_data"
expr_type <- "normalized_counts"
summary_type <- "mean"
sim_type <- "pearson"
n_resamples <- 10

#cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "EXC", "INH")

#CT_real_vs_rand_similarity(datinfo, data_type, expr_type, sim_type, cell_classes, summary_type, n_resamples)
