setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/14_eigencell_similarity/between_datasets")

source("eigencell_similarity_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
sim_type <- "pearson"

eigencell_paths <- list.files(path=file.path("../data", data_type, expr_type), pattern="eigencell", full.names=T)
expr_paths <- list.files(path=file.path("../data", data_type, expr_type), pattern="max_expr", full.names=T)

datasets1 <- sapply(strsplit(sapply(strsplit(eigencell_paths, "/"), "[", 5), "_eigencells"), "[", 1)
datasets2 <- sapply(strsplit(sapply(strsplit(expr_paths, "/"), "[", 5), "_cell_max"), "[", 1)

all.equal(datasets1, datasets2)

datinfo <- datinfo[match(datasets1, datinfo$Dataset),]

#eigencell_pairwise_similarity(datinfo, eigencell_paths, data_type, expr_type, gene_list, sim_type)

## Calc similarity stats:

eigen_sim_paths <- list.files(path=file.path("data", data_type, expr_type), pattern="_eigencell_vs", full.names=T)

datasets <- sapply(strsplit(sapply(strsplit(eigen_sim_paths, "/"), "[", 4), "_eigencell"), "[", 1)

datinfo <- datinfo[match(datasets, datinfo$Dataset),]

#eigencell_similarity_stats(datinfo, eigen_sim_paths, data_type, expr_type)

############################################# MNN #############################################

stats_list <- readRDS(list.files(path=file.path("data", data_type, expr_type), pattern="^eigencell_vs", full.names=T))

############################################# Real vs. rand #############################################