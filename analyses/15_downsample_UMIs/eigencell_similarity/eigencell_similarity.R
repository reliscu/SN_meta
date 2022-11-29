setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/eigencell_similarity")

source("eigencell_similarity_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
sim_type <- "pearson"

expr_list_paths <- list.files("../data/author_data/QC_counts", pattern="downsampled", full.names=T)
eigencell_list_paths <- list.files(path=file.path("data/author_data/", expr_type), pattern="downsampled", full.names=T)

datasets1 <- sapply(strsplit(sapply(strsplit(expr_list_paths, "/"), "[", 5), "_UMIs"), "[", 1)
datasets2 <- sapply(strsplit(sapply(strsplit(eigencell_list_paths, "/"), "[", 5), "_eigencells"), "[", 1)
all.equal(datasets1, datasets2)

datinfo <- datinfo[match(datasets1, datinfo$Dataset),]

#eigencell_similarity(atinfo, expr_list_paths, eigencell_list_paths, data_type, expr_type, sim_type, top_n=NULL)

## Calc eigencell similarity stats:

eigen_sim_paths <- list.files(path=file.path("data", data_type, expr_type), pattern="cell_vs_eigencell", full.names=T)

datasets <- sapply(strsplit(sapply(strsplit(eigen_sim_paths, "/"), "[", 4), "_cell_vs_eigencell"), "[", 1)

datinfo <- datinfo[match(datasets, datinfo$Dataset),]

eigencell_sim_list <- lapply(eigen_sim_paths, function(x) fread(x, data.table=F))

#eigencell_similarity_stats(datinfo, eigencell_sim_list, data_type, expr_type)


