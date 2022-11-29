setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/14_eigencell_similarity/between_datasets")

library(data.table)

source("eigencell_similarity_plot_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_cell_type_annotations/neuronal_subtype_annotations.csv")

data_type <- "author_data"
expr_type <- "normalized_counts"
sim_type <- "pearson"

df <- fread(list.files(path=file.path("data", data_type, expr_type), pattern="MNN", full.names=T), data.table=F)
