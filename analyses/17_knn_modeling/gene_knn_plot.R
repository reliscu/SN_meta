setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/13_gene_nearest_neighbors")

library(dplyr)
library(data.table)

source("gene_knn_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

pc_genes <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/protein_coding_19723_union_genes.csv", data.table=F)

data_type <- "author_data"
expr_type <- "normalized_counts"

model_res_paths <- list.files(path=file.path("data", data_type, expr_type), pattern="modeling", full.names=T)

