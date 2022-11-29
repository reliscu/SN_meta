setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/10_mega_matrix_fidelity")

source("mega_matrix_fidelity_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 

data_type <- "author_data"
expr_type <- "normalized_counts"
prop_scaled <- F
cell_classes <- c("ASC", "END", "NEU", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")

cellinfo <- fread("/home/rebecca/SCSN_meta_analysis/prelim_new/16_mega_matrix/data/author_data/normalized_counts/cellinfo_author_data_normalized_counts_421993_nuclei_30_datasets_8596_intersection_genes.csv", data.table=F)

expr <- readRDS("/home/rebecca/SCSN_meta_analysis/prelim_new/16_mega_matrix/data/author_data/normalized_counts/expr_author_data_normalized_counts_421993_nuclei_30_datasets_8596_intersection_genes.RDS")

# cell_class_fidelity(expr, cellinfo, data_type, expr_type, cell_classes, prop_scaled=F)