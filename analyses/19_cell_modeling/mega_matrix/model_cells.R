setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/19_cell_modeling/mega_matrix")

source("model_cells_fxns.R")

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_CT_annotations/neuronal_subtype_annotations.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
n_resamples <- 10
summary_type <- "mean"

cellinfo <- fread("/home/rebecca/SCSN_meta_analysis/prelim_new/16_mega_matrix/data/author_data/normalized_counts/cellinfo_author_data_normalized_counts_421993_nuclei_30_datasets_8596_intersection_genes.csv", data.table=F)

############################################# Model all cells with real and random CELL CLASS summary vectors ############################################# 

expr <- readRDS("/home/rebecca/SCSN_meta_analysis/prelim_new/16_mega_matrix/data/author_data/normalized_counts/expr_author_data_normalized_counts_421993_nuclei_30_datasets_8596_intersection_genes.RDS")

summary_vecs <- fread(Sys.glob(path=file.path("data", data_type, expr_type, paste0("cell_class_", summary_type, "*csv")))[1], data.table=F)

random_vecs <- readRDS(Sys.glob(path=file.path("data", data_type, expr_type, paste0("random_cell_class_", summary_type, "*"))))

# cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU") 
# cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "EXC", "INH") 
cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "EXC", "INH", "CGE", "MGE")

# model_class_real_vs_random(expr, cellinfo, neu_subtypes, summary_vecs, random_vecs, data_type, expr_type, summary_type, cell_classes, n_resamples)

# model_class_real_vs_random_sequentially(expr, cellinfo, neu_subtypes, summary_vecs, random_vecs, data_type, expr_type, summary_type, cell_classes, n_resamples)

############################################# Model left out dataset with CELL CLASS summary vectors ############################################# 

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/15_CT_summary_vecs/datinfo_SCSN_meta_analysis_summary_vecs.csv") %>% na_if("") 

summary_list <- readRDS(Sys.glob(path=file.path("data", data_type, expr_type, paste0("*", summary_type, "cells_left_out*"))))

cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU") 
# cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "EXC", "INH") 
# cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "EXC", "INH", "CGE", "MGE")

# model_class_leftout_dataset(datinfo, cellinfo, summary_list, data_type, expr_type, summary_type, cell_classes)