setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/19_cell_modeling/mega_matrix")

source("summary_cell_fxns.R")

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_CT_annotations/neuronal_subtype_annotations.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
n_resamples <- 10
summary_type <- "mean"
cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "EXC", "INH", "CGE", "MGE")

cellinfo <- fread("/home/rebecca/SCSN_meta_analysis/prelim_new/16_mega_matrix/data/author_data/normalized_counts/cellinfo_author_data_normalized_counts_421993_nuclei_30_datasets_8596_intersection_genes.csv", data.table=F)

expr <- readRDS("/home/rebecca/SCSN_meta_analysis/prelim_new/16_mega_matrix/data/author_data/normalized_counts/expr_author_data_normalized_counts_421993_nuclei_30_datasets_8596_intersection_genes.RDS")

############################################# Calc. class summary vectors: SCALED #############################################

scaled <- T

# calc_class_summary_vecs(expr, cellinfo, neu_subtypes, data_type, expr_type, summary_type, cell_classes, scaled)

# calc_class_random_vecs(expr, cellinfo, neu_subtypes, data_type, expr_type, summary_type, cell_classes, n_resamples)

############################################# Calc. class summary vectors: UNSCALED #############################################

scaled <- F

# calc_class_summary_vecs(expr, cellinfo, neu_subtypes, data_type, expr_type, summary_type, cell_classes, scaled)

############################################# Calc. class LEFT OUT dataset summary vectors: SCALED #############################################

scaled <- T

# calc_class_summary_vecs_leftout(expr, cellinfo, neu_subtypes, data_type, expr_type, summary_type, cell_classes, scaled)