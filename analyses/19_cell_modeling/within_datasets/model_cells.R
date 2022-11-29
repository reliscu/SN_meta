setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/19_cell_modeling/within_datasets")

source("model_cells_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") %>% dplyr::filter(!is.na(Cell_Class_Mean_Vecs))

data_type <- "author_data"
expr_type <- "normalized_counts"
n_resamples <- 10
summary_type <- "mean"

############################################# Real vs. random cell class modeling ############################################# 

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/14_neuron_subtypes_annotations/neuron_subtype_annotations_917_CTs_18_datasets.csv") %>% na_if("")

# cell_classes <- c("ASC", "MIC", "OG", "OPC", "NEU") 
# cell_classes <- c("ASC", "MIC", "OG", "OPC", "EXC", "INH") 
cell_classes <- c("ASC", "MIC", "OG", "OPC", "EXC", "CGE", "MGE") 

# model_class_real_vs_random(datinfo, data_type, expr_type, summary_type, cell_classes, neu_subtypes, n_resamples)

############################################# Real vs. random CT modeling ############################################# 

# model_CT_real_vs_random(datinfo, data_type, expr_type, summary_type, n_resamples)

############################################# Cell class vs. CT modeling ############################################# 

# cell_classes <- c("ASC", "MIC", "OG", "OPC", "NEU") 