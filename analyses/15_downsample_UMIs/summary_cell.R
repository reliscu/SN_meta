setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs")

source("summary_cell_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("")

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_CT_annotations/neuronal_subtype_annotations.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
n_resamples <- 10
summary_type <- "mean"
scaled <- T

expr_paths <- Sys.glob(paths=paste0("data/", data_type, "/QC_counts/*"))

datasets <- sapply(strsplit(sapply(strsplit(expr_paths, "/", fixed=T), "[", 4), "_UMI"), "[", 1)

datinfo <- datinfo[match(datasets, datinfo$Dataset),]

############################################# Calc. class summary vectors #############################################

## Create summary vectors only from cell classes we're going to model with: cell classes with 1+ nuclei, present in most datasets, neu = exc/inh

cell_classes <- c("ASC", "MIC", "OG", "OPC", "EXC", "INH", "MGE", "CGE") 

#calc_class_summary_vecs(expr, cellinfo, neu_subtypes, data_type, expr_type, summary_type, cell_classes)

#calc_class_random_vecs(expr, cellinfo, neu_subtypes, data_type, expr_type, summary_type, cell_classes, n_resamples)
