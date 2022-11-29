setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity")

source("fidelity_fxns.R")

datinfo <- read.csv("../15_CT_summary_vecs/datinfo_SCSN_meta_analysis_summary_vecs.csv") %>% na_if("") %>% dplyr::filter(!is.na(Author_Barcode_Annotations_Class))

data_type <- "author_data"
expr_type <- "normalized_counts"
cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "EXC", "INH")
n_resamples <- 10

prop_scaled <- T

# CT_fidelity(datinfo, data_type, expr_type, cell_classes, prop_scaled)
# cell_class_fidelity(datinfo, data_type, expr_type, cell_classes, prop_scaled)



## Add datinfo paths:

datinfo <- read.csv("../15_CT_summary_vecs/datinfo_SCSN_meta_analysis_summary_vecs.csv") %>% na_if("") 

file_path <- file.path(getwd(), Sys.glob(paste0("data/", data_type, "/", expr_type, "/*CT_fidelity_", data_type, "*")))
datasets <- sapply(strsplit(sapply(strsplit(file_path, "/"), "[", 10), "_CT"), "[", 1)
idx <- match(datasets, datinfo$Dataset)
datinfo$CT_Fidelity <- NA
datinfo$CT_Fidelity[idx] <- file_path

#fwrite(datinfo, file="datinfo_SCSN_meta_analysis_fidelity.csv")

datinfo <- read.csv("datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 

file_path <- file.path(getwd(), Sys.glob(paste0("data/", data_type, "/", expr_type, "/*CT_fidelity_prop*")))
datasets <- sapply(strsplit(sapply(strsplit(file_path, "/"), "[", 10), "_CT"), "[", 1)
idx <- match(datasets, datinfo$Dataset)
datinfo$CT_Fidelity_Prop_Scaled <- NA
datinfo$CT_Fidelity_Prop_Scaled[idx] <- file_path

#fwrite(datinfo, file="datinfo_SCSN_meta_analysis_fidelity.csv")
