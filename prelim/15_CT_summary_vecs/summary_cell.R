setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/15_CT_summary_vecs")

source("summary_cell_fxns.R")

datinfo <- read.csv("../13_CT_stats/datinfo_SCSN_meta_analysis_CT_stats.csv") %>% na_if("") %>% dplyr::filter(!is.na(Author_Barcode_Annotations_Class))

neu_subtypes <- read.csv("../14_neuron_subtypes_annotations/neuron_subtype_annotations_917_CTs_18_datasets.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
n_resamples <- 10

############################################# Calc. mean CELL CLASS summary vectors: SCALED #############################################

summary_type <- "mean"
scaled <- T

## Create summary vectors only from cell classes we're going to model with: cell classes with 1+ nuclei, present in most datasets, neu = exc/inh

cell_classes <- c("ASC", "MIC", "OG", "OPC", "EXC", "INH", "MGE", "CGE") 

# calc_class_summary_vecs(datinfo, data_type, expr_type, summary_type, cell_classes, neu_subtypes, scaled)
# calc_class_random_vecs(datinfo, data_type, expr_type, summary_type, cell_classes, neu_subtypes, n_resamples)

## Add datinfo paths:

datinfo <- read.csv("../13_CT_stats/datinfo_SCSN_meta_analysis_CT_stats.csv") %>% na_if("") 

file_path <- file.path(getwd(), Sys.glob(paste0("data/", data_type, "/", expr_type, "/*class_", summary_type, "*")))
file_path <- file_path[!grepl("unscaled", file_path)]
datasets <- sapply(strsplit(sapply(strsplit(file_path, "/"), "[", 10), "_cell"), "[", 1)
idx <- match(datasets, datinfo$Dataset)
datinfo$Cell_Class_Mean_Vecs <- NA
datinfo$Cell_Class_Mean_Vecs[idx] <- file_path

file_path <- file.path(getwd(), Sys.glob(paste0("data/", data_type, "/", expr_type, "/*class_random_", summary_type, "*")))
datasets <- sapply(strsplit(sapply(strsplit(file_path, "/"), "[", 10), "_cell"), "[", 1)
idx <- match(datasets, datinfo$Dataset)
datinfo$Cell_Class_Random_Mean_Vecs <- NA
datinfo$Cell_Class_Random_Mean_Vecs[idx] <- file_path

fwrite(datinfo, file="datinfo_SCSN_meta_analysis_summary_vecs.csv")

############################################# Calc. mean CELL CLASS summary vectors: UNSCALED #############################################

summary_type <- "mean"
scaled <- F

## Create summary vectors from all cell types with >1 nuclei:

cell_classes <- c("ASC", "MIC", "OG", "OPC", "EXC", "INH", "NEU", "END", "PER", "VSMC") 

# calc_CT_summary_vecs(datinfo, data_type, expr_type, summary_type, cell_classes, neu_subtypes, scaled)

## Add datinfo paths:

# datinfo <- read.csv("../17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 

file_path <- file.path(getwd(), Sys.glob(paste0("data/", data_type, "/", expr_type, "/*class_", summary_type, "*")))
file_path <- file_path[grepl("unscaled", file_path)]
datasets <- sapply(strsplit(sapply(strsplit(file_path, "/"), "[", 10), "_cell"), "[", 1)
idx <- match(datasets, datinfo$Dataset)
datinfo$Cell_Class_Mean_Vecs_Unscaled <- NA
datinfo$Cell_Class_Mean_Vecs_Unscaled[idx] <- file_path

# fwrite(datinfo, file="../17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv")

############################################# Calc. mean CT summary vectors: SCALED #############################################

summary_type <- "mean"
scaled <- T

## Create summary vectors only from cell classes we're going to model with: cell classes/cell types with 1+ nuclei, present in most datasets, neu = exc/inh
## Note: there may be fewer CTs reprsented in random vectors because they require a cell type to be from a cell class with >1 cell type

cell_classes <- c("ASC", "MIC", "OG", "OPC", "EXC", "INH") 

# calc_CT_summary_vecs(datinfo, data_type, expr_type, summary_type, cell_classes, neu_subtypes, scaled)
# calc_CT_random_vecs(datinfo, data_type, expr_type, summary_type, cell_classes, neu_subtypes, n_resamples)

## Add datinfo paths:

# datinfo <- read.csv("../17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 

file_path <- file.path(getwd(), Sys.glob(paste0("data/", data_type, "/", expr_type, "/*CT_", summary_type, "*")))
file_path <- file_path[!grepl("unscaled", file_path)]
datasets <- sapply(strsplit(sapply(strsplit(file_path, "/"), "[", 10), "_CT"), "[", 1)
idx <- match(datasets, datinfo$Dataset)
datinfo$CT_Mean_Vecs <- NA
datinfo$CT_Mean_Vecs[idx] <- file_path

file_path <- file.path(getwd(), Sys.glob(paste0("data/", data_type, "/", expr_type, "/*CT_random_", summary_type, "*")))
datasets <- sapply(strsplit(sapply(strsplit(file_path, "/"), "[", 10), "_CT"), "[", 1)
idx <- match(datasets, datinfo$Dataset)
datinfo$CT_Random_Mean_Vecs <- NA
datinfo$CT_Random_Mean_Vecs[idx] <- file_path

fwrite(datinfo, file="datinfo_SCSN_meta_analysis_summary_vecs.csv")

############################################# Calc. mean CT summary vectors: UNSCALED #############################################

summary_type <- "mean"
scaled <- F

## Create summary vectors from all cell classes with >1 nuclei:

cell_classes <- c("ASC", "MIC", "OG", "OPC", "EXC", "INH", "NEU", "END", "PER", "VSMC") 

# calc_class_summary_vecs(datinfo, data_type, expr_type, summary_type, cell_classes, neu_subtypes, scaled)

## Add datinfo paths:

# datinfo <- read.csv("../17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 

file_path <- file.path(getwd(), Sys.glob(paste0("data/", data_type, "/", expr_type, "/*CT_", summary_type, "*")))
file_path <- file_path[grepl("unscaled", file_path)]
datasets <- sapply(strsplit(sapply(strsplit(file_path, "/"), "[", 10), "_CT"), "[", 1)
idx <- match(datasets, datinfo$Dataset)
datinfo$CT_Mean_Vecs_Unscaled <- NA
datinfo$CT_Mean_Vecs_Unscaled[idx] <- file_path

# fwrite(datinfo, file="../17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv")

############################################# Calc. eigen CT summary vectors: SCALED #############################################

summary_type <- "eigen"
scaled <- T

## Create summary vectors from all cell types with >1 nuclei:

cell_classes <- c("ASC", "MIC", "OG", "OPC", "EXC", "INH", "NEU", "END", "PER", "VSMC") 

# calc_CT_summary_vecs(datinfo, data_type, expr_type, summary_type, cell_classes, neu_subtypes, scaled)

## Add datinfo paths:

# datinfo <- read.csv("../17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 

file_path <- file.path(getwd(), Sys.glob(paste0("data/", data_type, "/", expr_type, "/*CT_", summary_type, "*")))
datasets <- sapply(strsplit(sapply(strsplit(file_path, "/"), "[", 10), "_CT"), "[", 1)
idx <- match(datasets, datinfo$Dataset)
datinfo$CT_Eigen_Vecs <- NA
datinfo$CT_Eigen_Vecs[idx] <- file_path

# fwrite(datinfo, file="../17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv")
