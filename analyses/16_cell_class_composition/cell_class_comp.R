setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/16_cell_class_composition")

source("cell_class_comp_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"

cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "NEU", "OG", "OPC", "PER", "VSMC")
na_cts <- c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class")

## Restrict to cortical datasets:

datinfo <- datinfo %>% na_if("") %>% 
  dplyr::filter(grepl("C$", Region_Code)) %>% 
  dplyr::filter(!is.na(Author_Barcode_Annotations))

