setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/16_cell_class_composition")

source("cell_class_comp_plot_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"

cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "NEU", "OG", "OPC", "PER", "VSMC")
na_cts <- c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class")

sim_type <- "pearson"
clust_method <- "average"

df <- read.csv(list.files(path=file.path("data", data_type), full.names=T))
