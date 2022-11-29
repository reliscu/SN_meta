setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/1_summary_figures")

source("summary_plots_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "NEU", "OG", "OPC", "PER", "VSMC")
na_cts <- c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class")

## Restrict to cortical datasets:

datinfo <- datinfo %>% dplyr::filter(grepl("C$", Region_Code))

#trend_plots(datinfo)
