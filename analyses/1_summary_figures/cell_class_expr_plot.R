setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/1_summary_figures")

source("cell_class_expr_plot_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
which_genes <- "intersection"
top_n <- NULL

na_cts <- c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class")

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

## Restrict to cortical datasets:

datinfo <- datinfo %>% dplyr::filter(grepl("C$", Region_Code)) %>% dplyr::filter(!is.na(Author_Barcode_Annotations))
