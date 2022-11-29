setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/eigencell_similarity")

source("eigencell_calc_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")
na_cts <- c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class")

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/GRCh38_gencode_primary_SYMBOL_protein_coding.csv", data.table=F)[,1]

expr_list_paths <- list.files("../data/author_data/QC_counts", pattern="downsampled", full.names=T)

datasets <- sapply(strsplit(sapply(strsplit(expr_list_paths, "/"), "[", 5), "_UMIs"), "[", 1)

datinfo <- datinfo[match(datasets, datinfo$Dataset),]
