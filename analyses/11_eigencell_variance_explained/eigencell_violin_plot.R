setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/11_eigencell_variance_explained")

source("eigencell_violin_plot_fxn.R")

data_type <- "author_data"
expr_type <- "normalized_counts"

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

## Restrict to cortical datasets:

datinfo <- datinfo %>% dplyr::filter(grepl("C$", Region_Code))

ve_df <- fread(list.files(path=paste0("data/", data_type, "/", expr_type), full.names=T), data.table=F)

ve_df <- ve_df[is.element(ve_df$Dataset, datinfo$Dataset),]

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

n_genes <- length(gene_list)

plot_eigencell_VE(ve_df, datinfo, data_type, expr_type, n_genes)
