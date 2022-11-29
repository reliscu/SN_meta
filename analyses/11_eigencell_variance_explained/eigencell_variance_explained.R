setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/11_eigencell_variance_explained")

source("eigencell_variance_explained_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

## Restrict to cortical datasets:

datinfo <- datinfo %>% na_if("") %>% dplyr::filter(grepl("C$", Region_Code)) %>% dplyr::filter(!is.na(Author_Barcode_Annotations))

#eigencell_VE(datinfo, expr_type, gene_list, cell_classes, na_cts)
