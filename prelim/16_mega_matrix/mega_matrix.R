setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/16_mega_matrix")

source("mega_matrix_fxn.R")

datinfo <- read.csv("../15_CT_summary_vecs/datinfo_SCSN_meta_analysis_summary_vecs.csv") %>% na_if("") %>% dplyr::filter(!is.na(Author_Barcode_Annotations_Class))

gene_list <- fread("../8_intersection_genes/intersection_protein_coding_genes_33_cortical_datasets_8596_genes.csv", data.table=F)[,1]

## Restrict to cortical datasets:

datinfo <- datinfo %>% na_if("") %>% dplyr::filter(grepl("C$", Region_Code)) %>% dplyr::filter(!is.na(Author_Barcode_Annotations))

data_type <- "author_data"
expr_type <- "normalized_counts"
cell_classes <- c("NEU", "ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")

#mega_matrix(datinfo, data_type, expr_type, cell_classes, gene_list)