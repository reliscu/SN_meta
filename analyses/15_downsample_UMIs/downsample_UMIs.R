setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs")

source("downsample_UMIs_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/15_CT_summary_vecs/datinfo_SCSN_meta_analysis_summary_vecs.csv") %>% na_if("") %>% dplyr::filter(!is.na(Author_Barcode_Annotations_Class))

data_type <- "author_data"
expr_type <- "QC_counts"
n_umis_list <- c(400, 500, 600, 700, 800, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 10e3, 50e3, 100e3, 500e3, 800e3)
cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")

datinfo <- datinfo %>% dplyr::arrange(desc(Author_Median_UMIs_PC)) %>% dplyr::slice(1:5)
