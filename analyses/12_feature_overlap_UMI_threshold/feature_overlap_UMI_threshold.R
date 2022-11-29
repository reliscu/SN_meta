setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/12_feature_overlap_UMI_threshold")

library(data.table)

source("feature_overlap_UMI_threshold_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/15_CT_summary_vecs/datinfo_SCSN_meta_analysis_summary_vecs.csv") %>% na_if("") %>% dplyr::filter(!is.na(Author_Barcode_Annotations_Class))

data_type <- "author_data"
expr_type <- "QC_counts"

umi_cut_list <- c(0, .0005, .01, .1, 1, 10, 500) ## Based on SANITY paper analysis: https://www.nature.com/articles/s41587-021-00875-x/figures/3

jaccard <- F

# feature_overlap_heatmap(datinfo, data_type, expr_type, umi_cut_list, jaccard, pc_genes)
# 
# feature_overlap_boxplot(datinfo, data_type, expr_type, umi_cut_list, jaccard, pc_genes)
