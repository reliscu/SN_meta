setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/13_gene_nearest_neighbors")

source("gene_knn_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

pc_genes <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/protein_coding_19723_union_genes.csv", data.table=F)

cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")

data_type <- "author_data"
expr_type <- "normalized_counts"

n_resamples <- 100
n_neighbors <- c(1, 2, 5, 10, 20, 50)

## Restrict to cortical datasets:

datinfo <- datinfo %>% dplyr::filter(grepl("C$", Region_Code))

## Restrict to a subset of datasets for now:

datinfo <- datinfo %>% 
  dplyr::filter(
    is.element(Dataset, c(
      "Habib_2017_DroNc-seq_Prefrontal_cortex",
      "Velmeshev_2019_10x_Chromium_Prefrontal_cortex",
      "Tran_2020_10x_Chromium_V3_Dorsolateral_prefrontal_cortex",
      "Luo_2019_10x_Chromium_V3_BA46"
    ))
  )

# "Habib_2017_DroNc-seq_Prefrontal_cortex",
# "Lake_2018_snDrop-seq_BA6",
# "ABI_2019_SMART-Seq_v4_Anterior_cingulate_cortex"
# "Schirmer_2019_10x_Chromium_V2_Prefrontal_cortex",

#model_knn_vs_rand(datinfo, data_type, expr_type, n_neighbors, n_resamples, cell_classes, pc_genes)

