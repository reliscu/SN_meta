setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/16_cell_class_composition/sample_class_composition")

source("sample_cell_comp_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "QC_counts"

## For now restrict to cell classes present in most datasests:

cell_classes <- c("INH", "EXC", "ASC", "MIC", "OG", "OPC") # "END", PER", "VSMC

n_nuclei <- 1000
n_resamples <- 1

excl <- c(
  "Grubman_2019_10x_Chromium_V3_Entorhinal_cortex",
  "Luo_2019_10x_Chromium_V3_BA44-45",
  "Luo_2019_snmC2T-seq_BA10"
)

datinfo <- datinfo %>% na_if("") %>%
  dplyr::filter(grepl("C$", Region_Code)) %>%
  dplyr::filter(!is.na(Author_Barcode_Annotations)) %>% 
  dplyr::filter(!is.element(Dataset, excl))

cell_props <- read.csv("../data/author_data/cell_class_composition_30_datasets.csv")

## Take average proportion of each cell class from UNBIASED datasets:

idx <- which(
  datinfo$Unbiased_Sampling=="Y"
)

cell_props <- cell_props[is.element(cell_props$Dataset, datinfo$Dataset[idx]),]

prop_df <- cell_props %>%
  dplyr::group_by(Cell_Class) %>%
  dplyr::summarise(
    Prop=mean(Proportion)
  ) 

prop_df <- prop_df[match(cell_classes, prop_df$Cell_Class),]

props <- prop_df$Prop

## Should sum to 1:

props <- props/sum(props)

#sample_cell_comp(datinfo, cell_classes, props, n_nuclei)