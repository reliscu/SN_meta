setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/eigencell_variance_explained")

source("eigencell_violin_plot_fxn.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"

ve_paths <- list.files(path=file.path("data", data_type, expr_type), pattern="downsampled", full.names=T)

datasets <- sapply(strsplit(sapply(strsplit(ve_paths, "/"), "[", 4), "_CT_PC1"), "[", 1)

datinfo <- datinfo[match(datasets, datinfo$Dataset),]

ve <- lapply(1:length(ve_paths), function(i){
  temp <- fread(ve_paths[i], data.table=F)
  temp <-  temp[temp$Threshold!=1e6,]
  return(temp)
})

ve_df <- do.call(rbind, ve)

df <- ve_df

