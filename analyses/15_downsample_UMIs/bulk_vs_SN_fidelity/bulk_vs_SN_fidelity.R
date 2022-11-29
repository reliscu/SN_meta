setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/bulk_vs_SN_fidelity")

source("bulk_vs_SN_fidelity_plot.R")

source("/home/rebecca/code/misc/map_identifiers/map_identifiers_function.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
prop_scaled <- F

fid_pattern <- "fidelity_author"
if(prop_scaled){
  fid_pattern <- "prop"
}

sn_fid_paths <- list.files(path=paste0("data/", data_type, "/", expr_type), pattern=fid_pattern, full.names=T)

bulk_fid <- fread("/home/rebecca/kelley_hifi_genes_tableS3.csv", data.table=F)

bulk_kme <- fread("data/consensus_kME_kelley_neuron_enriched_9.5e-10_pval_2226_samples.csv", data.table=F)

## Prep bulk data:

colnames(bulk_fid) <- gsub("Neuron", "NEU", colnames(bulk_fid))
colnames(bulk_fid) <- gsub("Astrocyte", "ASC", colnames(bulk_fid))
colnames(bulk_fid) <- gsub("Endothelial", "END", colnames(bulk_fid))
colnames(bulk_fid) <- gsub("Microglia", "MIC", colnames(bulk_fid))
colnames(bulk_fid) <- gsub("Oligodendrocyte", "OG", colnames(bulk_fid))

## Map bulk gene symbols to latest identifiers:

bulk_kme <- mapAlias2Symbol(features=bulk_kme, unique_id_col=2, mapping_tables_dir="/home/rebecca/omicon/mapping_tables/", keep_all_mappings=F, fill_NAs=T)

bulk_fid <- mapAlias2Symbol(features=bulk_fid, unique_id_col=1, mapping_tables_dir="/home/rebecca/omicon/mapping_tables/", keep_all_mappings=F, fill_NAs=T)

## Combine SN datasets from cell class of interest:

datasets <- sapply(
  strsplit(sapply(strsplit(sn_fid_paths, "/"), "[", 4), "_cell_class"
  ), "[", 1)

cell_class <- "NEU"

class_fid_list <- lapply(1:length(datasets), function(i){
  sn_fid <- readRDS(sn_fid_paths[i])
  idx <- which(is.element(names(sn_fid), cell_class))
  class_fid <- sn_fid[[idx]]
  class_fid$Dataset <- datasets[i]
  return(class_fid)
})

class_sn_fid <- do.call(rbind, class_fid_list)

plot_bulk_vs_SN_fidelity(datinfo, class_sn_fid, bulk_fid, cell_class, data_type, expr_type, prop_scaled)
