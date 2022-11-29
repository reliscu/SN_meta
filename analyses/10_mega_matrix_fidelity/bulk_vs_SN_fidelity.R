setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/10_bulk_vs_SN_fidelity")

source("bulk_vs_SN_fidelity_plot.R")

source("/home/rebecca/code/misc/map_identifiers/map_identifiers_function.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")
  
bulk_fid <- fread("/home/rebecca/kelley_hifi_genes_tableS3.csv", data.table=F)

data_type <- "author_data"
expr_type <- "normalized_counts"
cell_classes <- c("NEU", "ASC", "MIC", "OG", "OPC") # "PER", "VSMC", "EXC", "INH"

prop_scaled <- F

fid_pattern <- "fidelity_author"
if(prop_scaled){
  fid_pattern <- "prop"
}

sn_fid <- fread(list.files(path=paste0("data/", data_type, "/", expr_type), pattern=fid_pattern, full.names=T), data.table=F)

## Format bulk data:

colnames(bulk_fid) <- gsub("Neuron", "NEU", colnames(bulk_fid))
colnames(bulk_fid) <- gsub("Astrocyte", "ASC", colnames(bulk_fid))
colnames(bulk_fid) <- gsub("Endothelial", "END", colnames(bulk_fid))
colnames(bulk_fid) <- gsub("Microglia", "MIC", colnames(bulk_fid))
colnames(bulk_fid) <- gsub("Oligodendrocyte", "OG", colnames(bulk_fid))

## Map bulk gene symbols to latest identifiers:

mapping_tables_dir <- "/home/rebecca/omicon/mapping_tables/"

bulk_fid <- mapAlias2Symbol(features=bulk_fid, unique_id_col=1, mapping_tables_dir, keep_all_mappings=F, fill_NAs=T)

plot_bulk_vs_SN_fidelity(sn_fid, bulk_fid, data_type, expr_type, cell_classes, prop_scaled)
