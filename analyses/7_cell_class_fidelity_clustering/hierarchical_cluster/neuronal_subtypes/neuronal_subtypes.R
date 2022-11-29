setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/7_cell_class_fidelity_clustering/hierarchical_cluster/neuronal_subtypes")

source("neuronal_subtypes_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
which_genes <- "intersection"
prop_metric <- "rho"

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

neu_subtypes <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim/14_cell_type_annotations/neuronal_subtype_annotations.csv") %>% na_if("")

prop_scaled <- F

pattern <- "fidelity_author"
if(prop_scaled){
  pattern <- "prop"
}

class_fid_paths <- list.files(path=paste0("../../data/", data_type, "/", expr_type), pattern=pattern, full.names=T)
class_expr_paths <- list.files(path=paste0("../../../6_cell_class_mean_expr_clustering/data/", data_type, "/", expr_type), full.names=T)

## Restrict to neuronal cell types:

class_fid_paths <- class_fid_paths[grep("EXC|INH", class_fid_paths)]
class_expr_paths <- class_expr_paths[grep("EXC|INH", class_expr_paths)]

sim_type <- "pearson"
top_n <- NULL

## Get list of intersection genes for datasets with cell annotations:

for(i in 1:length(class_fid_paths)){

  cell_class <- unlist(
    strsplit(sapply(strsplit(class_fid_paths[i], "/"), function(x) x[length(x)]), "_")[[1]]
    )[1]

  print(cell_class)

  class_fid <- fread(class_fid_paths[i], data.table=F)
  
  for(clust_method in c("complete", "ward.D2", "average")){
    plot_dendro(class_fid, datinfo, cell_class, data_type, expr_type, prop_scaled, which_genes, gene_list, sim_type, prop_metric, clust_method, top_n, class_expr, neu_subtypes)
  }
  

}
