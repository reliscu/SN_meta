setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/7_cell_class_fidelity_clustering/hierarchical_cluster")

source("hiearchical_cluster_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"
prop_metric <- "rho"
which_genes <- "intersection"

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/data/cortical_regions/intersection_protein_coding_genes_33_cortical_datasets_21_studies_8723_genes.csv", header=F, data.table=F)[,1]

prop_scaled <- F

pattern <- "fidelity_author"
if(prop_scaled){
  pattern <- "prop"
}

class_fid_paths <- list.files(path=paste0("../data/", data_type, "/", expr_type), pattern=pattern, full.names=T)
class_expr_paths <- list.files(path=paste0("../../6_cell_class_mean_expr_clustering/data/", data_type, "/", expr_type), full.names=T)

sim_type <- "pearson"
top_n <- NULL

for(i in 1:length(class_fid_paths)){
  
  cell_class <- strsplit(sapply(strsplit(class_fid_paths[i], "/"), function(x) x[length(x)]), "_")[[1]][1]
  
  cat("\n")
  
  print(class_fid_paths[i])
  print(class_expr_paths[i])
  
  class_fid <- fread(class_fid_paths[i], data.table=F)
  class_expr <- fread(class_expr_paths[i], data.table=F)

  for(clust_method in c("complete", "ward.D2", "average")){
    #for(sim_type in c("pearson", "prop")){
      plot_dendro(class_fid, datinfo, cell_class, data_type, expr_type, prop_scaled, which_genes, gene_list, sim_type, prop_metric, clust_method, top_n, class_expr)
   # }
  }

}
