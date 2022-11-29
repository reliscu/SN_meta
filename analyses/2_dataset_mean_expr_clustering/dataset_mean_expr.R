setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/2_dataset_mean_expr_clustering")

library(dplyr)
library(data.table)

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") 
 
data_type <- "author_data"
expr_type <- "normalized_counts"

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/prelim_new/8_intersection_genes/intersection_protein_coding_genes_33_cortical_datasets_8596_genes.csv", data.table=F)[,1]

mean_expr_list <- lapply(1:nrow(datinfo), function(i){
  mean_expr <- fread(datinfo$Author_Mean_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
  mean_expr <- mean_expr[match(gene_list, mean_expr[,1]),]
  if(!identical(gene_list, mean_expr[,1])){
    stop("!identical(gene_list, mean_expr[,1])")
  }
  return(mean_expr[,-c(1)])
})

mean_expr <- do.call(cbind, mean_expr_list)
colnames(mean_expr) <- datinfo$Dataset
mean_expr <- data.frame(SYMBOL=gene_list, mean_expr, check.names=F)

fwrite(mean_expr, file=paste0("data/", data_type, "/", expr_type, "/dataset_mean_expr_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies_", length(gene_list), "_genes.csv"))
