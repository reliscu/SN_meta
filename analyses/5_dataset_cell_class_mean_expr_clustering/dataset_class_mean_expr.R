setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/5_dataset_cell_class_mean_expr_clustering")

library(dplyr)
library(data.table)

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/17_fidelity/datinfo_SCSN_meta_analysis_fidelity.csv") %>% na_if("") %>% dplyr::filter(!is.na(Cell_Class_Mean_Vecs))

data_type <- "author_data"
expr_type <- "normalized_counts"
cell_classes <- c("ASC", "END", "EXC", "INH", "OG", "OPC", "MIC", "PER", "VSMC")

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/prelim_new/8_intersection_genes/intersection_protein_coding_genes_33_cortical_datasets_8596_genes.csv", data.table=F)[,1]

lapply(1:length(cell_classes), function(j){
  
  print(cell_classes[j])
  
  mean_expr_list <- lapply(1:nrow(datinfo), function(i){
    mean_expr <- fread(datinfo$Cell_Class_Mean_Vecs[i], data.table=F)
    mean_expr <- mean_expr[match(gene_list, mean_expr[,1]),]
    if(!identical(gene_list, mean_expr[,1])){
      stop("!identical(gene_list, mean_expr[,1])")
    }
    mean_expr <- data.frame(mean_expr[,is.element(colnames(mean_expr), cell_classes[j])])
    if(ncol(mean_expr)>0){
      colnames(mean_expr) <- datinfo$Dataset[i]
      return(mean_expr)
    }
  })
  mean_expr_list <- mean_expr_list[lengths(mean_expr_list)>0]
  mean_expr <- do.call(cbind, mean_expr_list)
  datinfo1 <- datinfo[is.element(datinfo$Dataset, colnames(mean_expr)),]
  mean_expr <- data.frame(SYMBOL=gene_list, mean_expr, check.names=F)
  
  fwrite(mean_expr, file=paste0("data/", data_type, "/", expr_type, "/dataset_", cell_classes[j], "_mean_expr_", n_distinct(datinfo1$Dataset), "_datasets_", n_distinct(datinfo1$Study), "_studies_", length(gene_list), "_genes.csv"))
  
})

