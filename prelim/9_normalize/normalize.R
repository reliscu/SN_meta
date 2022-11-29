setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/9_normalize")

library(dplyr)
library(data.table)
library(future.apply)

plan(multicore, workers=10)

source("/home/rebecca/code/misc/normalize_functions.R")

datinfo <- read.csv("../7_max_expr/datinfo_SCSN_meta_analysis_max_expr.csv") %>% na_if("")

## Restrict to intersection genes for normalization:

gene_list <- fread("../8_intersection_genes/intersection_protein_coding_genes_33_cortical_datasets_8596_genes.csv", data.table=F)[,1]

future_lapply(1:nrow(datinfo), FUN=function(i){ 
  
  expr <- fread(datinfo$Author_Counts_QC_PC_Genes[i], data.table=F) 
  expr <- expr[match(gene_list, expr$SYMBOL),]
  if(!identical(gene_list, expr$SYMBOL)){
    stop(paste("Not all genes are present in", datinfo$Dataset[i]))
  }
  expr_norm <- log2(normalize_fxn(expr[,-c(1)], scale_factor=1e4)+1)
  expr_out <- data.frame(SYMBOL=expr[,c(1)], expr_norm)
  
  file_path <- gsub("7_max_expr", "9_normalize", datinfo$Author_Counts_QC_PC_Genes[i])
  file_path <- gsub("[0-9]+_protein_coding", paste0(length(gene_list), "_intersection"), file_path)
  file_path <- gsub(".csv", "_lib_norm_log2.csv", file_path)
  
  fwrite(expr_out, file=file_path, col.names=T)
  
})

## Add datinfo paths:

datinfo <- read.csv("../7_max_expr/datinfo_SCSN_meta_analysis_max_expr.csv") %>% na_if("")

file_path <- gsub("7_max_expr", "9_normalize", datinfo$Author_Counts_QC_PC_Genes)
file_path <- gsub("[0-9]+_protein_coding", paste0(length(gene_list), "_intersection"), file_path)
file_path <- gsub(".csv", "_lib_norm_log2.csv", file_path)

# sum(sapply(file_path, file.exists)==F)

datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes <- file_path

fwrite(datinfo, file="datinfo_SCSN_meta_analysis_normalized.csv")
