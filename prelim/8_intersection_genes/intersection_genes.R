setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/8_intersection_genes")

library(dplyr)
library(data.table)

datinfo <- read.csv("../7_max_expr/datinfo_SCSN_meta_analysis_max_expr.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "QC_counts"

## Get intersection of nonzero protein coding genes in cortical datasets:

gene_list <- lapply(1:nrow(datinfo), function(i){
  
  print(datinfo$Dataset[i])
  expr <- fread(datinfo$Author_Counts_QC_PC_Genes[i], data.table=F) 
  return(expr$SYMBOL)
  
})

genes <- Reduce(intersect, gene_list)

fwrite(data.frame(SYMBOL=genes), file=paste0("intersection_protein_coding_genes_", length(gene_list), "_cortical_datasets_", length(genes), "_genes.csv"))
