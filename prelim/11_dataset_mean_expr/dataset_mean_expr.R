setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/11_dataset_mean_expr")

library(dplyr)
library(data.table)
library(future.apply)

plan(multicore, workers=10)

datinfo <- read.csv("../10_dataset_stats/datinfo_SCSN_meta_analysis_dataset_stats.csv") %>% na_if("")

data_type <- "author_data"
expr_type <- "normalized_counts"

future_lapply(1:nrow(datinfo), FUN=function(i){ 
  
  expr <- fread(datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
  genes <- expr[,c(1)]
  mean_expr <- rowMeans(expr)
  expr_out <- data.frame(SYMBOL=genes, Mean_Expr=mean_expr)
  
  file_path <- gsub("9_normalize", "11_dataset_mean_expr", datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes[i])
  file_path <- gsub("QC_counts", expr_type,  gsub(".csv", "_mean_expr.csv", file_path))
  
  fwrite(expr_out, file=file_path, col.names=T)
  
})

## Add datinfo paths:

datinfo <- read.csv("../10_dataset_stats/datinfo_SCSN_meta_analysis_dataset_stats.csv") %>% na_if("")

file_path <- gsub("9_normalize", "11_dataset_mean_expr", datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes)
file_path <- gsub("QC_counts", expr_type,  gsub(".csv", "_mean_expr.csv", file_path))

sum(sapply(file_path, file.exists)==F)

datinfo$Author_Mean_Normalized_Counts_QC_Intersection_PC_Genes <- file_path

fwrite(datinfo, file="datinfo_SCSN_meta_analysis_mean_expr.csv")