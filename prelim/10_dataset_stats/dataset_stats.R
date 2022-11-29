setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/10_dataset_stats")

library(dplyr)
library(data.table)

datinfo <- read.csv("../9_normalize/datinfo_SCSN_meta_analysis_normalized.csv") %>% na_if("")

## Use all PC genes to calculate stats (versus intersection genes)

data_type <- "author_data"
expr_type <- "QC_counts"

datinfo$Author_No.Nuclei <- NA
datinfo$Author_No.Cell_Types <- NA
datinfo$Author_Median_UMIs_PC<- NA
datinfo$Author_Median_Unique_Genes_PC <- NA

for(i in 1:nrow(datinfo)){ 
  
  cat("\n")
  print(datinfo$Dataset[i])
  
  expr <- fread(datinfo$Author_Counts_QC_PC_Genes[i], data.table=F)
  expr <- expr[,-c(1)]
  
  if(!is.na(datinfo$Author_Barcode_Annotations[i])){
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F)
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("!identical(colnames(expr), cellinfo$Cell_ID)")
    }
    datinfo$Author_No.Cell_Types[i] <- n_distinct(cellinfo$Cell_Type)
    
  } 
  
  datinfo$Author_No.Nuclei[i] <- ncol(expr)
  datinfo$Author_Median_UMIs_PC[i] <- median(colSums(expr))
  datinfo$Author_Median_Unique_Genes_PC[i] <- median(apply(expr, 2, function(x) sum(x>0)))
  
} ## for(i in 1:nrow(datinfo)){

datinfo <- datinfo %>% dplyr::relocate(c(Author_No.Nuclei, Author_No.Cell_Types, Author_Median_UMIs_PC, Author_Median_Unique_Genes_PC, Plot_Label), .before="Directory")

fwrite(datinfo, file="datinfo_SCSN_meta_analysis_dataset_stats.csv")
