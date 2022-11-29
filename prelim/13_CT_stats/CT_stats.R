setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/13_CT_stats")

library(dplyr)
library(data.table)
library(future.apply)

plan(multicore, workers=10)

datinfo <- read.csv("../12_cell_class_annotations/datinfo_SCSN_meta_analysis_cell_class.csv") %>% na_if("")

## Use all PC genes to calculate stats (versus intersection genes)

data_type <- "author_data"
expr_type <- "QC_counts"

future_lapply(1:nrow(datinfo), FUN=function(i){
  
  if(!is.na(datinfo$Author_Barcode_Annotations_Class[i])){
    
    cat("\n")
    print(datinfo$Dataset[i])
   
    expr <- fread(datinfo$Author_Counts_QC_PC_Genes[i], data.table=F)
    expr <- expr[,-c(1)]
    cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F) %>% na_if("")
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]

    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("!identical(colnames(expr), cellinfo$Cell_ID)")
    }
    
    cts <- unique(cellinfo$Cell_Type) %>% na.omit()
    
    stats_list <- lapply(1:length(cts), function(j){
      
      idx <- which(is.element(cellinfo$Cell_Type, cts[j]))
      expr_ct <- as.matrix(expr[,idx])
      n_umis <- median(colSums(expr_ct))
      n_genes <- median(apply(expr_ct, 2, function(x) sum(x>0)))
      
      return(
        data.frame(Cell_Type=cts[j], No.Nuclei=length(idx), Median_UMIs=n_umis, Median_Unique_Genes=n_genes)
      )
      
    }) ##  for(j in 1:length(cell_class)){
    
    df <- do.call(rbind, stats_list)
    
    ## Add cell class annotations:
    
    cellinfo <- cellinfo %>% 
      dplyr::arrange(Cell_Type) %>% 
      dplyr::filter(!duplicated(Cell_Type))
    
    idx <- match(df$Cell_Type, cellinfo$Cell_Type)
    df <- data.frame(df, Cell=Class=cellinfo$MO_Cell_Class[idx])

    file_path <- gsub("7_max_expr", "13_CT_stats", datinfo$Author_Counts_QC_PC_Genes[i])
    file_path <- gsub(".csv", "_CT_stats.csv", file_path)
    
    fwrite(df, file=file_path)
    
  } ## if(!is.na(datinfo$Author_Barcode_Annotations_Class[i])){
  
}) ## future_lapply(

## Add datinfo paths:

datinfo <- read.csv("../12_cell_class_annotations/datinfo_SCSN_meta_analysis_cell_class.csv") %>% na_if("")

file_path <- gsub("7_max_expr", "13_CT_stats", datinfo$Author_Counts_QC_PC_Genes)
file_path <- gsub(".csv", "_CT_stats.csv", file_path)

idx <- sapply(file_path, file.exists)==F

sum(idx)

datinfo$Author_Barcode_Annotations_Class[idx]

datinfo$Author_Counts_QC_PC_Genes_CT_Stats <- file_path

datinfo$Author_Counts_QC_PC_Genes_CT_Stats[idx] <- NA

fwrite(datinfo, file="datinfo_SCSN_meta_analysis_CT_stats.csv")
