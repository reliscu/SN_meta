setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/7_max_expr")

library(dplyr)
library(Matrix)
library(data.table)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=13)

datinfo <- read.csv("../5_QC/datinfo_SCSN_meta_analysis_QC.csv") %>% na_if("")

gene_list <- fread("/home/rebecca/SCSN_meta_analysis/prelim_new/6_gene_biotype/protein_coding_20051_union_genes.csv", data.table=F)[,1]

data_type <- "author_data"
expr_type <- "QC_counts"

future_lapply(1:nrow(datinfo), FUN=function(i){ 
  
  cat("\n")
  print(datinfo$Dataset[i])
  
  expr_path <- datinfo$Author_Counts_QC[i]
  
  if(grepl("mtx", expr_path)){
    
    expr <- readMM(expr_path)
    
  } else {
    
    expr <- fread(expr_path, data.table=F)
    
  }
  
  expr <- as(expr, "matrix")
  
  genes <- fread(datinfo$Author_Genes_Mapped[i], data.table=F)
  barcodes <- fread(datinfo$Author_Barcodes_QC[i], header=F, data.table=F)
  
  ## Restrict to PC genes:
  
  idx <- which(is.element(genes$SYMBOL, gene_list))
  genes <- genes[idx,]
  expr <- expr[idx,]
  
  ## Remove zero variance features:
  
  idx <- which(rowSums(expr)>0 )
  genes <- genes[idx,]
  expr <- expr[idx,]
  
  if(sum(duplicated(genes$SYMBOL))>0){
    
    ## If 2+ aliases map to an official symbol, keep the alias that matches the official symbol:
    
    dup_ids <- genes$SYMBOL[duplicated(genes$SYMBOL)]
    
    rm_idx <- c()
    for(j in 1:length(dup_ids)){
      temp <- genes[is.element(genes$SYMBOL, dup_ids[j]),]
      if(sum(temp$UNIQUE.ID==temp$SYMBOL)>0){
        nonmatching <- temp$UNIQUE.ID[temp$UNIQUE.ID!=temp$SYMBOL]
        rm_idx <- c(rm_idx, which(is.element(genes$UNIQUE.ID, nonmatching)))
      }
    }
    
    if(!is.null(rm_idx)){
      
      genes <- genes[-rm_idx,]
      expr <- expr[-rm_idx,]
      
    }
    
    ## For remaining duplicates, choose gene with highest mean expression
    
    if(sum(duplicated(genes$SYMBOL))>0){
      
      max_expr <- genes %>%
        dplyr::mutate(
          idx=1:n(), Mean_Expr=rowMeans(expr)
        ) %>%
        dplyr::group_by(SYMBOL) %>%
        dplyr::slice_max(
          order_by=Mean_Expr, with_ties=F
        )
      
      idx <- sort(max_expr$idx)
      expr <- expr[idx,]
      genes <- genes[idx,]
      
    } 
    
  } ## if(sum(duplicated(genes$SYMBOL))>0){
  
  colnames(expr) <- barcodes[,1]
  expr_out <- data.frame(SYMBOL=genes$SYMBOL, expr, check.names=F)
  
  fwrite(expr_out, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_expression_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr_out), "_protein_coding_genes.csv"), col.names=T)
  
}) ## future_lapply( 

## Add datinfo paths:

datinfo <- read.csv("../5_QC/datinfo_SCSN_meta_analysis_QC.csv") %>% na_if("")

file_path <- list.files(path=file.path("data", data_type, expr_type))

datasets <- sapply(strsplit(file_path, "_expression"), "[", 1)

file_path <- file_path[match(datinfo$Dataset, datasets)]

sum(sapply(datinfo$Author_Counts_QC_PC_Genes, file.exists)==F)

datinfo$Author_Counts_QC_PC_Genes <- file.path(getwd(), "data", data_type, expr_type, file_path)

fwrite(datinfo, file="datinfo_SCSN_meta_analysis_max_expr.csv")
