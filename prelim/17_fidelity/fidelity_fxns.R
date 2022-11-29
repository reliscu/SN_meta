library(dplyr)
library(Matrix)
library(data.table)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

source("/home/rebecca/SCSN_meta_analysis/code/fidelity_fxn.R")

cell_class_fidelity <- function(datinfo, data_type, expr_type, cell_classes=c("ASC", "END", "NEU", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"), prop_scaled=F){
  
  for(i in 1:nrow(datinfo)){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr <- fread(datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
    genes <- expr[,c(1)]
    expr <- expr[,-c(1)]
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
    cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
    
    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("Order of cells in expression matrix does not match metadata!")
    }
    
    ## Restrict to cell types in cell classes of interest:
    
    idx <- which(is.element(cellinfo$MO_Cell_Class, cell_classes))
    cellinfo <- cellinfo[idx,]
    expr <- expr[,idx]
    
    classes <- unique(cellinfo$MO_Cell_Class)
    
    if(sum(is.element(c("EXC", "INH"), classes))){
      classes <- unique(c("NEU", classes))
    }
    
    fid_list <- future_lapply(1:length(classes), FUN=function(j){
      
      idx <- which(is.element(cellinfo$MO_Cell_Class, classes[j]))
      if(classes[j]=="NEU"){
        idx <- which(is.element(cellinfo$MO_Cell_Class, c("NEU", "EXC", "INH")))
      }
      
      return(fidelity(expr, idx, prop_scaled))
      
    }) 
    
    df <- do.call(cbind.data.frame, fid_list)
    colnames(df) <- classes
    df <- data.frame(SYMBOL=genes, df)
    
    file_path <- paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_cell_class_fidelity_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv")
    
    if(prop_scaled){
      file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
    } 
    
    fwrite(df, file=file_path)
    
  } ## for(i in 1:nrow(datinfo)){
  
} 

CT_fidelity <- function(datinfo, data_type, expr_type, cell_classes=c("ASC", "END", "NEU", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"), prop_scaled=F){
  
  for(i in 1:nrow(datinfo)){ 
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr <- fread(datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
    genes <- expr[,c(1)]
    expr <- expr[,-c(1)]
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
    cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
    
    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("Order of cells in expression matrix does not match metadata!")
    }
    
    ## Restrict to cell types in cell classes of interest / with >1 nuclei:
    
    cellinfo <- cellinfo %>%
      dplyr::mutate(idx=1:n()) %>%
      dplyr::filter(
        is.element(cellinfo$MO_Cell_Class, cell_classes)
      ) %>%
      dplyr::group_by(Cell_Type) %>%
      dplyr::filter(n()>1) 

    expr <- expr[,cellinfo$idx]

    cts <- unique(cellinfo$Cell_Type)

    fid_list <- future_lapply(1:length(cts), FUN=function(j){
      idx <- which(is.element(cellinfo$Cell_Type, cts[j]))
      return(fidelity(expr, idx, prop_scaled))
    }) 
    
    df <- do.call(cbind.data.frame, fid_list)
    colnames(df) <- cts
    df <- data.frame(SYMBOL=genes, df, check.names=F)
    
    file_path <- paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_CT_fidelity_", data_type, "_", expr_type, "_", ncol(df)-1, "_CTs_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv")
    
    if(prop_scaled){
      file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
    } 
    
    fwrite(df, file=file_path)
    
  } ## for(i in 1:nrow(datinfo)){
  
}