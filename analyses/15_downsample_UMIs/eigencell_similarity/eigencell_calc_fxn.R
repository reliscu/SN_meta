library(dplyr)
library(Matrix)
library(data.table)
library(future.apply)

source("/home/rebecca/code/misc/normalize_functions.R")

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

calc_eigencell <- function(
  datinfo, 
  data_type=c("author_data"), 
  expr_type=c("counts", "QC_counts", "normalized_counts"), 
  cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"), 
  na_cts=c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class")
){
  
  for(i in 1:nrow(datinfo)){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr_list <- readRDS(expr_list_paths[i]) 
    
    ## Find intersection of nonzero genes:
    
    idx_list <- lapply(expr_list, function(expr) which(rowSums(expr)>0))
    idx <- Reduce(intersect, idx_list)
    expr_list <- lapply(expr_list, function(expr){
      expr[idx,]
    })
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F)

    ## Ensure cell order in expr_list matches cell metadata:
    
    expr <- expr_list[[1]]
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
    
    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("!identical(colnames(expr), cellinfo$Cell_ID)")
    }
    
    genes <- rownames(expr)
    
    rm(expr)
    
    ## Ensure cell order in expr_list matches cell metadata:
    
    match_cells <- unlist(lapply(expr_list, function(x) !identical(colnames(x), cellinfo$Cell_ID)))
    
    if(sum(match_cells)>0){
      stop("Cells in expr_list do not match cell metadata")
    }
    
    cts <- unique(
      cellinfo$Cell_Type[is.element(cellinfo$MO_Cell_Class, cell_classes)]
    )
    
    cts <- cts[!grepl(paste(na_cts, collapse="|"), cts)]
    
    df_list <- future_lapply(1:length(expr_list), FUN=function(i){
      
      expr <- log2(normalize_fxn(as(expr_list[[i]], "matrix"), scale_factor=1e4)+1)

      eigencell_list <- lapply(1:length(cts), function(j){
        
        ct_idx <- which(
          is.element(cellinfo$Cell_Type, cts[j])
        )
        
        if(length(ct_idx)>1){
          
          ct_expr <- expr[,ct_idx]
          ct_pcs <- prcomp(ct_expr, center=T, scale=T)
          
          eigencell <- data.frame(
            PC1=ct_pcs$x[,1], 
            No.Nuclei=length(ct_idx)
          )
          
          colnames(eigencell) <- paste0(
            gsub(" ", "_", cts[j]), "_", colnames(eigencell)
          )
          
          return(eigencell)
          
        } ## if(length(ct_idx)>1){
        
      }) ## eigencell_list <- future_lapply(
      
      null_entries <- unlist(lapply(eigencell_list, is.null))
      eigen_df <- do.call(cbind, eigencell_list[!null_entries])
      eigen_df$Dataset <- datinfo$Dataset[i]

      return(data.frame(SYMBOL=genes, eigen_df))
      
    }) ## df_list <- future_lapply(
    
    saveRDS(df_list, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_eigencells_UMIs_downsampled_", data_type, "_", expr_type, "_", length(genes), "_genes.RDS"))
    
  } ## for(i in 1:nrow(datinfo)){
  
} 
