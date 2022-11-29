library(dplyr)
library(Matrix)
library(data.table)
library(future.apply)
library(flexiblas)

source("/home/rebecca/code/misc/normalize_functions.R")

options(future.globals.maxSize=Inf)
plan(multicore, workers=2)

eigencell_VE <- function(
  datinfo, 
  expr_list_paths,
  data_type,
  expr_type
){
  
  for(i in 1:length(expr_list_paths)){
    
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
    
    n_genes <- nrow(expr)
    
    rm(expr)
    
    ## Ensure cell order in expr_list matches cell metadata:
    
    match_cells <- unlist(lapply(expr_list, function(x) !identical(colnames(x), cellinfo$Cell_ID)))
    
    if(sum(match_cells)>0){
      stop("Cells in expr_list do not match cell metadata")
    }
    
    cts <- unique(cellinfo$Cell_Type)
    
    ve_per_downsamp <- lapply(1:length(expr_list), function(j){

      expr <- log2(normalize_fxn(as(expr_list[[j]], "matrix"), scale_factor=1e4)+1)

      cts_ve <- future_lapply(cts, FUN=function(ct){
        
        print(ct)
        
        ct_idx <- which(
          is.element(cellinfo$Cell_Type, ct)
        )
        
        if(length(ct_idx)>1){
          
          ct_pcs <- prcomp(expr[,ct_idx], center=T, scale=T, retx=F)
          ve <- ct_pcs$sdev^2/sum(ct_pcs$sdev^2)
          
          return(data.frame(Cell_Type=ct, PC1_VE=ve[1], No.Nuclei=length(ct_idx)))
          
        } ## if(length(ct_idx)>1){

      }) ## cts_ve <- future_lapply(
      
      cts_ve_df <- do.call(rbind, cts_ve)
      cts_ve_df$Dataset <- datinfo$Dataset[i]
      cts_ve_df$Threshold <- names(expr_list)[j]
        
      return(cts_ve_df)

    }) ## ve_per_downsamp <- lapply(
    
    ve_per_downsamp <- do.call(rbind, ve_per_downsamp)
    
    fwrite(ve_per_downsamp, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_eigencell_VE_UMIs_downsampled_", data_type, "_", expr_type, "_", n_genes, "_genes.csv"))
    
  } ## for(i in 1:nrow(datinfo)){
  
} ## CT_PC1_VE
