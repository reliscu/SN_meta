library(dplyr)
library(Matrix)
library(data.table)
library(ggplot2)
library(future.apply)
library(tictoc)
library(propr) ## prop
library(WGCNA) ## cor (in parallel)
library(dismay) ## kendall_zi
library(pcaPP) ## cor.fk
library(proxyC) ## simil (in parallel)
library(RcppParallel) ## setThreadOptions

setThreadOptions(numThreads=10)
options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

source("/home/rebecca/code/misc/normalize_functions.R")

cell_nearest_neighbors <- function(
  datinfo,
  expr_list_paths, 
  sim_type=c("pearson", "cosine", "prop"),
  data_type=c("author_data"), 
  expr_type=c("normalized_counts")
){
  
  for(i in 1:length(expr_list_paths)){
    
    cat("\n")
    
    print(datinfo$Dataset[i])
    print(expr_list_paths[i])
    
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
    
    match_cells <- unlist(lapply(expr_list, function(x) !identical(colnames(x), cellinfo$Cell_ID)))
    
    if(sum(match_cells)>0){
      stop("Cells in expr_list do not match cell metadata")
    }
    
    ## Calculate number of nonzero genes per cell per threshold:
    
    n_nonzeros_list <- future_lapply(expr_list, FUN=function(expr){

      n_nonzeros <- apply(expr, 2, function(x) sum(x>0))

      return(n_nonzeros/nrow(expr))

    })

    ## Calculate cell-cell similarity for each threshold:
    
    if(sim_type!="prop"){
      
      sim_list <- lapply(expr_list, function(expr){
        
        expr <- log2(normalize_fxn(as(expr, "matrix"), scale_factor=1e4)+1)
        
        if(sim_type=="pearson"){
          
          sim_mat <- WGCNA::cor(expr)
          
        } else { ## cosine
          
          sim_mat <- proxyC::simil(expr, margin=2, method="cosine")
          
        }
  
        diag(sim_mat) <- 0
        
        return(sim_mat)
        
      })
      
    } else { ## if(sim_type=="pearson"){
      
      sim_list <- future_lapply(expr_list, FUN=function(expr){
        
        expr <- log2(normalize_fxn(as(expr, "matrix"), scale_factor=1e4)+1)
        sim_mat <- propr(expr, metric="rho", symmetrize=T, p=1)@matrix
        diag(sim_mat) <- 0
        
        return(sim_mat)
        
      }, future.seed=T)
      
    } ## if(sim_type=="pearson"){} else {

    thresholds <- names(expr_list)
    rm(expr_list)
    
    nn_sim_list <- future_lapply(1:length(sim_list), FUN=function(j){

      nn_idx <- apply(sim_list[[j]], 2, which.max)
      
      nn_df <- data.frame(
        Cell_ID=colnames(sim_list[[j]]),
        Neighbor_ID=colnames(sim_list[[j]])[nn_idx],
        Threshold=thresholds[j],
        Frac_Nonzero=n_nonzeros_list[[j]]
      )
      
      return(nn_df)
      
    }) ## nn_sim_list <- future_lapply(
    
    rm(sim_list)
    
    nn_df <- do.call(rbind, nn_sim_list)
    nn_df$No.Genes <- n_genes
    nn_df$Dataset <- sapply(strsplit(
      sapply(strsplit(expr_list_paths[i], "/"
    ), "[", 5), "_UMIs"), "[", 1)
    
    fwrite(nn_df, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_UMIs_downsampled_nearest_neighbor_", toupper(sim_type), "_similarity_", data_type, "_", expr_type, ".csv"))
    
  } ## for(i in 1:length(expr_list_paths)){
  
}
