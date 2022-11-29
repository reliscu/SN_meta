library(data.table)
library(dplyr)
library(gmodels) ## fast.prcomp()
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=5)

calc_class_summary_vecs <- function(expr, cellinfo, neu_subtypes, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes, scaled=F){
  
  if(!identical(colnames(expr), cellinfo$Cell_ID)){
    stop("Order of cells in expression matrix does not match metadata!")
  }
  
  ## Add subtype info:
  
  neu_subtypes$Label <- make.names(neu_subtypes$Label)
  idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
  neu_subtypes$Class_Level1[idx] <- NA
  temp <- make.names(paste0(cellinfo$Dataset, "_", gsub(" ", "_", cellinfo$Cell_Type)))
  neu_subtypes <- neu_subtypes[match(temp, neu_subtypes$Label),]
  cellinfo$MO_Cell_Class[!is.na(neu_subtypes$Class_Level1)] <- neu_subtypes$Class_Level1[!is.na(neu_subtypes$Class_Level1)]
  
  ## Calculate cell class summary vectors:
  
  if(scaled){
    expr <- scale(expr)
  }
  
  summary_list <- future_lapply(1:length(cell_classes), FUN=function(j){

    cat("\n")
    print(cell_classes[j])
    
    idx <- which(is.element(cellinfo$MO_Cell_Class, cell_classes[j]))
    
    if(is.element("NEU", cell_classes[j])){
      idx <- which(is.element(cellinfo$MO_Cell_Class, c("NEU", "EXC", "INH", "CGE", "MGE")))
    } else if(is.element("INH", cell_classes[j])) {
      idx <- which(is.element(cellinfo$MO_Cell_Class, c("INH", "CGE", "MGE")))
    }

    if(summary_type=="mean"){
      
      summary_vec <- rowMeans(expr[,idx])
      
    } else {
      
      summary_vec <- fast.prcomp(expr[,idx], center=F)$x[,1]
      
    }
    
    return(summary_vec)
    
  }) ## summary_vec_list <- lapply(
  
  df <- do.call(cbind, summary_list)
  colnames(df) <- cell_classes
  df <- data.frame(SYMBOL=rownames(expr), df)
  
  file_path <- paste0("data/", data_type, "/", expr_type, "/cell_class_", summary_type, "cells_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv")
  
  if(scaled==F){
    file_path <- gsub(paste0(expr_type, "_"), paste0("unscaled_", expr_type, "_"), file_path)
  }
  
  fwrite(df, file=file_path)
  
}

calc_class_random_vecs <- function(expr, cellinfo, neu_subtypes, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes, n_resamples=10){
  
  if(!identical(colnames(expr), cellinfo$Cell_ID)){
    stop("Order of cells in expression matrix does not match metadata!")
  }
  
  ## Add subtype info:
  
  neu_subtypes$Label <- make.names(neu_subtypes$Label)
  idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
  neu_subtypes$Class_Level1[idx] <- NA
  temp <- make.names(paste0(cellinfo$Dataset, "_", gsub(" ", "_", cellinfo$Cell_Type)))
  neu_subtypes <- neu_subtypes[match(temp, neu_subtypes$Label),]
  cellinfo$MO_Cell_Class[!is.na(neu_subtypes$Class_Level1)] <- neu_subtypes$Class_Level1[!is.na(neu_subtypes$Class_Level1)]
  
  ## Calculate random cell summary vectors (same # cells as each cell class):
  
  expr <- scale(expr)
  
  random_list <- future_lapply(1:length(cell_classes), FUN=function(j){
    
    cat("\n")
    print(cell_classes[j])
    
    idx <- which(is.element(cellinfo$MO_Cell_Class, cell_classes[j]))
    
    if(is.element("NEU", cell_classes[j])){
      idx <- which(is.element(cellinfo$MO_Cell_Class, c("NEU", "EXC", "INH", "CGE", "MGE")))
    } else if(is.element("INH", cell_classes[j])) {
      idx <- which(is.element(cellinfo$MO_Cell_Class, c("INH", "CGE", "MGE")))
    }
    
    rand_list <- lapply(1:n_resamples, function(n){

      ## Sample equal # cells per major cell class:
      
      neu_idx <- which(is.element(cellinfo$MO_Cell_Class, c("EXC", "INH", "CGE", "MGE")))
      cellinfo1 <- cellinfo
      cellinfo1$MO_Cell_Class[neu_idx] <- "NEU"
      
      # No. samples per cell class:
      
      n_samples <- ceiling(length(idx)/n_distinct(cellinfo1$MO_Cell_Class))
      
      temp <- cellinfo1 %>%
        dplyr::mutate(idx=1:n()) %>%
        dplyr::group_by(MO_Cell_Class) %>%
        dplyr::slice_sample(n=n_samples, replace=T) %>%
        arrange(idx)

      if(summary_type=="mean"){
        
        summary_vec <- rowMeans(expr[,temp$idx])
        
      } else {
        
        summary_vec <- fast.prcomp(expr[,temp$idx], center=F)$x[,1]
        
      }
      
      return(summary_vec)
      
    }) ##  rand_list <- lapply(
    
    df <- do.call(cbind, rand_list)
    colnames(df) <- paste0(cell_classes[j], "_Random", 1:n_resamples)
    return(data.frame(SYMBOL=rownames(expr), df))
    
  }, future.seed=T) ## summary_vec_list <- lapply(
  
  names(random_list) <- cell_classes
  
  saveRDS(random_list, file=paste0("data/", data_type, "/", expr_type, "/random_cell_class_", summary_type, "cells_", data_type, "_", expr_type, "_", n_resamples, "_resamples_", ncol(expr), "_nuclei_", nrow(expr), "_genes.RDS"))
  
}

calc_class_summary_vecs_leftout <- function(expr, cellinfo, neu_subtypes, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes, scaled=F){
  
  if(!identical(colnames(expr), cellinfo$Cell_ID)){
    stop("Order of cells in expression matrix does not match metadata!")
  }
  
  ## Add subtype info:
  
  neu_subtypes$Label <- make.names(neu_subtypes$Label)
  idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
  neu_subtypes$Class_Level1[idx] <- NA
  temp <- make.names(paste0(cellinfo$Dataset, "_", gsub(" ", "_", cellinfo$Cell_Type)))
  neu_subtypes <- neu_subtypes[match(temp, neu_subtypes$Label),]
  cellinfo$MO_Cell_Class[!is.na(neu_subtypes$Class_Level1)] <- neu_subtypes$Class_Level1[!is.na(neu_subtypes$Class_Level1)]
  
  ## Calculate cell class summary vectors:
  
  if(scaled){
    expr <- scale(expr)
  }
  
  datasets <- unique(cellinfo$Dataset)
  
  leftout_list <- lapply(1:length(datasets), function(i){
    
    cat("\n")
    print(datasets[i])
    
    ## Remove cells from left out dataset:
    
    idx <- which(!is.element(cellinfo$Dataset, datasets[i]))
    cellinfo1 <- cellinfo[idx,]
    expr1 <- expr[,idx]
    
    summary_list <- future_lapply(1:length(cell_classes), FUN=function(j){
      
      cat("\n")
      print(cell_classes[j])
      
      idx <- which(is.element(cellinfo1$MO_Cell_Class, cell_classes[j]))
      
      if(is.element("NEU", cell_classes[j])){
        idx <- which(is.element(cellinfo1$MO_Cell_Class, c("NEU", "EXC", "INH", "CGE", "MGE")))
      } else if(is.element("INH", cell_classes[j])) {
        idx <- which(is.element(cellinfo1$MO_Cell_Class, c("INH", "CGE", "MGE")))
      }
      
      if(summary_type=="mean"){
        
        summary_vec <- rowMeans(expr1[,idx])
        
      } else {
        
        summary_vec <- fast.prcomp(expr1[,idx], center=F)$x[,1]
        
      }
      
      return(summary_vec)
      
    }) ## summary_vec_list <- lapply(
    
    names(summary_list) <- cell_classes
    
    return(data.frame(SYMBOL=rownames(expr), do.call(cbind, summary_list)))
    
  })
  
  names(leftout_list) <- datasets
  
  file_path <- paste0("data/", data_type, "/", expr_type, "/cell_class_", summary_type, "cells_left_out_dataset_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.RDS")
  
  if(scaled==F){
    file_path <- gsub(paste0(expr_type, "_"), paste0("unscaled_", expr_type, "_"), file_path)
  }
  
  saveRDS(leftout_list, file=file_path)
  
}
