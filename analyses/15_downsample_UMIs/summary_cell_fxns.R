library(dplyr)
library(data.table)
library(future.apply)
library(gmodels) # fast.prcomp()

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

source("/home/rebecca/code/misc/normalize_functions.R")

calc_class_summary_vecs <- function(datinfo, expr_paths, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "CGE", "MGE"), neu_subtypes, scaled=F){
  
  neu_subtypes$Label <- make.names(neu_subtypes$Label)
  idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
  neu_subtypes$Class_Level1[idx] <- NA
  
  for(i in 1:nrow(datinfo)){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr_list <- readRDS(expr_paths[i])
    
    ## Restrict to genes with nonzero expression in all datasets:
    
    gene_list <- lapply(expr_list, function(expr){
      return(expr[rowSums(expr[,-c(1)])>0,1])
    })
    gene_list <- Reduce(intersect, gene_list)
    
    summary_vec_list <- lapply(1:length(expr_list), function(n){
      
      expr <- expr_list[[n]]
      
      idx <- match(gene_list, expr[,c(1)])
      genes <- expr[idx, c(1)]
      expr <- expr[idx, -c(1)]
      
      cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
      cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
      cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
      
      if(!identical(colnames(expr), cellinfo$Cell_ID)){
        stop("Order of cells in expression matrix does not match metadata!")
      }
    
      ## Add subtype info:
      
      temp <- make.names(paste0(datinfo$Dataset[i], "_", gsub(" ", "_", cellinfo$Cell_Type)))
      neu_subtypes1 <- neu_subtypes[match(temp, make.names(neu_subtypes$Label)),]
      idx <- which(!is.na(neu_subtypes1$Class_Level1))
      cellinfo$MO_Cell_Class[idx] <- neu_subtypes1$Class_Level1[idx]
      
      ## Restrict to cell classes present in dataset (and with >1 nuclei):
      
      idx <- which(is.element(cellinfo$MO_Cell_Class, cell_classes))
      n_class_nuclei <- table(cellinfo$MO_Cell_Class[idx])
      classes <- names(n_class_nuclei)[n_class_nuclei>1]
      classes <- classes[is.element(classes, cell_classes)]
      
      ## Restrict cells being sampled from:
      
      idx <- which(is.element(cellinfo$MO_Cell_Class, classes))
      cellinfo <- cellinfo[idx,]
      expr <- expr[,idx]
      
      ## Include broader cell classes if cells have more specific annotations:
      
      if(sum(is.element(c("EXC", "INH"), classes))>0){
        classes <- unique(c(classes, "NEU"))
        if(sum(is.element(c("CGE", "MGE"), classes))>0){
          classes <- unique(c(classes, "INH"))
        }
      }
      
      expr <- log2(normalize_fxn(expr, scale_factor=1e4)+1)
      
      ## Calculate cell type summary vectors:
      
      if(scaled){
        expr <- scale(expr)
      }  
      
      summary_list <- future_lapply(1:length(classes), FUN=function(j){
        
        if(is.element("NEU", classes[j])){
          idx <- which(is.element(cellinfo$MO_Cell_Class, c("NEU", "EXC", "INH", "CGE", "MGE")))
        } else if(is.element("INH", classes[j])) {
          idx <- which(is.element(cellinfo$MO_Cell_Class, c("INH", "CGE", "MGE")))
        } else {
          idx <- which(is.element(cellinfo$MO_Cell_Class, classes[j]))
        }

        if(summary_type=="mean"){
          
          summary_vec <- rowMeans(expr[,idx])
          
        } else {
          
          summary_vec <- fast.prcomp(expr[,idx], center=F)$x[,1]
          
        }
        
        return(summary_vec)
        
      }) ## summary_vec_list <- lapply(
      names(summary_list) <- classes
      
      # df <- data.frame(SYMBOL=genes, do.call(cbind, summary_list))
      # colnames(df) <- c("SYMBOL", classes)
      
      return(data.frame(SYMBOL=genes, do.call(cbind, summary_list)))
      
    }) ## lapply(1:length(expr_list)
    
    names(summary_vec_list) <- names(expr_list)
    
    file_path <- paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_UMIs_downsampled_cell_class_", summary_type, "cells_", data_type, "_", expr_type, "_", length(gene_list), "_genes.RDS")

    if(scaled==F){
      file_path <- gsub(paste0(expr_type, "_"), paste0("unscaled_", expr_type, "_"), file_path)
    }
    
    saveRDS(summary_vec_list, file=file_path)
    
  } ## for(i in 1:nrow(datinfo)){
  
} 

calc_class_random_vecs <- function(datinfo, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "CGE", "MGE"), neu_subtypes, n_resamples=10){
  
  neu_subtypes$Label <- make.names(neu_subtypes$Label)
  idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
  neu_subtypes$Class_Level1[idx] <- NA
  
  for(i in 1:nrow(datinfo)){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr_list <- readRDS(expr_paths[i])
    
    ## Restrict to genes with nonzero expression in all datasets:
    
    gene_list <- lapply(expr_list, function(expr){
      return(expr[rowSums(expr[,-c(1)])>0,1])
    })
    gene_list <- Reduce(intersect, gene_list)
    
    random_vec_list <- lapply(1:length(expr_list), function(n){
      
      expr <- expr_list[[n]]
      genes <- expr[,c(1)]
      expr <- expr[,-c(1)]
      
      idx <- match(gene_list, expr[,c(1)])
      genes <- genes[idx]
      expr <- expr[idx,]
      
      cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
      cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
      cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
      
      if(!identical(colnames(expr), cellinfo$Cell_ID)){
        stop("Order of cells in expression matrix does not match metadata!")
      }
      
      ## Add subtype info:
      
      temp <- make.names(paste0(datinfo$Dataset[i], "_", gsub(" ", "_", cellinfo$Cell_Type)))
      neu_subtypes1 <- neu_subtypes[match(temp, make.names(neu_subtypes$Label)),]
      idx <- which(!is.na(neu_subtypes1$Class_Level1))
      cellinfo$MO_Cell_Class[idx] <- neu_subtypes1$Class_Level1[idx]
      
      ## Restrict to cell classes present in dataset (and with >1 nuclei):
      
      idx <- which(is.element(cellinfo$MO_Cell_Class, cell_classes))
      n_class_nuclei <- table(cellinfo$MO_Cell_Class[idx])
      classes <- names(n_class_nuclei)[n_class_nuclei>1]
      classes <- classes[is.element(classes, cell_classes)]
      
      ## Restrict cells being sampled from:
      
      idx <- which(is.element(cellinfo$MO_Cell_Class, classes))
      cellinfo <- cellinfo[idx,]
      expr <- expr[,idx]
      
      ## Include broader cell classes if cells have more specific annotations:
      
      if(sum(is.element(c("EXC", "INH"), classes))>0){
        classes <- unique(c(classes, "NEU"))
        if(sum(is.element(c("CGE", "MGE"), classes))>0){
          classes <- unique(c(classes, "INH"))
        }
      }
      
      expr <- log2(normalize_fxn(expr, scale_factor=1e4)+1)
      expr <- scale(expr)
      
      ## Calculate cell type summary vectors:

      rand_list <- future_lapply(1:length(classes), FUN=function(j){
        
        if(is.element("NEU", classes[j])){
          idx <- which(is.element(cellinfo$MO_Cell_Class, c("NEU", "EXC", "INH", "CGE", "MGE")))
        } else if(is.element("INH", classes[j])) {
          idx <- which(is.element(cellinfo$MO_Cell_Class, c("INH", "CGE", "MGE")))
        } else {
          idx <- which(is.element(cellinfo$MO_Cell_Class, classes[j]))
        }
        
        random_class_list <- lapply(1:n_resamples, function(k){
          
          ## Sample equal # cells per major cell class:
          
          cellinfo1 <- cellinfo
          neu_idx <- which(is.element(cellinfo$MO_Cell_Class, c("EXC", "INH", "CGE", "MGE")))
          cellinfo1$MO_Cell_Class[neu_idx] <- "NEU"
          
          n_samples <- ceiling(length(idx)/n_distinct(cellinfo1$MO_Cell_Class))
          
          temp <- cellinfo1 %>%
            dplyr::mutate(idx=1:n()) %>%
            dplyr::group_by(MO_Cell_Class) %>%
            dplyr::slice_sample(n=n_samples, replace=T)
          
          if(summary_type=="mean"){
            
            summary_vec <- rowMeans(expr[,temp$idx])
            
          } else {
            
            summary_vec <- fast.prcomp(expr[,temp$idx], center=F)$x[,1]
            
          }
          
          return(summary_vec)
          
        })
        
        rand_vecs <- do.call(cbind, random_class_list)
        colnames(rand_vecs) <- rep(classes[j], ncol(rand_vecs))
        
        return(data.frame(SYMBOL=genes, rand_vecs, check.names=F))
        
      }, future.seed=T) ## rand_list <- future_lapply(
      
      names(rand_list) <- classes
      
    })
    
    names(random_vec_list) <- names(expr_list)
    
    saveRDS(summary_vec_list, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_UMIs_downsampled_cell_class_random_", summary_type, "cells_", data_type, "_", expr_type, "_", nrow(expr), "_genes.RDS"))
    
  } ## for(i in 1:nrow(datinfo)){
  
} 
