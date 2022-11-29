library(dplyr)
library(data.table)
library(future.apply)
library(gmodels) # fast.prcomp()

options(future.globals.maxSize=Inf)
plan(multicore, workers=5)

calc_CT_summary_vecs <- function(datinfo, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC", "NEU"), scaled=F){
  
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
    
    ## Calculate cell type summary vectors:
    
    if(scaled){
      expr <- scale(expr)
    }
    
    ## Restrict to cell types with at least 2 cells:
    
    n_ct_nuclei <- table(cellinfo$Cell_Type)
    cts <- names(n_ct_nuclei)[n_ct_nuclei>1]
  
    summary_list <- future_lapply(1:length(cts), FUN=function(j){
      
      idx <- which(is.element(cellinfo$Cell_Type, cts[j]))
      
      if(summary_type=="mean"){
        
        summary_vec <- rowMeans(expr[,idx])
        
      } else {
        
        summary_vec <- fast.prcomp(expr[,idx], center=F)$x[,1]
        
      }
      
      return(summary_vec)
      
    }) 
    
    names(summary_list) <- cts
    
    df <- data.frame(SYMBOL=genes, do.call(cbind, summary_list))
    
    file_path <- paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_CT_", summary_type, "cells_", data_type, "_", expr_type, "_", length(cts), "_CTS_", nrow(expr), "_genes.csv")
    
    if(scaled==F){
      file_path <- gsub(paste0(expr_type, "_"), paste0("unscaled_", expr_type, "_"), file_path)
    }
    
    fwrite(df, file=file_path)
    
  }
  
} 

calc_CT_random_vecs <- function(datinfo, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC", "NEU"), n_resamples=10){
  
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
    
    ## Restrict to cell classes of interest:
    
    idx <- which(is.element(cellinfo$MO_Cell_Class, cell_classes))
    cellinfo <- cellinfo[idx,]
    expr <- expr[,idx]

    ## Restrict to cell classes with >1 cell types (containing cell types with >1 nuclei):
    
    temp <- cellinfo %>% 
      dplyr::group_by(Cell_Type) %>% 
      dplyr::filter(n()>1) %>%
      dplyr::group_by(MO_Cell_Class) %>% 
      dplyr::summarise(
        No.CTs=n_distinct(Cell_Type)
      ) %>%
      dplyr::filter(No.CTs>1)
    
    classes <- unique(temp$MO_Cell_Class)
    
    if(length(classes)>0){
      
      ## Calculate cell type summary vectors:
      
      expr <- scale(expr)
      
      rand_class_list <- lapply(1:length(classes), FUN=function(j){
        
        print(classes[j])
        
        ## Restrict to cell class cell types with >1 nuclei:
        
        temp <- cellinfo %>%
          dplyr::mutate(idx=1:n()) %>%
          dplyr:::filter(
            is.element(cellinfo$MO_Cell_Class, classes[j])
          ) %>%
          dplyr::group_by(Cell_Type) %>%
          dplyr::filter(n()>1)
        
        idx <- temp$idx
        expr1 <- expr[,idx]
        cellinfo1 <- cellinfo[idx,]
        
        cts <- unique(cellinfo1$Cell_Type)
        
        rand_ct_list <- future_lapply(1:length(cts), FUN=function(k){
        
          idx <- which(is.element(cellinfo1$Cell_Type, cts[k]))

          rand_list <- lapply(1:n_resamples, function(n){
            
            ## Sample same number of random cells as there are cell type cells
            ## Enforce ~equal # cells per cell type:

            n_samples <- ceiling(length(idx)/n_distinct(cellinfo1$Cell_Type))
            
            temp <- cellinfo1 %>%
              dplyr::mutate(idx=1:n()) %>%
              dplyr::group_by(Cell_Type) %>%
              dplyr::slice_sample(n=n_samples, replace=T)
          
            if(summary_type=="mean"){
              
              summary_vec <- rowMeans(expr1[,temp$idx])
              
            } else {
              
              summary_vec <- fast.prcomp(expr1[,temp$idx], center=F)$x[,1]
              
            }
            
            return(summary_vec)
            
          })
          
          rand_vecs <- do.call(cbind, rand_list)
          
          colnames(rand_vecs) <- rep(cts[k], ncol(rand_vecs))
          
          return(data.frame(SYMBOL=genes, rand_vecs, check.names=F))
          
        }, future.seed=T) ## rand_list <- future_lapply(
        
        return(rand_ct_list)
        
      }) ## rand_list <- lapply(
      
      names(rand_class_list) <- classes
      
      saveRDS(rand_class_list, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_CT_random_", summary_type, "cells_", data_type, "_", expr_type, "_", n_resamples, "_resamples_", nrow(expr), "_genes.RDS"))
      
    } ## if(length(classes)>0){
    
  } ## for(i in 1:nrow(datinfo)){
  
}

calc_class_summary_vecs <- function(datinfo, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "CGE", "MGE"), neu_subtypes, scaled=F){
  
  neu_subtypes$Label <- make.names(neu_subtypes$Label)
  idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
  neu_subtypes$Class_Level1[idx] <- NA

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
      
    }) 
    
    df <- data.frame(genes, do.call(cbind.data.frame, summary_list))
    colnames(df) <- c("SYMBOL", classes)
    
    file_path <- paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_cell_class_", summary_type, "cells_", data_type, "_", expr_type, "_", nrow(expr), "_genes.csv")
    
    if(scaled==F){
      file_path <- gsub(paste0(expr_type, "_"), paste0("unscaled_", expr_type, "_"), file_path)
    }
    
    fwrite(df, file=file_path)

  } ## for(i in 1:nrow(datinfo)){
  
} 

calc_class_random_vecs <- function(datinfo, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes=c("ASC", "EXC", "INH", "MIC", "OG", "OPC", "NEU", "CGE", "MGE"), neu_subtypes, n_resamples=10){
  
  neu_subtypes$Label <- make.names(neu_subtypes$Label)
  idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
  neu_subtypes$Class_Level1[idx] <- NA
  
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

    ## Add subtype info:
    
    temp <- make.names(paste0(datinfo$Dataset[i], "_", gsub(" ", "_", cellinfo$Cell_Type)))
    neu_subtypes1 <- neu_subtypes[match(temp, make.names(neu_subtypes$Label)),]
    idx <- which(!is.na(neu_subtypes1$Class_Level1))
    cellinfo$MO_Cell_Class[idx] <- neu_subtypes1$Class_Level1[idx]
    
    ## Restrict to cell classes present in dataset (with >1 nuclei):
    
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
    
    ## Calculate cell type summary vectors:
    
    expr <- scale(expr)   
    
    rand_list <- future_lapply(1:length(classes), FUN=function(j){
      
      if(is.element("NEU", classes[j])){
        idx <- which(is.element(cellinfo$MO_Cell_Class, c("NEU", "EXC", "INH", "CGE", "MGE")))
      } else if(is.element("INH", classes[j])) {
        idx <- which(is.element(cellinfo$MO_Cell_Class, c("INH", "CGE", "MGE")))
      } else {
        idx <- which(is.element(cellinfo$MO_Cell_Class, classes[j]))
      }
    
      random_class_list <- lapply(1:n_resamples, function(n){

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
      
      rand_vecs <- do.call(cbind.data.frame, random_class_list)
      colnames(rand_vecs) <- rep(classes[j], ncol(rand_vecs))

      return(data.frame(SYMBOL=genes, rand_vecs, check.names=F))
      
    }, future.seed=T) ## rand_list <- future_lapply(
    
    names(rand_list) <- classes
    
    saveRDS(rand_list, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_cell_class_random_", summary_type, "cells_", data_type, "_", expr_type, "_", n_resamples, "_resamples_", nrow(expr), "_genes.RDS"))
    
  } ## for(i in 1:nrow(datinfo)){
  
} 
