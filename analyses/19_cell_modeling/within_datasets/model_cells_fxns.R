library(data.table)
library(dplyr)
library(gmodels) ## fast.prcomp()
#library(caret) ## createFolds()
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

model_class_real_vs_random <- function(datinfo, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes, neu_subtypes, n_resamples=10){

  neu_subtypes$Label <- make.names(neu_subtypes$Label)
  idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
  neu_subtypes$Class_Level1[idx] <- NA
  
  for(i in 1:nrow(datinfo)){ 
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr <- fread(datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
    summary_vecs <- fread(datinfo$Cell_Class_Mean_Vecs[i], data.table=F) 
    random_vecs <- readRDS(datinfo$Cell_Class_Random_Mean_Vecs[i]) 
    
    genes <- expr[,c(1)]
    expr <- expr[,-c(1)]
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
    classes <- cellinfo$MO_Cell_Class
    
    if(!is.element("EXC", classes)&is.element("EXC", cell_classes)){
      next
    }
    if(!is.element("INH", classes)&is.element("INH", cell_classes)){
      next
    }
    
    cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
    
    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("Order of cells in expression matrix does not match metadata!")
    }
    if(!identical(genes, summary_vecs$SYMBOL)){
      stop("Order of genes in expression matrix does not match summary cells!")
    }
    if(!identical(random_vecs[[1]]$SYMBOL, summary_vecs$SYMBOL)){
      stop("Order of genes in summary cells does not match random cells!")
    }
    
    ## Don't model NEU cells (cells not classified as EXC/INH) since they won't be modeled when we add EXC/INH predictors:
    
    idx <- which(!is.element(cellinfo$MO_Cell_Class, "NEU"))
    expr <- expr[,idx]
    cellinfo <- cellinfo[idx,]
    
    if(is.element("NEU", cell_classes)){
      
      idx <- which(is.element(cellinfo$MO_Cell_Class, c("INH", "EXC")))
      cellinfo$MO_Cell_Class[idx] <- "NEU"
      
    } else {
      
      ## Add subtype info:
      
      temp <- make.names(paste0(datinfo$Dataset[i], "_", gsub(" ", "_", cellinfo$Cell_Type)))
      neu_subtypes1 <- neu_subtypes[match(temp, make.names(neu_subtypes$Label)),]
      idx <- which(!is.na(neu_subtypes1$Class_Level1))
      cellinfo$MO_Cell_Class[idx] <- neu_subtypes1$Class_Level1[idx]
      
    }
    
    expr <- scale(expr)
  
    summary_vecs <- summary_vecs[,is.element(colnames(summary_vecs), cell_classes)]
  
    classes <- colnames(summary_vecs)
    
    if(!is.element("CGE", classes)&is.element("CGE", cell_classes)){
      next
    }
    
    ## Only model cells in classes represented by predictors:
    
    if(is.element("INH", classes)){
      idx <- which(
        is.element(cellinfo$MO_Cell_Class, unique(c(classes, "CGE", "MGE")))
      )
    } else if(is.element("CGE", classes)){
      idx <- which(
        is.element(cellinfo$MO_Cell_Class, unique(c(classes, "INH")))
      )
    } else {
      idx <- which(is.element(cellinfo$MO_Cell_Class, classes))
    }
    
    expr <- expr[,idx]
    cellinfo <- cellinfo[idx,]
    
    real_pred <- t(future_apply(expr, 2, FUN=function(cell_expr){
      mdl <- summary(lm(cbind(cell_expr, summary_vecs)))
      return(c(Real_R2=mdl$adj.r.squared, mdl$coefficients[-c(1),1]))
    }))
    
    real_coef <- real_pred[,-c(1)]
    colnames(real_coef) <- paste0("Real_", colnames(real_coef))
    
    random_pred <- lapply(1:n_resamples, function(n){
      rand_vecs <- do.call(cbind, lapply(random_vecs, function(x) x[,n+1]))
      rand_vecs <- rand_vecs[,is.element(colnames(rand_vecs), cell_classes)]
      rand_pred <- t(future_apply(expr, 2, FUN=function(cell_expr){
        mdl <- summary(lm(cbind.data.frame(cell_expr, rand_vecs)))
        return(c(Random_R2=mdl$adj.r.squared, mdl$coefficients[-c(1),1]))
      }))
      return(rand_pred)
    }) ## random_pred <- lapply(
    
    random_df <- do.call(cbind, random_pred)
    random_r2 <- random_df[,grep("R2", colnames(random_df))]
    random_coef <- random_df[,!grepl("R2", colnames(random_df))]
    random_coef <- random_coef[,order(colnames(random_coef))]
    colnames(random_coef) <- paste0("Random_", colnames(random_coef))
    
    r2_df <- data.frame(cellinfo, Real_R2=real_pred[,1], random_r2)
    coeff_df <- data.frame(cellinfo, real_coef, random_coef)

    fwrite(r2_df, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_cell_class_modeling_variance_explained_", paste(sort(classes), collapse="_"), "_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv"))
    
    fwrite(coeff_df, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_cell_class_modeling_coefficients_", paste(sort(classes), collapse="_"), "_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv"))
    
  }
  
}

model_CT_real_vs_random <- function(datinfo, data_type, expr_type, summary_type=c("mean", "eigen"), n_resamples=10){

  for(i in 1:nrow(datinfo)){ 
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr <- fread(datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
    summary_vecs <- fread(datinfo$CT_Mean_Vecs[i], data.table=F) 
    random_vecs <- readRDS(datinfo$CT_Random_Mean_Vecs[i]) 

    genes <- expr[,c(1)]
    expr <- expr[,-c(1)]
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
    cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
    
    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("Order of cells in expression matrix does not match metadata!")
    }
    if(!identical(genes, summary_vecs$SYMBOL)){
      stop("Order of genes in expression matrix does not match summary cells!")
    }
    if(!identical(random_vecs[[1]][[1]]$SYMBOL, summary_vecs$SYMBOL)){
      stop("Order of genes in summary cells does not match random cells!")
    }
    
    ## Restrict to real cell type predictors to those represented in random predictors:
    
    random_cts <- unlist(lapply(random_vecs, function(x) lapply(x, function(y) colnames(y)[2])))
    summary_vecs <- summary_vecs[,is.element(colnames(summary_vecs), random_cts)]
    
    ## Restrict to cell types represented by predictors:
    
    idx <- which(is.element(cellinfo$Cell_Type, colnames(summary_vecs)))
    expr <- expr[,idx]
    cellinfo <- cellinfo[idx,]
    
    expr <- scale(expr)
    
    real_pred <- t(future_apply(expr, 2, FUN=function(cell_expr){
      mdl <- summary(lm(cbind(cell_expr, summary_vecs)))
      return(c(Real_R2=mdl$adj.r.squared, mdl$coefficients[-c(1),1]))
    }))
    
    real_coef <- real_pred[,-c(1)]
    colnames(real_coef) <- paste0("Real_", colnames(real_coef))
    
    random_vecs <- unlist(random_vecs, recursive=F, use.names=F)

    random_pred <- lapply(1:n_resamples, function(n){
      rand_vecs <- do.call(cbind, lapply(random_vecs, function(x) x[,n+1]))
      colnames(rand_vecs) <- unlist(lapply(random_vecs, function(x) colnames(x)[n+1]))
      rand_pred <- t(future_apply(expr, 2, FUN=function(cell_expr){
        mdl <- summary(lm(cbind.data.frame(cell_expr, rand_vecs)))
        return(c(Random_R2=mdl$adj.r.squared, mdl$coefficients[-c(1),1]))
      }))
      return(rand_pred)
    }) ## random_r2 <- lapply(
    
    random_df <- do.call(cbind, random_pred)
    random_r2 <- random_df[,grep("R2", colnames(random_df))]
    random_coef <- random_df[,!grepl("R2", colnames(random_df))]
    random_coef <- random_coef[,order(colnames(random_coef))]
    colnames(random_coef) <- paste0("Random_", colnames(random_coef))
    
    r2_df <- data.frame(cellinfo, Real_R2=real_pred[,1], random_r2)
    coeff_df <- data.frame(cellinfo, real_coef, random_coef)
    
    fwrite(r2_df, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_CT_modeling_variance_explained_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv"))
    
    fwrite(coeff_df, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_CT_modeling_coefficients_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv"))
    
  }
  
}

model_class_vs_CT <- function(datinfo, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes, neu_subtypes){
  
  fwrite(r2_df, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_cell_class_vs_CT_modeling_variance_explained_", paste(sort(classes), collapse="_"), "_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv"))
  
  ## Make sure cell types are from the same cell class as cell classes being modeled
  
}