library(dplyr)
library(data.table)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

model_class_real_vs_random <- function(expr, cellinfo, neu_subtypes, summary_vecs, random_vecs, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes, n_resamples=10){
  
  if(!identical(colnames(expr), cellinfo$Cell_ID)){
    stop("Order of cells in expression matrix does not match metadata!")
  }
  if(!identical(rownames(expr), summary_vecs$SYMBOL)){
    stop("Order of genes in expression matrix does not match summary cells!")
  }
  if(!identical(random_vecs[[1]]$SYMBOL, summary_vecs$SYMBOL)){
    stop("Order of genes in summary cells does not match random cells!")
  }
  
  expr <- scale(expr)
  
  ## Restrict predictors:
  
  summary_vecs <- summary_vecs[,is.element(colnames(summary_vecs), cell_classes)]
  
  real_pred <- t(future_apply(expr, 2, FUN=function(cell_expr){
    mdl <- summary(lm(cbind(cell_expr, summary_vecs)))
    return(c(Real_R2=mdl$adj.r.squared, mdl$coefficients[-c(1),1]))
  }))
  
  real_r2 <- real_pred[,1]
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
  }) ## random_r2 <- lapply(
  
  random_df <- do.call(cbind, random_pred)
  random_r2 <- random_df[,grep("R2", colnames(random_df))]
  random_coef <- random_df[,!grepl("R2", colnames(random_df))]
  random_coef <- random_coef[,order(colnames(random_coef))]
  colnames(random_coef) <- paste0("Random_", colnames(random_coef))
  
  ## Maken sure cell class annotations match current predictors:
  
  if(sum(is.element(c("INH", "EXC"), cell_classes))==0){
    
    idx <- which(is.element(cellinfo$MO_Cell_Class, c("INH", "EXC")))
    cellinfo$MO_Cell_Class[idx] <- "NEU"
    
  } else {
    
    idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
    neu_subtypes$Class_Level1[idx] <- NA
    temp <- make.names(paste0(cellinfo$Dataset, "_", gsub(" ", "_", cellinfo$Cell_Type)))
    neu_subtypes <- neu_subtypes[match(temp, make.names(neu_subtypes$Label)),]
    cellinfo$MO_Cell_Class[!is.na(neu_subtypes$Class_Level1)] <- neu_subtypes$Class_Level1[!is.na(neu_subtypes$Class_Level1)]
    
  }
  
  r2_df <- data.frame(cellinfo, Real_R2=real_r2, random_r2)
  coeff_df <- data.frame(cellinfo, real_coef, random_coef)
  
  fwrite(r2_df, file=paste0("data/", data_type, "/", expr_type, "/modeling_variance_explained_", paste(cell_classes, collapse="_"), "_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv"))
  
  fwrite(coeff_df, file=paste0("data/", data_type, "/", expr_type, "/modeling_coefficients_", paste(cell_classes, collapse="_"), "_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv"))
  
}

model_class_real_vs_random_sequentially <- function(expr, cellinfo, neu_subtypes, summary_vecs, random_vecs, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes, n_resamples=10){
  
  idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
  neu_subtypes$Class_Level1[idx] <- NA
  temp <- make.names(paste0(cellinfo$Dataset, "_", gsub(" ", "_", cellinfo$Cell_Type)))
  neu_subtypes <- neu_subtypes[match(temp, make.names(neu_subtypes$Label)),]
  
  if(!identical(colnames(expr), cellinfo$Cell_ID)){
    stop("Order of cells in expression matrix does not match metadata!")
  }
  if(!identical(rownames(expr), summary_vecs$SYMBOL)){
    stop("Order of genes in expression matrix does not match summary cells!")
  }
  if(!identical(random_vecs[[1]]$SYMBOL, summary_vecs$SYMBOL)){
    stop("Order of genes in summary cells does not match random cells!")
  }
  
  ## Restrict predictors:
  
  summary_vecs <- summary_vecs[,is.element(colnames(summary_vecs), cell_classes)]
  
  ## Order cell classes by # nuclei:
  
  exc_inh <- table(cellinfo$MO_Cell_Class[is.element(cellinfo$MO_Cell_Class, c("EXC", "INH"))])
  cge_mge <- table(neu_subtypes$Class_Level1[!is.na(neu_subtypes$Class_Level1)])
  cellinfo$MO_Cell_Class[is.element(cellinfo$MO_Cell_Class, c("EXC", "INH"))] <- "NEU"
  
  predictors <- c(table(cellinfo$MO_Cell_Class), exc_inh, cge_mge)
  predictors <- names(rev(sort(predictors[predictors>0])))
  predictors <- predictors[is.element(predictors, cell_classes)]
  
  summary_vecs <- summary_vecs[,match(predictors, colnames(summary_vecs))]
  random_vecs <- random_vecs[match(predictors, names(random_vecs))]
  
  expr <- scale(expr)
  
  pred_list <- lapply(1:length(predictors), function(j){
    
    cat("\n")
    print(predictors[1:j])
    
    real_pred <- data.frame(future_apply(expr, 2, FUN=function(cell_expr){
      mdl <- summary(lm(cbind.data.frame(cell_expr, summary_vecs[,1:j])))
      return(mdl$adj.r.squared)
    }))
    
    colnames(real_pred) <- paste(colnames(summary_vecs)[1:j], collapse="_")

    random_pred <- lapply(1:n_resamples, function(n){
      rand_vecs <- do.call(cbind, lapply(random_vecs[1:j], function(x) x[,n+1]))
      rand_pred <- data.frame(future_apply(expr, 2, FUN=function(cell_expr){
        mdl <- summary(lm(cbind.data.frame(cell_expr, rand_vecs)))
        return(mdl$adj.r.squared)
      }))
      return(rand_pred)
    }) ## random_r2 <- lapply(
    
    random_pred <- do.call(cbind, random_pred)
    colnames(random_pred) <- rep(paste0("Random_", paste(colnames(summary_vecs)[1:j], collapse="_")), n_resamples)
    
    r2_df <- data.frame(cellinfo, real_pred, random_pred, check.names=F)
    
  })

  saveRDS(pred_list, file=paste0("data/", data_type, "/", expr_type, "/modeling_variance_explained_sequential_predictors_", paste(predictors, collapse="_"), "_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.RDS"))
  
}

model_class_left_out_dataset <- function(datinfo, cellinfo, summary_list, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes){
  
  datinfo <- datinfo[match(names(summary_list), datinfo$Dataset),]
  
  mdl_list <- lapply(1:nrow(datinfo), function(i){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    summary_vecs <- summary_list[[i]]
    expr <- fread(datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
    cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
    cellinfo <- cellinfo[match(colnames(expr)[-c(1)], cellinfo$Cell_ID),]
    
    if(!identical(colnames(expr)[-c(1)], cellinfo$Cell_ID)){
      stop("Order of cells in expression matrix does not match metadata!")
    }
    if(!identical(expr[,c(1)], summary_vecs$SYMBOL)){
      stop("Order of genes in expression matrix does not match summary cells!")
    }
    
    cellinfo <- cellinfo %>% dplyr::select(Cell_ID, Cell_Type, MO_Cell_Class)
    summary_vecs <- summary_vecs[,-c(1)]
    expr <- scale(expr[,-c(1)])
    
    mdl_df <- data.frame(t(future_apply(expr, 2, FUN=function(cell_expr){
      mdl <- summary(lm(cbind(cell_expr, summary_vecs)))
      return(c(R2=mdl$adj.r.squared, mdl$coefficients[-c(1),1]))
    })))
    
    ## Maken sure cell class annotations match current predictors:
    
    if(sum(is.element(c("INH", "EXC"), cell_classes))==0){
      
      idx <- which(is.element(cellinfo$MO_Cell_Class, c("INH", "EXC")))
      cellinfo$MO_Cell_Class[idx] <- "NEU"
      
    } else {
      
      idx <- which(!is.element(neu_subtypes$Class_Level1, cell_classes))
      neu_subtypes$Class_Level1[idx] <- NA
      temp <- make.names(paste0(cellinfo$Dataset, "_", gsub(" ", "_", cellinfo$Cell_Type)))
      neu_subtypes <- neu_subtypes[match(temp, make.names(neu_subtypes$Label)),]
      cellinfo$MO_Cell_Class[!is.na(neu_subtypes$Class_Level1)] <- neu_subtypes$Class_Level1[!is.na(neu_subtypes$Class_Level1)]
      
    }
  
    return(data.frame(Dataset=datinfo$Dataset[i], cellinfo, mdl_df))
    
  })
  
  names(mdl_list) <- datinfo$Dataset
  
  saveRDS(mdl_list, file=paste0("data/", data_type, "/", expr_type, "/modeling_", paste(cell_classes, collapse="_"), "_", summary_type, "cell_predictors_left_out_dataset_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.RDS"))
  
  
}

# model_CT_real_vs_random <- function(expr, cellinfo, summary_vecs, random_vecs, data_type, expr_type, summary_type=c("mean", "eigen"), cell_classes, n_resamples=10){
#   
#   if(!identical(colnames(expr), cellinfo$Cell_ID)){
#     stop("Order of cells in expression matrix does not match metadata!")
#   }
#   if(!identical(rownames(expr), summary_vecs$SYMBOL)){
#     stop("Order of genes in expression matrix does not match summary cells!")
#   }
#   if(!identical(random_vecs[[1]]$SYMBOL, summary_vecs$SYMBOL)){
#     stop("Order of genes in summary cells does not match random cells!")
#   }
#   
#   expr <- scale(expr)
#   
#   ## Confirm all cell types with predictors are represented in data:
#   
#   summary_vecs <- summary_vecs[,is.element(colnames(summary_vecs), cell_classes)]
#   
# }
