library(dplyr)
library(Matrix)
library(data.table)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=5)

source("/home/rebecca/code/misc/rank_percentile.R")

source("/home/rebecca/SCSN_meta_analysis/code/fidelity_fxn.R")

cell_class_fidelity <- function(expr, cellinfo, data_type, expr_type, cell_classes=c("ASC", "END", "NEU", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"), prop_scaled=F){
  
  if(!identical(colnames(expr), cellinfo$Cell_ID)){
    stop("Order of cells in expression matrix does not match metadata!")
  }

  fid_list <- future_lapply(1:length(cell_classes), FUN=function(j){
    
    idx <- which(is.element(cellinfo$MO_Cell_Class, cell_classes[j]))
    if(cell_classes[j]=="NEU"){
      idx <- which(is.element(cellinfo$MO_Cell_Class, c("NEU", "EXC", "INH")))
    }
    
    fid <- fidelity(expr, idx, prop_scaled)
    expr_percentile <- rank_percentile(rowMeans(expr[,idx]))
      
    return(data.frame(Fidelity=fid, Mean_Expr_Percentile=expr_percentile))
    
  }) 
  
  df <- do.call(cbind.data.frame, fid_list)
  colnames(df) <- paste(classes, colnames(df))
  df <- data.frame(SYMBOL=genes, df)
  
  file_path <- paste0("data/", data_type, "/", expr_type, "/cell_class_fidelity_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", nrow(expr), "_genes.csv")
  
  if(prop_scaled){
    file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
  } 
  
  fwrite(df, file=file_path)
  
}