library(dplyr)
library(data.table)
library(future.apply)
library(gmodels) # fast.prcomp()

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

eigencell_VE <- function(datinfo, data_type, expr_type, cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")){
  
  df_list <- lapply(1:nrow(datinfo), function(i){
    
    print(datinfo$Dataset[i])
    
    expr <- fread(datinfo$Author_Mean_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
    expr <- expr[,-c(1)]
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F)
    cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
    
    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("Order of cells in expression matrix does not match metadata!")
    }
    
    ## Restrict to cell classes of interest and cell types with at least 2 nuclei:
    
    temp <- cellinfo %>%
      dplyr::mutate(idx=1:n()) %>%
      dplyr:::filter(
        is.element(cellinfo$MO_Cell_Class, cell_classes)
      ) %>%
      dplyr::group_by(Cell_Type) %>%
      dplyr::filter(n()>1)
    
    cts <- unique(temp$Cell_Type)
    
    ve_list <- future_lapply(1:length(cts), FUN=function(j){
      
      idx <- which(is.element(cellinfo$Cell_Type, cts[j]))
      ct_pcs <- fast.prcomp(expr[,idx], center=T, scale=T, retx=F)
      return(data.frame(Cell_Type=cts[j], PC1_VE=ct_pcs$sdev^2/sum(ct_pcs$sdev^2)[1]))
    
    }) 
    
    names(ve_list) <- cts
    
    return(data.frame(Dataset=datinfo$Dataset[i], do.call(rbind, ve_list), check.names=F))
    
  }) 
  
  df <- do.call(rbind, df_list)
  
  fwrite(df, file=paste0("data/", data_type, "/", expr_type, "/eigencell_variance_explained_", data_type, "_", expr_type, "_", length(df_list), "_datasets_", n_distinct(datinfo$Study), "_studies.csv"))
  
} ## eigencell_VE <- function(
