library(dplyr)
library(data.table)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

CT_vs_class_similarity <- function(datinfo, data_type, expr_type, sim_type=c("pearson"), summary_type){
  
  for(i in 1:nrow(datinfo)){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr <- fread(datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
    ct_vecs <- fread(datinfo$CT_Mean_Vecs_Unscaled[i], data.table=F)
    class_vecs <- fread(datinfo$Cell_Class_Mean_Vecs_Unscaled[i], data.table=F)
    
    genes <- expr[,c(1)]
    expr <- expr[,-c(1)]
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
    cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
    
    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("Order of cells in expression matrix does not match metadata!")
    }
    if(!identical(genes, ct_vecs$SYMBOL)){
      stop("Order of genes in expression matrix does not match cell type summary cells!")
    }
    if(!identical(genes, class_vecs$SYMBOL)){
      stop("Order of genes in expression matrix does not match cell class summary cells!")
    }
    
    idx <- is.element(cellinfo$Cell_Type, colnames(ct_vecs))
    cell_classes <- unique(cellinfo$MO_Cell_Class[idx])

    sim_list <- future_lapply(cell_classes, FUN=function(cell_class){

      ## Require at least 2 cell types within a given cell class:
      
      cts <- unique(cellinfo$Cell_Type[is.element(cellinfo$MO_Cell_Class, cell_class)])
      idx <- which(is.element(colnames(ct_vecs), cts))
      
      if(length(idx)>1){
        
        print(cell_class)
        
        ct_vecs1 <- data.frame(ct_vecs[,idx], check.names=F)
        class_vecs1 <- class_vecs[is.element(colnames(class_vecs), cell_class)]
        
        ## Calculate similarity of cell class nuclei to CELL TYPE vs. CELL CLASS summary vectors:
        
        cts <- colnames(ct_vecs1)
        idx <- which(is.element(cellinfo$Cell_Type, cts))
        expr1 <- expr[,idx]
        cellinfo1 <- cellinfo[idx,] %>% dplyr::select(Cell_ID, Cell_Type, MO_Cell_Class)
        
        ct_sim <- apply(ct_vecs1, 2, function(ct_vec){
          
          if(sim_type=="pearson"){
            
            ct_sim <- apply(expr1, 2, function(cell_expr) cor(ct_vec, cell_expr))
            
          } else {
            ## To do?
          }
          
          return(ct_sim)
          
        }) ## ct_sim <- apply(
        
        if(sim_type=="pearson"){
          
          class_sim <- apply(expr1, 2, function(cell_expr) cor(class_vecs1[,1], cell_expr))
          
        } else {
          ## To do?
        }
        
        temp <- data.frame(cellinfo1, ct_sim, class_sim, check.names=F)
        colnames(temp)[ncol(temp)] <- cell_class
        df <- reshape2::melt(temp, value.name="Similarity", variable.name="Distro")
        
        return(data.frame(Dataset=datinfo$Dataset[i], df, check.names=F))
        
      } ## if(length(idx)>1){
      
    }, future.seed=T) ## class_list <- lapply(
    
    names(sim_list) <- cell_classes
    sim_list <- sim_list[lengths(sim_list)>0]
    
    saveRDS(sim_list, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_CT_vs_cell_class_similarity_", data_type, "_", expr_type, "_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies_", n_resamples, "_resamples.RDS"))
    
  }
  
}
