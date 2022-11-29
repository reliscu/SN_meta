library(dplyr)
library(Matrix)
library(data.table)
library(future.apply)
library(openblasctl)

openblas_set_num_threads(10)
options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

source("/home/rebecca/code/misc/rank_percentile.R")
source("/home/rebecca/code/misc/normalize_functions.R")

cell_class_fidelity <- function(
  datinfo, 
  expr_list,
  data_type, 
  expr_type,
  cell_classes=c("ASC", "END", "NEU", "MIC", "OG", "OPC", "PER", "VSMC"),
  prop_scaled=F
){
  
  for(i in 1:length(expr_list_paths)){
    
    cat("\n")
    
    print(datinfo$Dataset[i])
    print(expr_list_paths[[i]])
    
    expr_list <- readRDS(expr_list_paths[[i]])
    
    ## Find intersection of nonzero genes:
    
    idx_list <- lapply(expr_list, function(expr) which(rowSums(expr)>0))
    idx <- Reduce(intersect, idx_list)
    expr_list <- lapply(expr_list, function(expr){
      expr[idx,]
    })
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F) %>% na_if("")
    
    ## Ensure cell order in expr_list matches cell metadata:
    
    expr <- expr_list[[1]]
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]
    
    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("!identical(colnames(expr), cellinfo$Cell_ID)")
    }
    
    genes <- rownames(expr)
    
    rm(expr)
    
    match_cells <- unlist(lapply(expr_list, function(x) !identical(colnames(x), cellinfo$Cell_ID)))
    
    if(sum(match_cells)>0){
      stop("Cells in expr_list do not match cell metadata")
    }
    
    expr_counts_list <- expr_list
    
    ## Normalize:
    
    expr_list <- future_lapply(expr_list, FUN=function(expr){
      log2(normalize_fxn(expr, scale_factor=1e4)+1)
    })
    
    ## Calculate fidelity for each cell class:
    
    fidelity_class_list <- vector(mode="list", length=length(cell_classes))
    names(fidelity_class_list) <- cell_classes
    
    for(j in 1:length(cell_classes)){
      
      cell_class <- cell_classes[j]
      
      cat("\n")
      print(cell_class)
      
      if(cell_class=="NEU"){
        cellinfo$MO_Cell_Class[is.element(cellinfo$MO_Cell_Class, c("INH", "EXC"))] <- "NEU"
        cellinfo$MO_Cell_Class <- toupper(cellinfo$MO_Cell_Class)
      }
      
      idx <- which(is.element(cellinfo$MO_Cell_Class, cell_class))
      
      ## Subset to cell class cells for each threshold:
      
      class_expr <- lapply(expr_list, function(expr){
        expr[,idx]
      })
      
      ## Calculate fidelity for each threshold:
      
      class_sensitivity <- future_lapply(class_expr, FUN=function(expr){
        apply(expr, 1, function(x) sum(x>0)/length(x))
      })
      class_total_expr <- lapply(class_expr, rowSums)
      total_expr <- lapply(expr_list, function(expr) rowSums(expr))
      class_specificity <- lapply(1:length(expr_list), function(i){
        class_total_expr[[i]]/total_expr[[i]]
      }) 
      
      class_fidelity <- lapply(1:length(expr_list), function(i) class_specificity[[i]]*class_sensitivity[[i]]) 
      
      if(prop_scaled){
        
        class_all_gene_expr <- lapply(class_total_expr, sum)
        
        ## Probability a transcript came from the gene of interest:
        
        class_prop <- lapply(1:length(expr_list), function(i){ 
          rowSums(class_expr[[i]])/class_all_gene_expr[[i]]
        })
        
        class_fidelity <- lapply(1:length(expr_list), function(i) class_fid[[i]]*class_prop[[i]]) 
        
      } ## if(prop_scaled){
      
      fidelity_df_list <- lapply(1:length(expr_list), function(i){
        
        class_mean_expr_perc <- rank_percentile(
          rowMeans(class_expr[[i]])
        )
        class_mean_umis <- rowMeans(expr_counts_list[[i]])
        
        fidelity_df <- data.frame(
          SYMBOL=genes, 
          Fidelity=class_fidelity[[i]], 
          Mean_Expr_Percentile=class_mean_expr_perc, 
          Mean_UMIs=class_mean_umis, 
          No.Nuclei=length(idx),
          Threshold=names(expr_list)[i],
          Dataset=datinfo$Dataset[i]
        )
        
        return(fidelity_df)
        
      }) ## fidelity_df_df_list <- lapply(
      
      fidelity_df <- do.call(rbind, fidelity_df_list)
      
      fidelity_class_list[[j]] <- fidelity_df
      
    } ## for(j in 1:length(cell_classes)){
    
    file_path <- paste0(
      "data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_cell_class_fidelity_", data_type, "_", expr_type, "_", length(genes), "_genes.RDS"
    )
    
    if(prop_scaled){
      file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
    }
    
    saveRDS(fidelity_class_list, file=file_path)
    
  }
  
} ## cell_class_fidelity


