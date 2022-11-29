library(dplyr)
library(Matrix)
library(data.table)
library(future.apply)
library(lsa) ## cosine()

source("/home/rebecca/code/misc/normalize_functions.R")

options(future.globals.maxSize=Inf)
plan(multicore, workers=5)

eigencell_similarity <- function(datinfo, expr_list_paths, eigencell_list_paths, data_type, expr_type, sim_type, top_n){
  
  for(i in 1:length(eigencell_list_paths)){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr_list <- readRDS(expr_list_paths[i])
    
    eigencell_list <- readRDS(eigencell_list_paths[i])
    
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
    
    genes <- rownames(expr)
    
    rm(expr)
    
    ## Ensure cell order in expr_list matches cell metadata:
    
    match_cells <- unlist(lapply(expr_list, function(x) !identical(colnames(x), cellinfo$Cell_ID)))
    
    if(sum(match_cells)>0){
      stop("Cells in expr_list do not match cell metadata")
    }
    
    sim_list <- future_lapply(1:length(expr_list), function(j){
      
      expr <- expr_list[[j]]
      
      eigencells <- eigencell_list[[j]]
      eigencells <- eigencells[,c(1, grep("PC1", colnames(eigencells)))]
      colnames(eigencells) <- gsub("_PC1", "", colnames(eigencells))
      
      ## Restrict expression data to cell types represented by eigencells:
      
      idx <- which(
        is.element(make.names(gsub(" ", "_", cellinfo$Cell_Type)), colnames(eigencells))
      )
      
      expr <- expr[,idx]
      cellinfo1 <- cellinfo[idx,]

      if(!identical(colnames(expr), cellinfo1$Cell_ID)){
        stop("!identical(colnames(expr), cellinfo$Cell_ID)")
      }
      
      expr <- log2(normalize_fxn(as(expr, "matrix"), scale_factor=1e4)+1)
      
      ## Make sure order of features match between cells and eigencells:
      
      idx <- match(rownames(expr), eigencells[,c(1)])
      eigencells <- eigencells[idx,]
      
      eigencell_sim <- apply(eigencells[,-c(1)], 2, function(eigencell){
        apply(expr, 2, function(cell_expr){
          if(sim_type=="pearson"){
            stats::cor(eigencell, cell_expr)
          } else {
            as.numeric(cosine(eigencell, cell_expr))
          } 
        })
      })
      
      colnames(eigencell_sim) <- colnames(eigencells)[-c(1)]
      
      df <- data.frame(
        Dataset=datinfo$Dataset[i],
        Threshold=names(expr_list)[j],
        Cell_ID=cellinfo1$Cell_ID,
        Cell_Type=cellinfo1$Cell_Type,
        Cell_Class=cellinfo1$MO_Cell_Class,
        eigencell_sim
      )
      
    }) ## lapply(1:length(expr_list), function(j){
    
    sim_df <- do.call(rbind, sim_list)
    
    fwrite(sim_df, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_cell_vs_eigencell_", toupper(sim_type), "_similarity_UMIs_downsampled_", data_type, "_", expr_type, "_", length(genes), "_genes.csv"))
    
  } ## for(i in 1:length(eigencell_list_paths)){
  
} ## eigencell_similarity <- function(

eigencell_similarity_stats <- function(datinfo, eigencell_sim_list, data_type, expr_type){
  
  stat_list <- lapply(eigencell_sim_list, function(eigencell_sim){
    
    cat("\n")
    print(eigencell_sim$Dataset[1])
    
    eigencell_sim$Cell_Type <- make.names(
      gsub(" ", "_", eigencell_sim$Cell_Type)
    )
    
    max_eigencell <- lapply(1:nrow(eigencell_sim), function(i){
      ## Sign of PC1 is arbitrary -- we assume all cells are positively correlated to their eigencell:
      cell_sim <- abs( 
        eigencell_sim[i,-c(1:5)]
      ) 
      idx <- which.max(cell_sim)
      max_eigencell <- names(cell_sim)[idx]
      max_sim <- as.numeric(cell_sim[idx])
      idx <- which(
        is.element(colnames(cell_sim), eigencell_sim$Cell_Type[i])
      )
      true_sim <- as.numeric(cell_sim[idx])
      return(data.frame(
        Max_Eigencell=max_eigencell, 
        Max_Sim=max_sim, 
        True_Sim=true_sim
      ))
    }) ## max_eigencell <- lapply(
    
    df <- do.call(rbind, max_eigencell)
    
    return(
      data.frame(eigencell_sim[,c(1:5)], df)
    )
    
  }) ## future_lapply(1:length(eigen_sim_paths)
  
  saveRDS(stat_list, file=paste0("data/", data_type, "/", expr_type, "/", "cell_vs_eigencell_", toupper(sim_type), "_similarity_stats_UMIs_downsampled_", data_type, "_", expr_type, "_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies.RDS"))
  
}