library(dplyr)
library(Matrix)
library(data.table)
library(future.apply)
library(lsa) ## cosine()

options(future.globals.maxSize=Inf)
plan(multicore, workers=12)

CT_pairwise_similarity <- function(datinfo, summary_vec_paths, random_vec_paths, expr_type, sim_type, cell_classes, summary_type=c("mean", "eigen"), n_resamples=10){
  
  for(i in 1:length(summary_vec_paths)){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    summary_vecs <- fread(summary_vec_paths[i], data.table=F)
    summary_vecs <- summary_vecs[,c(1, grep("_PC1", colnames(summary_vecs)))]
    colnames(summary_vecs) <- gsub("_PC1", "", colnames(summary_vecs))
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F)
    
    cellinfo <- cellinfo %>% dplyr::arrange(Cell_Type) %>% dplyr::filter(!duplicated(Cell_Type))
    cellinfo$Cell_Type <- make.names(gsub(" ", "_", cellinfo$Cell_Type))
    idx <- match(
      colnames(summary_vecs)[-c(1)], cellinfo$Cell_Type
    )
    cellinfo <- cellinfo[idx,]
    
    if(!identical(colnames(summary_vecs)[-c(1)], cellinfo$Cell_Type)){
      stop("!identical(colnames(summary_vecs)[-c(1)], cellinfo$Cell_Type")
    }

    summary_vec_sim_list <- lapply(1:length(summary_vec_paths), function(j){
      
      if(i!=j){
        
        summary_vecs1 <- fread(summary_vec_paths[j], data.table=F)
        summary_vecs1 <- summary_vecs1[,c(1, grep("_PC1", colnames(summary_vecs1)))]
        colnames(summary_vecs1) <- gsub("_PC1", "", colnames(summary_vecs1))
        
        ## Make sure order of features match between cells and summary_vecs:
        
        idx <- match(summary_vecs[,c(1)], summary_vecs1[,c(1)])
        summary_vecs1 <- summary_vecs1[idx,]
        
        summary_vec_sim <- future_apply(summary_vecs1[,-c(1)], 2, FUN=function(summary_vec1){
          apply(summary_vecs[,-c(1)], 2, function(summary_vec){
            if(sim_type=="pearson"){
              stats::cor(summary_vec1, summary_vec)
            } else {
              as.numeric(cosine(summary_vec1, summary_vec))
            } 
          })
        }) ## summary_vec_sim <- future_apply(
        
        return(data.frame(
          Dataset1=datinfo$Dataset[i],
          Dataset2=datinfo$Dataset[j],
          Cell_Type=rownames(summary_vec_sim),
          Cell_Class=cellinfo$MO_Cell_Class,
          summary_vec_sim
        ))
        
      }
      
    }) ## for(j in 1:length(summary_vec_paths)){
    
    summary_vec_sim_list <- summary_vec_sim_list[!unlist(lapply(summary_vec_sim_list, is.null))]
    
    saveRDS(summary_vec_sim_list, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_", summary_type, "cell_vs_", summary_type, "cell_", toupper(sim_type), "_similarity_", data_type, "_", expr_type, "_", nrow(summary_vecs), "_genes.RDS"))
    
  } ## for(i in 1:length(summary_vec_paths)){
  
} ## summary_vec_pairwise_similarity <- function(

summary_vec_similarity_stats <- function(datinfo, sim_paths, data_type, expr_type){
  
  stat_list <- lapply(1:length(sim_paths), function(i){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    summary_vec_sim_list <- readRDS(sim_paths[i])
    
    df_list <- lapply(1:length(summary_vec_sim_list), function(j){
      
      summary_vec_sim <- summary_vec_sim_list[[j]]
      
      max_summary_vec <- lapply(1:nrow(summary_vec_sim), function(i){
        ## Sign of PC1 is arbitrary -- we assume all cells are positively correlated to other summary_vecs:
        cell_sim <- abs(summary_vec_sim[i,-c(1:4)]) 
        idx <- which.max(cell_sim)
        max_summary_vec <- names(cell_sim)[idx]
        max_sim <- as.numeric(cell_sim[idx])
        return(data.frame(
          Max_summary_vec=max_summary_vec, 
          Max_Sim=max_sim
        ))
      }) ## max_summary_vec <- lapply(
      
      max_df <- do.call(rbind, max_summary_vec)
      return(data.frame(
        summary_vec_sim[,c(1:4)], max_df
      ))
      
    })
    df <- do.call(rbind, df_list)
    
  }) ## future_lapply(1:length(sim_paths)
  
  names(stat_list) <- datinfo$Dataset
  
  file_path <- paste0("data/", data_type, "/", expr_type, "/summary_vec_vs_summary_vec_", toupper(sim_type), "_similarity_stats_", data_type, "_", expr_type, "_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies.RDS")
  
  saveRDS(stat_list, file=file_path)
  
}

summary_vec_MNN <- function(stats_list, sim_type, data_type, expr_type){
  
  ## For each dataset: Look at dataset 2:
  
  nn_list <- lapply(1:length(stats_list), function(i){
    
    sim <- stats_list[[i]]
    
    cat("\n")
    print(sim$Dataset1[i])
    
    datasets2 <- unique(sim$Dataset2)
    
    nn_list <- lapply(1:length(datasets2), function(j){
      
      ## Restrict to one comparison dataset:
      
      sim1 <- sim[is.element(sim$Dataset2, datasets2[j]),]
      
      idx <- which(
        is.element(names(stats_list), sim1$Dataset2[1])
      )
      
      sim2 <- stats_list[[idx]]
      
      idx <- which(
        is.element(sim2$Dataset2, sim1$Dataset1[1])
      )
      
      sim2 <- sim2[idx,]
      
      nn_dat_list <- lapply(1:nrow(sim1), function(i){
        
        ## What is the nearest neighbor for cell type from dataset1 in dataset2?
        summary_vec <- sim1[i,]
        nn1 <- summary_vec$Max_summary_vec
        ## What is the nearest neighbor cell type's nearest neighbor?
        idx <- which(
          is.element(sim2$Cell_Type, nn1)
        )
        nn2 <- sim2$Max_summary_vec[idx]
        
        return(data.frame(
          Dataset1_CT=summary_vec$Cell_Type,
          Dataset1_Class=summary_vec$Cell_Class,
          Dataset2_Class=sim2$Cell_Class[idx],
          Dataset1_NN=nn1,
          Dataset2_NN=nn2,
          Sim=summary_vec$Max_Sim
        ))
        
      }) ## nn_list <- lapply(1:nrow(sim1)
      
      nn_dat <- do.call(rbind, nn_dat_list)
      
      return(data.frame(
        Dataset1=sim1$Dataset1[1],
        Dataset2=sim2$Dataset1[1],
        Dataset1_No.CTs=nrow(sim1),
        Dataset2_No.CTs=nrow(sim2),
        nn_dat
      ))
      
    }) ## lapply(1:length(datasets), function(j){
    
    df <- do.call(rbind, nn_list)
    
    return(df)
  
  }) ## nn_list <- lapply(

  df <- do.call(rbind, nn_list)
  
  file_path <- paste0("data/", data_type, "/", expr_type, "/", summary_type, "cell_MNN_", toupper(sim_type), "_similarity_", data_type, "_", expr_type, "_", length(stats_list), "_datasets.csv")

  fwrite(df, file=file_path)
  
} ## summary_vec_MNN <- function(

