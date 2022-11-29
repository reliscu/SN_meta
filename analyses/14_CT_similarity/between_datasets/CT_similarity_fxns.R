library(dplyr)
library(data.table)
library(WGCNA)

CT_real_vs_rand_similarity <- function(datinfo, data_type, expr_type, sim_type=c("pearson"), cell_classes=c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "EXC", "INH"), summary_type, n_resamples){
  
  ## Make sure order of features match across datasets in real summary vectors:
  
  gene_list <- fread(datinfo$CT_Mean_Vecs[1], data.table=F)[,1]

  all_genes <- unlist(lapply(2:nrow(datinfo), function(i, gene_list){
    return(all.equal(gene_list, fread(datinfo$CT_Mean_Vecs[i], data.table=F)[,1]))
  }, gene_list))
  
  if(sum(!all_genes)>0){
    stop("Cell type summary vector features do not match")
  }
  
  ## Make sure order of features match across datasets in random summary vectors:

  gene_list <- readRDS(datinfo$CT_Random_Mean_Vecs[1])[[1]][[1]][,c(1)]
  
  all_genes <- unlist(lapply(1:nrow(datinfo), function(i, gene_list){
    temp <- Reduce(merge, readRDS(datinfo$CT_Random_Mean_Vecs[i]))
    mismatch <- unlist(lapply(grep("SYMBOL", colnames(temp)), function(n) !identical(temp[,n], temp[,1])))
    if(sum(mismatch)>0){
      stop("Order of genes in resamples do not match")
    }
    temp <- readRDS(datinfo$CT_Random_Mean_Vecs[i])[[1]][[1]]
    return(all.equal(gene_list, temp[,1]))
  }, gene_list))
  
  if(sum(!all_genes)>0){
    stop("Cell type random vector features do not match")
  }
  
  real_sim <- lapply(cell_classes, function(cell_class){
    
    cat("\n")
    print(cell_class)
    
    class_list <- lapply(1:nrow(datinfo), function(i){

      cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
      summary_vecs <- fread(datinfo$CT_Mean_Vecs[i], data.table=F)[,-c(1)]
      random_vecs <- readRDS(datinfo$CT_Random_Mean_Vecs[i])[[cell_class]]
      
      if(!is.null(random_vecs)){
        
        print(datinfo$Dataset[i])
        
        ## Restrict to cell types represented in both real and random:
        
        random_cts <- unlist(lapply(random_vecs, function(x) colnames(x)[2]))
        summary_vecs <- summary_vecs[,is.element(colnames(summary_vecs), random_cts)]
        
        ## Get cell types in current cell class:
        
        cellinfo <- cellinfo %>% 
          dplyr::filter(MO_Cell_Class==cell_class) %>%
          dplyr::filter(is.element(Cell_Type, colnames(summary_vecs))) ## Restrict to cell types represented by summary vectors
        
        idx <- which(
          is.element(colnames(summary_vecs), cellinfo$Cell_Type)
        )
        summary_vecs1 <- data.frame(summary_vecs[,idx])
        colnames(summary_vecs1) <- paste(datinfo$Dataset[i], colnames(summary_vecs)[idx])
        
        return(summary_vecs1)
        
      }
      
    }) ## class_list <- lapply(
    class_list <- class_list[lengths(class_list)>0]
    class_mat <- do.call(cbind, class_list)
    
    if(sim_type=="pearson"){
      
      sim_mat <- stats::cor(class_mat)
      
    } else {
      
      ## To do?
      
    }
    
  
    #sim_mat[lower.tri(sim_mat)] <- NA
    sim_df <- reshape2::melt(sim_mat, value.name="Similarity") %>% na.omit()
    
    dataset1 <- sapply(strsplit(as.character(sim_df$Var1), " "), "[", 1)
    dataset2 <- sapply(strsplit(as.character(sim_df$Var2), " "), "[", 1)
    study1 <- datinfo$Study[match(dataset1, datinfo$Dataset)]
    study2 <- datinfo$Study[match(dataset2, datinfo$Dataset)]
    
    sim_df <- sim_df %>% 
      dplyr::mutate(
        Dataset1=dataset1, Dataset2=dataset2,
        Study1=study1, Study2=study2
      ) %>%
      dplyr::filter(Study1!=Study2)

    return(sim_df)
    
  }) ## sim_list <- lapply(
  
  rand_sim <- lapply(cell_classes, function(cell_class){
    
    cat("\n")
    print(cell_class)
    
    rand_list <- lapply(1:nrow(datinfo), function(i){
      
      cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
      random_list <- readRDS(datinfo$CT_Random_Mean_Vecs[i])[[cell_class]]

      if(length(random_list)>0){
        
        print(datinfo$Dataset[i])
        
        random_list <- lapply(random_list, function(x) x[,-c(1)])
        random_vecs <- do.call(cbind, random_list)

        colnames(random_vecs) <- paste(datinfo$Dataset[i], colnames(random_vecs))
        return(random_vecs)
        
      } ## if(is.element(cell_class, vec_classes)){
      
    }) ## class_list <- lapply(
    rand_list <- rand_list[lengths(rand_list)>0]
    rand_mat <- do.call(cbind, rand_list)
    
    if(sim_type=="pearson"){
      
      sim_mat <- WGCNA::cor(rand_mat)
      
    } else {
      
      ## To do?
      
    }
    
    #sim_mat[lower.tri(sim_mat)] <- NA
    sim_df <- reshape2::melt(sim_mat, value.name="Similarity") %>% na.omit()
    
    dataset1 <- sapply(strsplit(as.character(sim_df$Var1), " "), "[", 1)
    dataset2 <- sapply(strsplit(as.character(sim_df$Var2), " "), "[", 1)
    study1 <- datinfo$Study[match(dataset1, datinfo$Dataset)]
    study2 <- datinfo$Study[match(dataset2, datinfo$Dataset)]
    
    sim_df <- sim_df %>% 
      dplyr::mutate(
        Dataset1=dataset1, Dataset2=dataset2,
        Study1=study1, Study2=study2
      ) %>%
      dplyr::filter(Study1!=Study2)
    
    return(sim_df)
    
  }) ## sim_list <- lapply(
  
  names(real_sim) <- names(rand_sim) <- cell_classes
  
  saveRDS(real_sim, file=paste0("data/", data_type, "/", expr_type, "/pairwise_CT_", summary_type, "cell_", toupper(sim_type), "_similarity_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies.RDS"))
  saveRDS(rand_sim, file=paste0("data/", data_type, "/", expr_type, "/pairwise_CT_random_", summary_type, "cell_", toupper(sim_type), "_similarity_", n_resamples, "_resamples_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies.RDS"))
  
}