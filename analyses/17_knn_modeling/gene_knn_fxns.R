library(dplyr)
library(data.table)
library(Matrix)
library(gtools) ## permute()
library(WGCNA)
library(future.apply)
library(openblasctl)
library(tictoc)

openblas_set_num_threads(10)
options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

model_knn_vs_rand <- function(
  datinfo, 
  data_type, 
  expr_type, 
  n_neighbors=c(1, 2, 5, 10, 20, 50), 
  n_resamples, 
  cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"), 
  pc_genes=NULL
){
  
  for(i in 1:nrow(datinfo)){
    
    tic(datinfo$Dataset[i])
    
    expr_path <- datinfo$Author_Normalized_Expression[i]
    
    if(grepl("mtx", expr_path)){
      
      expr <- readMM(expr_path)
      
    } else {
      
      expr <- fread(expr_path, data.table=F)
      
    }
    
    expr <- as(expr, "matrix")
    
    genes <- fread(datinfo$Author_Genes_Mapped[i], data.table=F) %>% na_if("")
    barcodes <- fread(datinfo$Author_Barcodes_QC[i], header=F, data.table=F) %>% na_if("")
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F) %>% na_if("")
    
    ## Subset to cell class cells:
    
    if(nrow(barcodes)<ncol(expr)){
      barcodes <- fread(datinfo$Author_Barcodes[i], header=F, data.table=F)
    }
    
    cellinfo <- cellinfo[match(barcodes[,1], cellinfo$Cell_ID),]
    
    if(!identical(barcodes[,1], cellinfo$Cell_ID)){
      stop("!identical(barcodes[,1], cellinfo$Cell_ID)")
    }
    
    idx <- which(
      is.element(cellinfo$MO_Cell_Class, cell_classes)
    )
    cellinfo <- cellinfo[idx,]
    expr <- expr[,idx]
    
    if(!is.null(pc_genes)){
      
      ## Subset to protein coding genes:
      
      idx <- which(
        is.element(genes$SYMBOL, pc_genes$SYMBOL)
      )
      expr <- expr[idx,]
      genes <- genes[idx,]
      
    }
    
    ## Remove zero variance genes:
    
    idx <- rowSums(expr)>0
    genes <- genes[idx,]
    expr <- expr[idx,]
    
    if(sum(duplicated(genes$SYMBOL))>0){
      
      ## Use gene with highest mean expression:
      
      mean_expr <- genes %>%
        dplyr::mutate(
          Mean_Expr=rowMeans(expr), idx=1:n()
        ) %>% 
        dplyr::group_by(SYMBOL) %>%
        dplyr::slice_max(
          Mean_Expr
        )
      
      idx <- mean_expr$idx
      expr <- expr[idx,]
      genes <- genes[idx,]
      
    }
    
    sim_mat <- WGCNA::cor(t(expr))
    diag(sim_mat) <- 0
    
    ## Find max(k) nearest neighbors per gene:
    
    neighbor_idx <- apply(sim_mat, 2, function(x){
      order(-x)[1:max(n_neighbors)]
    })
    
    rm(sim_mat)
    
    ## Model each gene with its top k neighbors
    
    expr <- t(scale(t(expr)))

    nn_r2_per_k <- future_lapply(n_neighbors, FUN=function(k){
      unlist(lapply(1:nrow(expr), function(i){
        idx <- neighbor_idx[1:k,i]
        if(k==1){
          summary(lm(expr[i,]~expr[idx,]))$adj.r.squared
        } else {
          summary(lm(expr[i,]~colMeans(expr[idx,])))$adj.r.squared
        }
      }))
    }) ## for(k in n_neighbors){
    
    names(nn_r2_per_k) <- as.character(n_neighbors)
    
    ## Sample max(k) random "neighbors" per gene n times:
    
    rand_idx_list <- vector(mode="list", length=n_resamples)

    for(n in 1:n_resamples){
      rand_idx_sample <- lapply(1:nrow(expr), FUN=function(i){
        idx <- 1:nrow(expr) 
        sample(idx[idx!=i], size=max(k))
      })
      rand_idx_list[[n]] <- do.call(cbind, rand_idx_sample)
    } 
    
    ## Model each gene with k random neighbors n times:

    rand_r2_per_k <- future_lapply(n_neighbors, FUN=function(k){
      
      rand_r2_resample <- lapply(1:n_resamples, function(i){
        
        rand_idx <- rand_idx_list[[i]]
        
        r2_per_gene <- unlist(lapply(1:nrow(expr), function(j){
          rand_idx_gene <- rand_idx[1:k, j]
          if(k==1){
            summary(lm(expr[j,]~expr[rand_idx_gene,]))$adj.r.squared
          } else {
            summary(lm(expr[j,]~colMeans(expr[rand_idx_gene,])))$adj.r.squared
          }
        }))
        
        return(r2_per_gene)
        
      }) ## rand_r2_resample <- lapply(
      
      return(do.call(cbind, rand_r2_resample))
      
    }) ## rand_r2_per_k <- future_lapply(
    
    names(rand_r2_per_k) <- as.character(n_neighbors)
   
    ## Calc fraction of random r2 >= real r2:
    
    emp_pval_per_k <- vector(mode="list", length=length(n_neighbors))
    names(emp_pval_per_k) <- as.character(n_neighbors)
    
    for(k in n_neighbors){
      nn_r2 <- nn_r2_per_k[[as.character(k)]]
      rand_r2 <- rand_r2_per_k[[as.character(k)]]
      emp_pval_per_k[[as.character(k)]] <- unlist(
        lapply(1:length(nn_r2), function(i) sum(rand_r2[1,]>=nn_r2[i])/n_resamples
        ))
    }
    
    emp_pval_df <- do.call(cbind, emp_pval_per_k)
    colnames(emp_pval_df) <- paste0("k_", n_neighbors)
    
    emp_pval_df <- data.frame(
      SYMBOL=genes$SYMBOL, 
      emp_pval_df
    )
    
    fwrite(emp_pval_df, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_gene_knn_vs_random_modeling_empirical_pvals_", n_resamples, "_resamples_", data_type, "_", expr_type, ".csv"))
    
    toc()
    
  } ## for(i in 1:nrow(datinfo)){
  
} ## knn_model_vs_rand <- function(

nn_vs_chance_sim <- function(
  datinfo, 
  data_type, 
  expr_type, 
  n_resamples, 
  cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"),  
  pc_genes=NULL
){

  for(i in 1:nrow(datinfo)){

    print(datinfo$Dataset[i])

    expr_path <- datinfo$Author_Normalized_Expression[i]

    if(grepl("mtx", expr_path)){

      expr <- readMM(expr_path)

    } else {

      expr <- fread(expr_path, data.table=F)

    }

    expr <- as(expr, "matrix")

    genes <- fread(datinfo$Author_Genes_Mapped[i], data.table=F) %>% na_if("")
    barcodes <- fread(datinfo$Author_Barcodes_QC[i], header=F, data.table=F) %>% na_if("")
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F) %>% na_if("")
    
    ## Subset to cell class cells:
    
    if(nrow(barcodes)<ncol(expr)){
      barcodes <- fread(datinfo$Author_Barcodes[i], header=F, data.table=F)
    }
    
    cellinfo <- cellinfo[match(barcodes[,1], cellinfo$Cell_ID),]
    
    if(!identical(barcodes[,1], cellinfo$Cell_ID)){
      stop("!identical(barcodes[,1], cellinfo$Cell_ID)")
    }
    
    idx <- which(
      is.element(cellinfo$MO_Cell_Class, cell_classes)
    )
    cellinfo <- cellinfo[idx,]
    expr <- expr[,idx]
    
    if(!is.null(pc_genes)){
      
      ## Subset to protein coding genes:
      
      idx <- which(
        is.element(genes$SYMBOL, pc_genes$SYMBOL)
      )
      expr <- expr[idx,]
      genes <- genes[idx,]
      
    }
    
    ## Remove zero variance genes:
    
    idx <- rowSums(expr)>0
    genes <- genes[idx,]
    expr <- expr[idx,]
    
    if(sum(duplicated(genes$SYMBOL))>0){
      
      ## Use gene with highest mean expression:
      
      mean_expr <- genes %>%
        dplyr::mutate(
          Mean_Expr=rowMeans(expr), 
          idx=1:n()
        ) %>% 
        dplyr::group_by(SYMBOL) %>%
        dplyr::slice_max(
          Mean_Expr
        )
      
      idx <- mean_expr$idx
      expr <- expr[idx,]
      genes <- genes[idx,]
      
    }

    sim_mat <- WGCNA::cor(t(expr))
    diag(sim_mat) <- 0

    ## Find nearest neighbor per each gene

    neighbor_vec <- apply(sim_mat, 2, max)

    rm(sim_mat)

    ## Permute counts for each gene and find nearest neighbor:

    perm_neighbor_list <- lapply(1:n_resamples, function(i){
      sim_perm <- WGCNA::cor(apply(expr, 1, permute))
      return(apply(sim_perm, 2, max))
    })

    perm_neighbor_df <- do.call(cbind, perm_neighbor_list)

    ## Calc fraction of perm sim >= real sim:

    emp_pval <- unlist(lapply(1:length(neighbor_vec), function(x){
      sum(perm_neighbor_df[i,]>=neighbor_vec[i])/n_resamples
    }))

    emp_pval_df <- data.frame(
      SYMBOL=genes$SYMBOL,
      Empirical_Pval=emp_pval
    )

    fwrite(emp_pval_df, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_gene_nn_vs_chance_similarity_empirical_pvals_", n_resamples, "_resamples_", data_type, "_", expr_type, "_", nrow(datinfo), "_datasets.csv"))

  } ## for(i in 1:nrow(datinfo)){

} ## nn_vs_chance <- function(
