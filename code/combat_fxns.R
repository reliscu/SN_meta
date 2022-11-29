library(sva)
library(data.table)

combat_fxn <- function(datinfo, expr, data_type, expr_type){
  
  datinfo <- datinfo %>% na_if("") %>% dplyr::filter(is.element(Dataset, colnames(expr)))

  datinfo <- datinfo[match(colnames(expr)[-c(1)], datinfo$Dataset),]
  
  ## Remove genes that are NA in ANY dataset:
  
  expr <- expr[apply(expr[,-c(1)], 1, function(x) sum(is.na(x)))==0,]
  
  ## Remove zero variance genes:
  
  expr <- expr[apply(expr[,-c(1)], 1, function(x) sum(x))>0,]
  
  batch_vec <- as.factor(datinfo$Platform)
  
  expr_combat <- sva::ComBat(expr[,-c(1)], batch=batch_vec)
  
  expr <- data.frame(SYMBOL=expr[,c(1)], expr_combat)
  
  fwrite(expr, file=paste0("data/", data_type, "/", expr_type, "/dataset_mean_expr_COMBAT_", data_type, "_", expr_type, "_", ncol(expr_combat), "_datasets_", nrow(expr), "_genes.csv"))
  
}

combat_class_fxn <- function(datinfo, expr_paths, data_type, expr_type, gene_list=NULL){
  
  for(i in 1:length(expr_paths)){
    
    cell_class <- strsplit(
      strsplit(expr_paths[i], "/", fixed=T)[[1]][5], "_"
    )[[1]][1]
    
    cat("\n")
    print(cell_class)
    
    expr <- fread(expr_paths[i], data.table=F)
    idx <- grep(paste(datinfo$Dataset, collapse="|"), colnames(expr))
    expr <- expr[,c(1, idx)]
    
    labs <- rep(datinfo$Dataset, lengths(lapply(
      datinfo$Dataset, function(x) grep(x, colnames(expr))
    )))
    datinfo1 <- datinfo[match(labs, datinfo$Dataset),]
    
    if(!identical(labs, datinfo1$Dataset)){
      stop("!identical(labs, datinfo1$Dataset)")
    }
    
    ## Remove genes that are NA in ANY dataset:
    
    expr <- expr[apply(expr[,-c(1)], 1, function(x) sum(is.na(x)))==0,]
    
    ## Remove zero variance genes:
    
    expr <- expr[apply(expr[,-c(1)], 1, function(x) sum(x))>0,]
    
    if(!is.null(gene_list)){
      
      expr <- expr[is.element(expr$SYMBOL, gene_list),]
      
    }
    
    batch_vec <- as.factor(paste(datinfo1$Study, datinfo1$Platform))
    
    expr_combat <- sva::ComBat(expr[,-c(1)], batch=batch_vec)
    
    expr <- data.frame(SYMBOL=expr[,c(1)], expr_combat)

    fwrite(expr, file=paste0("data/", data_type, "/", expr_type, "/", cell_class, "_mean_expr_COMBAT_", data_type, "_", expr_type, "_", ncol(expr_combat), "_datasets_", nrow(expr), "_genes.csv"))
    
  }
  
}