library(dplyr)
library(data.table)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

downsample_UMIs <- function(datinfo, data_type, expr_type, n_umis_list){
  
  for(i in 1:nrow(datinfo)){

    cat("\n")
    print(datinfo$Dataset[i])
    
    expr <- fread(datinfo$Author_Counts_QC_PC_Genes[i], data.table=F)
    genes <- expr[,c(1)]
    expr <- expr[,-c(1)]
    
    downsamp_list <- lapply(1:length(n_umis_list), function(j){
    
      print(paste("No. UMIs:", n_umis_list[j]))
       
      expr_list <- future_lapply(1:ncol(expr), FUN=function(i){
        
        ## Create list of "transcripts" per gene:

        replace <- F
        if(sum(expr[,i])<n_umis_list[j]){
          replace <- T
        }
        
        downsample <- table(sample(rep(genes, times=expr[,i]), size=n_umis_list[j], replace=replace))
        downsample <- downsample[match(genes, names(downsample))]
        
        return(data.frame(as.numeric(downsample)))
  
      }, future.seed=T) ## expr_list <- future_lapply(
      
      expr1 <- do.call(cbind, expr_list)
      expr1[is.na(expr1)] <- 0
      colnames(expr1) <- colnames(expr)
      
      return(data.frame(SYMBOL=genes, expr1))
      
    })
    
    names(downsamp_list) <- as.character(n_umis_list)
    
    saveRDS(downsamp_list, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_UMIs_downsampled_", data_type, "_", expr_type, ".RDS"))
    
  } ## for(i in 1:nrow(datinfo)){
  
}

