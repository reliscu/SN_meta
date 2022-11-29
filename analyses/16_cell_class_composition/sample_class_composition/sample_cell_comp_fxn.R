library(dplyr)
library(data.table)
library(Matrix)
library(future.apply)
library(openblasctl)

openblas_set_num_threads(10)
options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

sample_cell_prop <- function(
  datinfo, 
  cell_classes=c("INH", "END", "EXC", "ASC", "MIC", "OG", "OPC", "PER", "VSMC"), 
  props, ## Vector of desired cell class proportions of the same length as 'cell_classes'
  n_nuclei, ## Total # nuclei to sample
  n_resamples
){
  
  for(i in 1:nrow(datinfo)){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr_path <- datinfo$Author_Counts_QC[i]
    
    if(grepl("mtx", expr_path)){
      
      expr <- readMM(expr_path)
      
    } else {
      
      expr <- fread(expr_path, data.table=F)
      
    }
    
    expr <- as(expr, "matrix")
    
    genes <- fread(datinfo$Author_Genes_Mapped[i], header=F, data.table=F) %>% na_if("")
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
    
    n_nuclei_per_class <- ceiling(n_nuclei*props)
    
    ## Sample each cell class to desired proportion:
 
    expr_list <- future_lapply(1:n_resamples, FUN=function(n){
      
      cell_idx <- unlist(lapply(1:length(cell_classes), function(j){
        
        class_idx <- which(
          is.element(cellinfo$MO_Cell_Class, cell_classes[j])
        )
        replace <- F
        if(length(class_idx)<n_nuclei_per_class[j]){
          replace <- T
        }
        
        return(sample(class_idx, size=n_nuclei_per_class[j], replace=replace))
        
      })) ## class_idx <- unlist(
      
      expr <- expr[,cell_idx]
      colnames(expr) <- cellinfo$Cell_ID[cell_idx]
      
      return(expr)
      
    }, future.seed=T) ## expr_list <- future_lapply(
    
    saveRDS(expr_list, file=paste0("data/", data_type, "/", expr_type, "/", datinfo$Dataset[i], "_sampled_cell_class_composition_", data_type, "_", expr_type, "_", n_nuclei, "_nuclei_", n_resamples, "_resamples.RDS"))
    
  } ## for(i in 1:nrow(datinfo)){
  
} ## sample_cell_prop <- function(