library(dplyr)
library(Matrix)
library(data.table)

mega_matrix <- function(datinfo, data_type, expr_type, cell_classes=c("ASC", "END", "EXC", "INH", "NEU", "MIC", "OG", "OPC", "PER", "VSMC"), gene_list){
  
  cellinfo_list <- vector(mode="list", length=nrow(datinfo))
  expr_list <- vector(mode="list", length=nrow(datinfo))
  
  for(i in 1:nrow(datinfo)){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    expr <- fread(datinfo$Author_Normalized_Counts_QC_Intersection_PC_Genes[i], data.table=F)
    expr <- expr[match(gene_list, expr[,c(1)]),]
    expr <- expr[,-c(1)]
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations_Class[i], data.table=F)
    cellinfo$Cell_ID <- make.names(cellinfo$Cell_ID)
    cellinfo <- cellinfo[match(colnames(expr), cellinfo$Cell_ID),]

    if(!identical(colnames(expr), cellinfo$Cell_ID)){
      stop("!identical(colnames(expr), cellinfo$Cell_ID)")
    }
    
    ## Subset to cell classes of interest:
    
    idx <- which(is.element(cellinfo$MO_Cell_Class, cell_classes))
    expr <- expr[,idx]
    cellinfo <- cellinfo[idx,]
    
    cellinfo <- cellinfo %>% 
      dplyr::select(
        Cell_ID, Cell_Type, MO_Cell_Class
      ) %>% 
      dplyr::mutate(Dataset=datinfo$Dataset[i]) %>%
      dplyr::relocate(Dataset, .before="Cell_ID")
    
    expr_list[[i]] <- as.matrix(expr)
    cellinfo_list[[i]] <- cellinfo
    
  } ## for(j in 1:length(cell_classes)){
  
  cellinfo <- do.call(rbind, cellinfo_list)
  expr <- do.call(cbind, expr_list)
  rownames(expr) <- gene_list

  fwrite(cellinfo, file=paste0("data/", data_type, "/", expr_type, "/cellinfo_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", n_distinct(datinfo$Dataset), "_datasets_", nrow(expr), "_intersection_genes.csv"))
  
  saveRDS(expr, file=paste0("data/", data_type, "/", expr_type, "/expr_", data_type, "_", expr_type, "_", ncol(expr), "_nuclei_", n_distinct(datinfo$Dataset), "_datasets_", nrow(expr), "_intersection_genes.RDS"))
  
} ## cell_class_fidelity
