setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/5_QC")

library(dplyr)
library(data.table)
library(Matrix)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

source("/home/rebecca/code/misc/rank_percentile.R")

## Restrict to cortical datasets:

datinfo <- read.csv("../4_map_identifiers/datinfo_SCSN_meta_analysis_mapped_genes.csv") %>% na_if("") %>% dplyr::filter(grepl("C$", Region_Code)) 

qc_crit <- paste0("QC_counts_>.05%tile_genes_expressed_&_<=5%_mito_reads")

na_cts <- c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class", "d5")

## Imposing QC requirements:

# Filter:
## 1) Doublets / droplets / unassigned cell types
## 2) Nuclei in bottom .5 %tile of # genes expressed
## 3) Nuclei with >5% of UMIs derived from mitochondrial genes

future_lapply(1:nrow(datinfo), FUN=function(i){
  
  expr_path <- datinfo$Author_Counts[i]
  
  if(grepl("mtx", expr_path)){
    
    ext <- ".mtx"
    expr <- readMM(expr_path)
    
  } else {
    
    ext <- ".csv"
    expr <- fread(expr_path, data.table=F)
    
  }
  
  barcodes <- fread(datinfo$Author_Barcodes[i], header=F, data.table=F) %>% na_if("")
  genes <- fread(datinfo$Author_Genes_Mapped[i], data.table=F) %>% na_if("")
  
  if(!is.na(datinfo$Author_Barcode_Annotations[i])){
    
    ## Filter doublets / droplets / unassigend cell types:
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F) %>% na_if("")
    cellinfo <- cellinfo[match(barcodes[,1], cellinfo$Cell_ID),]
    cellinfo$Cell_Type[is.na(cellinfo$Cell_Type)] <- "Unassigned"
      
    if(!identical(barcodes[,1], cellinfo$Cell_ID)){
      stop("!identical(barcodes[,1], cellinfo$Cell_ID)")
    }
    
    idx <- which(!grepl(paste(na_cts, collapse="|"), cellinfo$Cell_Type))
    expr <- expr[,idx]
    barcodes <- data.frame(barcodes[idx,])
    
  }

  mito_umis <- apply(expr[grep("^MT-", genes[,2]),], 2, sum)
  n_genes_per_cell <- apply(expr, 2, function(x) sum(x>0))
  
  # Nuclei in bottom .5%tile of # genes expressed:
  
  low_expr_cells <- which(n_genes_per_cell<quantile(n_genes_per_cell, .005))
  
  ## Nuclei with >5% of UMIs derived from mitochondrial genes:
  
  low_qual_cells <- which((mito_umis/colSums(expr))>.05)
  
  rm_idx <- c(unique(low_expr_cells, low_qual_cells))
  
  expr <- expr[,-rm_idx]
  barcodes <- data.frame(barcodes[-rm_idx,])
  
  expr_path <- paste0("data/author_data/", datinfo$Dataset[i], "_expression_", qc_crit, ext)
  
  if(grepl("mtx", ext)){
    writeMM(expr, file=expr_path)
  } else {
    fwrite(expr, file=expr_path)
  }
  
  fwrite(barcodes, file=paste0("data/author_data/", datinfo$Dataset[i], "_barcodes_", qc_crit, ".csv"), col.names=F)
  
}) 

## Add datinfo paths:

datinfo <- read.csv("../4_map_identifiers/datinfo_SCSN_meta_analysis_mapped_genes.csv") %>% na_if("") %>% dplyr::filter(grepl("C$", Region_Code)) 

qc_crit <- paste0("QC_counts_>.05%tile_genes_expressed_&_<=5%_mito_reads")

ext <- sapply(strsplit(datinfo$Author_Counts, ".", fixed=T), function(x) x[length(x)])

datinfo$Author_Counts_QC <- paste0(getwd(), "/data/author_data/", datinfo$Dataset, "_expression_", qc_crit, ".", ext)

datinfo$Author_Barcodes_QC <- paste0(getwd(), "/data/author_data/", datinfo$Dataset, "_barcodes_", qc_crit, ".csv")

fwrite(datinfo, file="datinfo_SCSN_meta_analysis_QC.csv")
