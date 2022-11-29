setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/4_map_identifiers")

library(dplyr)
library(biomaRt)

source("/home/rebecca/code/misc/map_identifiers/map_identifiers_function.R")

tables_path <- "/home/rebecca/omicon/mapping_tables"

datinfo <- read.csv("../3_init_datinfo/datinfo_SCSN_meta_analysis_init.csv") %>% na_if("")

############################################# Separate datasets by ENSEMBL / ENTREZ / SYMBOL identifiers ############################################# 

ens_list <- c()
entr_list <- c()
sym_list <- c()

for(i in 1:nrow(datinfo)){
  genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
  if(length(grep("ENS", genes[1,]))>0){
    ens_list <- c(ens_list, datinfo$Dataset[i]) 
  } else if(ncol(genes)>1){
    entr_list <- c(entr_list, datinfo$Dataset[i])
  } else {
    sym_list <- c(sym_list, datinfo$Dataset[i])
  }
}

############################################# Look at situations where mapping produces duplicate symbol identifiers ############################################# 

# idx <- which(is.element(datinfo$Dataset, sym_list))
# 
# i <- 3
# 
# for(i in idx){
#   genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
#   genes_mapped <- mapAlias2Symbol(features=genes, unique_id_col=1, tables_path, keep_all=F, fill_NAs=T)
#   
#   temp <- genes_mapped[is.element(genes_mapped$SYMBOL, genes_mapped$SYMBOL[duplicated(genes_mapped$SYMBOL)]),]
#   temp <- temp %>% dplyr::arrange(SYMBOL)
#   dup_ids <- unique(temp$SYMBOL)
#   for(j in 1:length(dup_ids)){
#     temp1 <- temp[is.element(temp$SYMBOL, dup_ids[j]),]
#     matches <- sum(temp1$UNIQUE.ID==temp1$SYMBOL)
#     if(matches==0){
#       cat("\n")
#       print(temp1)
#     }
#   }
# 
#   expr <- fread(datinfo$Author_Counts[i], data.table=F)
#   
#   gene <- temp$SYMBOL[5]
#   print(gene)
#   idx1 <- which(genes_mapped$SYMBOL==gene)
#   cor(t(expr[idx1,]))
# 
# }

############################################# Map symbol identifiers ############################################# 

idx <- which(is.element(datinfo$Dataset, sym_list))

tbl <- fread("/home/rebecca/omicon/mapping_tables/homo_sapiens_mapping_tables/2022-11-02/SYMBOL/Homo_sapiens_SYMBOL_ALIAS.csv", data.table=F)

for(i in idx){
  genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
  genes_mapped <- mapAlias2Symbol(features=genes, unique_id_col=1, tables_path, keep_all=F, fill_NAs=T)
  ## If a mapping is duplicated and the unique identifier is an official gene symbol, keep original unique identifier:
  # dup_ids <- genes_mapped$SYMBOL[duplicated(genes_mapped$SYMBOL)]
  # idx1 <- which(is.element(toupper(genes_mapped$UNIQUE.ID), tbl$SYMBOL)&is.element(genes_mapped$SYMBOL, dup_ids))
  # genes_mapped$SYMBOL[idx1] <- genes_mapped$UNIQUE.ID[idx1]
  ## Save:
  file_path <- paste0(getwd(), "/data/", datinfo$Dataset[i], "_genes_mapped.csv")
  datinfo$Author_Genes_Mapped[i] <- file_path
  fwrite(genes_mapped, file=file_path)
}

############################################# Map entrez identifiers #############################################

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="https://asia.ensembl.org")

idx <- which(is.element(datinfo$Dataset, entr_list))

# ## What percent of identifiers are mapped with Annotation Db?
# 
# 
# for(i in idx){
#   genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
#   genes_mapped <- map2Any(features=genes, unique_id="ENTREZID", map_to="SYMBOL", unique_id_col=1, platform=NULL, tables_path, keep_all=F)
#   print(sum(!is.na(genes_mapped$SYMBOL))/nrow(genes_mapped))
# }
# 
# # [1] 0.9189157
# 
# ## What percent are mapped with biomaRt?
# 
# for(i in idx){
#   genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
#   genes_mapped <- getBM(attributes=c("entrezgene_id", "hgnc_symbol"), filters="entrezgene_id", values=genes[,1], mart=mart)
#   genes_mapped <- genes_mapped[genes_mapped[,2]!="",]
#   print(nrow(genes_mapped)/nrow(genes))
# }
# 
# # [1] 0.5105507

## Annotation Db fares far beter

for(i in idx){
  genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
  genes_mapped <- map2Any(features=genes, unique_id="ENTREZID", map_to="SYMBOL", unique_id_col=1, platform=NULL, tables_path, keep_all=F)
  genes_mapped <- genes_mapped[,c(1,2)]
  file_path <- paste0(getwd(), "/data/", datinfo$Dataset[i], "_genes_mapped.csv")
  datinfo$Author_Genes_Mapped[i] <- file_path
  fwrite(genes_mapped, file=file_path)
}

############################################# Map ensembl identifiers ############################################# 

idx <- which(is.element(datinfo$Dataset, ens_list))

# ## What percent of identifiers are mapped with Annotation Db?
# 
# 
# for(i in idx){
#   genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
#   if(length(grep(".", genes[,1], fixed=T))>0){
#     genes[,1] <- sapply(strsplit(genes[,1], ".", fixed=T), "[", 1)
#   }
#   genes_mapped <- map2Any(features=genes, unique_id="ENSEMBL", map_to="SYMBOL", unique_id_col=1, platform=NULL, tables_path, keep_all_mappings=F)
#   print(sum(!is.na(genes_mapped$SYMBOL))/nrow(genes_mapped))
# }
# 
# # [1] 0.6911023
# # ...
# 
# ## What percent are mapped with biomaRt?
# 
# for(i in idx){
#   genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
#   if(length(grep(".", genes[,1], fixed=T))>0){
#     genes[,1] <- sapply(strsplit(genes[,1], ".", fixed=T), "[", 1)
#   }
#   genes_mapped <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes[,1], mart=mart)
#   ## Make sure order of genes matches original order:
#   genes_mapped <- genes_mapped %>% na_if("")
#   genes_mapped <- genes_mapped[match(genes[,1], genes_mapped[,1]),]
#   genes_mapped[is.na(genes_mapped[,1]),1] <- genes[is.na(genes_mapped[,1]),1]
#   if(!identical(genes_mapped[,1], genes[,1])){
#     stop("!identical(genes_mapped[,1], genes[,1])")
#   }
#   ## Now map to symbols:
#   genes_mapped <- mapAlias2Symbol(features=genes_mapped, unique_id_col=2, tables_path, keep_all_mappings=F, fill_NAs=T)
#   genes_mapped <- genes_mapped[,c(3,2)]
#   print(sum(!is.na(genes_mapped[,2]))/nrow(genes))
# }
# 
# # [1] 0.7321778
# # ...

## biomaRt fares beter

for(i in idx){
  genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
  if(length(grep(".", genes[,1], fixed=T))>0){
    genes[,1] <- sapply(strsplit(genes[,1], ".", fixed=T), "[", 1)
  }
  genes_mapped <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes[,1], mart=mart) %>% na_if("")
  ## Make sure order of genes matches original order:
  genes_mapped <- genes_mapped[match(genes[,1], genes_mapped[,1]),]
  genes_mapped[is.na(genes_mapped[,1]),1] <- genes[is.na(genes_mapped[,1]),1]
  if(!identical(genes_mapped[,1], genes[,1])){
    stop("!identical(genes_mapped[,1], genes[,1])")
  }
  ## Now map to symbols:
  genes_mapped <- mapAlias2Symbol(features=genes_mapped, unique_id_col=2, tables_path, keep_all=F, fill_NAs=T)
  # ## If a symbol is duplicated and the original unique identifier is an official gene symbol, keep original unique identifier:
  # dup_ids <- genes_mapped$SYMBOL[duplicated(genes_mapped$SYMBOL)]
  # idx1 <- which(is.element(toupper(genes_mapped$UNIQUE.ID), tbl$SYMBOL)&is.element(genes_mapped$SYMBOL, dup_ids))
  # genes_mapped$SYMBOL[idx1] <- genes_mapped$UNIQUE.ID[idx1]
  genes_mapped <- genes_mapped[,c(3,2)]
  colnames(genes_mapped) <- c("UNIQUE.ID", "SYMBOL")
  ## Save:
  file_path <- paste0(getwd(), "/data/", datinfo$Dataset[i], "_genes_mapped.csv")
  datinfo$Author_Genes_Mapped[i] <- file_path
  fwrite(genes_mapped, file=file_path)
}

############################################# Make sure order of mapped identifiers matches original order ############################################# 

for(i in 1:nrow(datinfo)){
  genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
  if(length(intersect(grep(".", genes[,1], fixed=T), grep("ENS", genes[,1])))>0){
    genes[,1] <- sapply(strsplit(genes[,1], ".", fixed=T), "[", 1)
  }
  genes_mapped <- fread(datinfo$Author_Genes_Mapped[i], data.table=F)
  if(!identical(genes[,1], genes_mapped$UNIQUE.ID)){
    cat("")
    print(i)
    print(datinfo$Dataset[i])
  }
}

datinfo <- fwrite(datinfo, file="datinfo_SCSN_meta_analysis_mapped_genes.csv")
