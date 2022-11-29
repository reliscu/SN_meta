setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/4_map_identifiers")

library(dplyr)
library(data.table)

source("/home/rebecca/code/misc/map_identifiers/map_identifiers_function.R")

datinfo <- read.csv("../3_init_datinfo/datinfo_SCSN_meta_analysis_init.csv") %>% na_if("")

## Check what type of identifiers each dataset has:

for(i in 1:nrow(datinfo)){
  genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
  cat("\n")
  print(datinfo$Dataset[i])
  print(genes[1,])
  if(sum(grepl(".AS", genes[,1], fixed=T))>0){
    stop()
  }
}

## Separate datasets into those with ENSEMBL / ENTREZ / SYMBOL identifiers:

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

## Look at datasets with symbol identifiers:

idx <- which(is.element(datinfo$Dataset, sym_list))

for(i in idx){
  genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
  print(i)
  print(datinfo$Dataset[i])
  print(head(genes, 1))
  cat("\n")
}

## How many datasets with only symbol identfiers have duplicate identifiers?

for(i in idx){
  genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
  if(sum(duplicated(genes[,1]))>0){
    print(i)
    print(datinfo$Dataset[i])
  }
}

## Only 1.

## What are these duplicate identifiers?

for(i in idx){
  genes <- fread(datinfo$Author_Genes[i], data.table=F, header=F)
  if(sum(duplicated(genes[,1]))>0){
    print(genes[,1][duplicated(genes[,1])])
    print(datinfo$Dataset[i])
  }
}

## Mostly related to sex chromosomes. Let's remove rows associated with these genes. (Doing this in the prep data file for this dataset).


