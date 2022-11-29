setwd("/home/rebecca/SCSN_meta_analysis/gene_sets/SCSN_gene_sets")

library(dplyr)

#current_legend <- read.csv("/home/shared/genesets/OurSets/MyGeneSetsLEGENDold5.csv")
max_set_no <- max(as.numeric(gsub("[A-Z]", "", current_legend$SetID)))

# sets <- list.files(path="/home/shared/genesets/OurSets/", pattern="MOSET")
# set_no <- as.numeric(gsub("[A-Z]", "", gsub(".csv", "", sets)))
# rm_idx <- which(set_no>max_set_no)
# file.remove(list.files(path="/home/shared/genesets/OurSets/", pattern="MOSET", full.names=T)[rm_idx])

new_legend <- data.frame(
  SetID=NA,
  SetName=NA,
  Category="Cell type",
  SetSize=NA,
  Species="Homo sapiens",
  PubMed=NA,
  Description=NA
)

#######################################################################################################
############################################# Agarwal CTX #############################################
#######################################################################################################

pmid <- "32826893"
author <- "agarwal_CTX_2020"

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/cortex/marker_genes_cortex_DE_all_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

## All cell types 

j <- max_set_no
for(i in 1:length(marker_list)) {
  new_legend[i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all cortex cell type clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## All excitatory cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/cortex/marker_genes_cortex_DE_excitatory_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_EX_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to cortex excitatory neuron clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## All inhibitory cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/cortex/marker_genes_cortex_DE_inhibitory_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_IN_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to cortex inhibitory neuron clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}


######################################################################################################
############################################# Agarwal SN #############################################
######################################################################################################

author <- "agarwal_SN_2020"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/SN/marker_genes_SN_DE_all_clusters_tableS3.csv") %>% na_if("")
marker_df$Cluster <- gsub("GABAergic neurons", "GABA neurons", marker_df$Cluster)
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all substantia nigra (SN) cell type clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Astrocyte cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/SN/marker_genes_SN_astrocyte_DE_astrocyte_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_ASC_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to substantia nigra astrocyte clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Neuron cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/SN/marker_genes_SN_neuron_DE_neuron_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_NEU_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to substantia nigra neuron clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

######################################################################################################
############################################# Habib 2017 #############################################
######################################################################################################

pmid <- "28846088"
author <- "habib_2017"

## All cell types 

clust_labels <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/cluster_label_mapping_tableS8.csv")
marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/marker_genes_DE_all_clusters_tableS8.csv")
ids2label <- clust_labels$Name[match(marker_df$Cluster.ID, clust_labels$Cell.ID)]
marker_df$Cluster.ID <- ids2label
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster.ID[!is.na(marker_df$Cluster.ID)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts
names(marker_list)[grep("NSC", names(marker_list))] <- "ependymal"

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S8)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## GABAergic cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/marker_genes_GABAergic_DE_GABAergic_clusters_tableS9.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster.ID[!is.na(marker_df$Cluster.ID)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- paste0("GABA", cts)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_GABA_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to GABAergic clusters (table S9)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Excitatory cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/marker_genes_PFC_glutamatergic_DE_glutamatergic_clusters_tableS9.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster.ID[!is.na(marker_df$Cluster.ID)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- paste0("GLUT", cts)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_GLUT_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "PFC cell type genes differentially expressed with respect to glutamatergic clusters (table S9)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#####################################################################################################
############################################# Lake 2018 #############################################
#####################################################################################################

pmid <- "29227469"
author <- "lake_2018"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/marker_genes_DE_all_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Excitatory cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/marker_genes_excitatory_DE_excitatory_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_EX_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to excitatory neuron clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Inhbibitory cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/marker_genes_inhibitory_DE_inhibitory_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_IN_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to inhibitory neuron clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Oligo cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/marker_genes_oligo_DE_oligo_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_OG_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to oligodendrocyte and OPC clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

author <- "lake_2018_CER"

## Cer astrocyte cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/marker_genes_cerebellum_astrocyte_DE_astrocyte_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_ASC_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to cerebellum astrocyte clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Cerebellum cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/marker_genes_cerebellum_DE_cerebellar_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to cerebellum cell type clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Cer OPC cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/marker_genes_cerebellum_OPC_DE_OPC_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_OPC_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to cerebellum OPC clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

###################################################################################################
############################################# Li 2018 #############################################
###################################################################################################

pmid <- "30545854"
author <- "li_2018"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/li_2018/marker_genes_tableS8.csv")
marker_df$Cell.type <- gsub("_adult", "", marker_df$Cell.type)
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cell.type[!is.na(marker_df$Cell.type)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cell.type%in%ct])
names(marker_list) <- cts
names(marker_list)[grep("INNT", names(marker_list), ignore.case=T)] <- "IN"

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S8)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#######################################################################################################
############################################# Mathys 2019 #############################################
#######################################################################################################

pmid <- "31042697"
author <- "mathys_2019"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/mathys_2019/marker_genes_DE_subtype_clusters_tableS6.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$subpopulation[!is.na(marker_df$subpopulation)])
marker_list <- lapply(cts, function(ct) marker_df$gene.name[marker_df$subpopulation%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_subtype_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to subtype clusters (table S6)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#########################################################################################################
############################################# Schirmer 2019 #############################################
#########################################################################################################

pmid <- "31316211"
author <- "schirmer_2019"

## All cell types 

marker_df <- fread("/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/marker_genes_DE_all_clusters_tableS6.csv", data.table=F)
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,1][!is.na(marker_df[,1])])
marker_list <- lapply(cts, function(ct) marker_df[,3][marker_df[,1]%in%ct])
names(marker_list) <- cts
names(marker_list) <- gsub("EN", "EX", names(marker_list))

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S6)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#########################################################################################################
############################################# Tran 2020 ACC #############################################
#########################################################################################################

pmid <- "34582785"
author <- "tran_ACC_2020"

## All cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_ACC_top40_tableS5.csv")
marker_df <- marker_df[,grepl("_1vAll", colnames(marker_df))]
colnames(marker_df) <- gsub("_1vAll", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to anterior cingulate cortex (ACC) cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## All cell types pairwise

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_ACC_top40_tableS5.csv")
marker_df <- marker_df[,grepl("_pw", colnames(marker_df))]
colnames(marker_df) <- gsub("_pw", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_pairwise", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed in pairwise comparison of anterior cingulate cortex (ACC) cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#########################################################################################################
############################################# Tran 2020 AMY #############################################
#########################################################################################################

pmid <- "34582785"
author <- "tran_AMY_2020"

## All cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_amygdala_top40_tableS5.csv") %>% na_if("")
marker_df <- marker_df[,grepl("_1vAll", colnames(marker_df))]
colnames(marker_df) <- gsub("_1vAll", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to amygdala cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## All cell types pairwise

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_amygdala_top40_tableS5.csv")
marker_df <- marker_df[,grepl("_pw", colnames(marker_df))]
colnames(marker_df) <- gsub("_pw", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_pairwise", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed in pairwise comparison of amygdala cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

###########################################################################################################
############################################# Tran 2020 DLPFC #############################################
###########################################################################################################

pmid <- "34582785"
author <- "tran_DLPFC_2020"

## All cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_DLPFC_top40_tableS5.csv")
marker_df <- marker_df[,grepl("_1vAll", colnames(marker_df))]
colnames(marker_df) <- gsub("_1vAll", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to dorsolateral prefrontal cortex (DLPFC) cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## All cell types pairwise

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_DLPFC_top40_tableS5.csv")
marker_df <- marker_df[,grepl("_pw", colnames(marker_df))]
colnames(marker_df) <- gsub("_pw", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_pairwise", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed in pairwise comparison of dorsolateral prefrontal cortex (DLPFC) cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#########################################################################################################
############################################# Tran 2020 HIP #############################################
#########################################################################################################

pmid <- "34582785"
author <- "tran_HIP_2020"

## All cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_hippocampus_top40_tableS5.csv")
marker_df <- marker_df[,grepl("_1vAll", colnames(marker_df))]
colnames(marker_df) <- gsub("_1vAll", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to hippocampus cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## All cell types pairwise

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_hippocampus_top40_tableS5.csv")
marker_df <- marker_df[,grepl("_pw", colnames(marker_df))]
colnames(marker_df) <- gsub("_pw", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_pairwise", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed in pairwise comparison of hippocampus cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#########################################################################################################
############################################# Tran 2020 NAc #############################################
#########################################################################################################

pmid <- "34582785"
author <- "tran_NAC_2020"

## All cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_NAc_top40_tableS5.csv")
marker_df <- marker_df[,grepl("_1vAll", colnames(marker_df))]
colnames(marker_df) <- gsub("_1vAll", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to nucleus accumbens (NAc) cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## All cell types pairwise

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/marker_genes_NAc_top40_tableS5.csv")
marker_df <- marker_df[,grepl("_pw", colnames(marker_df))]
colnames(marker_df) <- gsub("_pw", "", colnames(marker_df))
marker_df <- reshape2::melt(as.matrix(marker_df))
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df[,2][!is.na(marker_df[,3])])
cts <- cts[!is.na(cts)]
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,3][marker_df[,2]%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_pairwise", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed in pairwise comparison of nucleus accumbens (NAc) cell type clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

##########################################################################################################
############################################# Velmeshev 2019 #############################################
##########################################################################################################

pmid <- "31097668"
author <- "velmeshev_2019"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/marker_genes_DE_all_clusters_tableS3.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cell.type[!is.na(marker_df$Cell.type)])
marker_list <- lapply(cts, function(ct) marker_df$Gene.name[marker_df$Cell.type%in%ct])
names(marker_list) <- cts
names(marker_list)[grep("^L", names(marker_list))] <- paste0(names(marker_list)[grep("^L", names(marker_list))], "-Neu")

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}


######################################################################################################
############################################# Yang COVID #############################################
######################################################################################################

pmid <- "34153974"

## PFC cell types 

author <- "yang_PFC_2021"
marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/marker_genes_PFC_DE_all_clusters_tableS2.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type differentially expressed with respect to prefrontal cortex (PFC) cell type clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Choroid cell types 

author <- "yang_CP_2021"
marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/marker_genes_choroid_DE_all_clusters_tableS4.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- cts
names(marker_list) <- gsub("-HI", "", names(marker_list), ignore.case=T)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to choroid plexus (CP) cell type clusters (table S4)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}


###################################################################################################
############################################# Yang AD #############################################
###################################################################################################

pmid <- "35165441"
author <- "yang_2021"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/marker_genes_DE_all_clusters_tableS2.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cell.Type[!is.na(marker_df$Cell.Type)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cell.Type%in%ct])
names(marker_list) <- cts
names(marker_list) <- gsub(",", "", names(marker_list), fixed=T)
names(marker_list) <- gsub("/ ", "/", names(marker_list), fixed=T)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Myeloid cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/marker_genes_myeloid_DE_myeloid_clusters_tableS5.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$cluster%in%ct])
names(marker_list) <- cts
names(marker_list)[length(marker_list)] <- "pericyte meningeal macrophage"

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_myeloid_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to myeloid clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## BEC cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/marker_genes_BEC_DE_BEC_clusters_tableS5.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cell.subtype[!is.na(marker_df$Cell.subtype)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cell.subtype%in%ct])
names(marker_list) <- cts
names(marker_list)[length(marker_list)] <- "proteostatic"

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_BEC_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to blood vascular endothelial cell (BEC) clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Fibroblast cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/marker_genes_fibroblast_DE_fibroblast_clusters_tableS5.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cell.subtype[!is.na(marker_df$Cell.subtype)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cell.subtype%in%ct])
names(marker_list) <- cts
names(marker_list)[length(marker_list)] <- "arachnoid_pial_pericyte"

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_fibro_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to fibroblast clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Mural cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/marker_genes_mural_DE_mural_clusters_tableS5.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cell.subtype[!is.na(marker_df$Cell.subtype)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cell.subtype%in%ct])
names(marker_list) <- cts
names(marker_list) <- gsub(" ", "", names(marker_list))

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_mural_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to mural clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#####################################################################################################
############################################# Nagy 2020 #############################################
#####################################################################################################

pmid <- "32341540"
author <- "nagy_2020"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/nagy_2020/marker_genes_table_DE_all_clusters_S6.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Cluster%in%ct])
names(marker_list) <- cts
names(marker_list) <- gsub(",", "", names(marker_list), fixed=T)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

########################################################################################################
############################################# Grubman 2019 #############################################
########################################################################################################

pmid <- "31768052"
author <- "grubman_2019"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/grubman_2019/marker_genes_DE_all_clusters_tableS2.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$group[!is.na(marker_df$group)])
marker_list <- lapply(cts, function(ct) marker_df$geneID[marker_df$group%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

########################################################################################################
############################################# Absinta 2021 #############################################
########################################################################################################

pmid <- "34497421"
author <- "absinta_2021"

## All cell types 

clust_labels <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/cluster_label_mapping_fig1.csv") 
marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/marker_genes_DE_all_clusters_tableS3.csv") %>% na_if("")
ids2label <- clust_labels$Cell_type[match(marker_df$cluster, clust_labels$Cluster)]
marker_df$cluster <- ids2label
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Astrocyte cell types 

clust_labels <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/cluster_label_mapping_astrocyte_supp_info.csv")
marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/marker_genes_astrocyte_DE_astrocyte_clusters_tableS7.csv") %>% na_if("")
ids2label <- clust_labels$Cell_type[match(marker_df$cluster, clust_labels$Cluster)]
marker_df$cluster <- ids2label
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- gsub(" ", "_", cts)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_ASC_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to astrocyte clusters (table S7)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Immune cell types 

clust_labels <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/cluster_label_mapping_immune_supp_info.csv")
marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/marker_genes_immune_DE_immune_clusters_tableS4.csv") %>% na_if("")
ids2label <- clust_labels$Cell_type[match(marker_df$cluster, clust_labels$Cluster)]
marker_df$cluster <- ids2label
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- gsub(" ", "_", cts)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_MG_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to immune clusters (table S4)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Oligo cell types 

clust_labels <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/cluster_label_mapping_oligodendrocyte_supp_info.csv")
marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/marker_genes_oligo_DE_oligo_clusters_tableS5.csv") %>% na_if("")
ids2label <- clust_labels$Cell_type[match(marker_df$cluster, clust_labels$Cluster)] 
marker_df$cluster <- ids2label
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- gsub(" ", "_", cts)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_OG_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to oligodendrocyte lineage clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Vascular cell types 

clust_labels <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/cluster_label_mapping_vasc_supp_info.csv")
marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/marker_genes_vascular_DE_vascular_clusters_tableS13.csv") %>% na_if("")
ids2label <- clust_labels$Cell_type[match(marker_df$cluster, clust_labels$Cluster)]
marker_df$cluster <- ids2label
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- gsub(" ", "_", cts)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_VASC_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to vascular clusters (table S13)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

####################################################################################################
############################################# Lau 2020 #############################################
####################################################################################################

pmid <- "32989152"
author <- "lau_2020"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/lau_2020/marker_genes_DE_all_clusters_tableS2.csv") %>% na_if("")
cts <- unique(marker_df$Target_cell_type[!is.na(marker_df$Target_cell_type)])
marker_list <- lapply(cts, function(ct) marker_df$Gene[marker_df$Target_cell_type%in%ct])
names(marker_list) <- cts
names(marker_list)[grep("Mic", names(marker_list))] <- "Microglia"

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}
####################################################################################################
############################################# Morabito 2021 ########################################
####################################################################################################

pmid <- "34239132"
author <- "morabito_2021"

## All cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/morabito_2021/marker_genes_DE_all_clusters_dataS1c.csv") %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- cts
names(marker_list) <- gsub(".", "_", names(marker_list), fixed=T)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (data S1c)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Subtype clusters

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/morabito_2021/marker_genes_DE_subtype_clusters_dataS1d.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- cts
names(marker_list) <- gsub("ODC", "oligo", names(marker_list))
names(marker_list) <- gsub("MG", "microglia", names(marker_list))
names(marker_list) <- gsub("ASC", "astrocyte", names(marker_list))
names(marker_list) <- gsub(".", "_", names(marker_list), fixed=T)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_subtype_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to subtype clusters (data S1d)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

###############################################################################################
############################################# Luo 2019 ########################################
###############################################################################################

pmid <- "35419551"
author <- "luo_2019"

## Neuron cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/marker_genes_mCH_neuron_DE_neuron_clusters_tableS8.csv") %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene.Name[marker_df$Cluster%in%ct])
names(marker_list) <- cts
names(marker_list)[grep("-$", names(marker_list))] <- gsub("_-", "", names(marker_list)[grep("-$", names(marker_list))])

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_neuron_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to neuron clusters (table S8)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Non-neuron cell types 

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/marker_genes_mCH_non-neuron_DE_non-neuronal_clusters_tableS8.csv") %>% na_if("")
cts <- unique(marker_df$Cluster[!is.na(marker_df$Cluster)])
marker_list <- lapply(cts, function(ct) marker_df$Gene.Name[marker_df$Cluster%in%ct])
names(marker_list) <- cts
names(marker_list) <- gsub("NonN_", "", names(marker_list))

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_glial_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to non-neuronal clusters (table S8)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#################################################################################################
############################################# Jakel 2019 ########################################
#################################################################################################

pmid <- "30747918"
author <- "jakel_2019"

## All cell types

marker_files <- list.files(path="/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/", pattern="tableS3", full.names=T)
ct_labels <- list.files(path="/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/", pattern="tableS3")
ct_labels <- sapply(strsplit(ct_labels, "marker_genes_"), "[", 2)
ct_labels <- sapply(strsplit(ct_labels, "_DE"), "[", 1)
for(i in 1:length(marker_files)) {
  marker_df1 <- read.csv(marker_files[i])
  marker_df1 <- marker_df1[,!grepl("^X", colnames(marker_df1))]
  marker_df1$ct <- ct_labels[i]
  if(i==1) {
    marker_df <- marker_df1
  } else {
    marker_df <- rbind(marker_df, marker_df1)
  }
}
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$ct[!is.na(marker_df$ct)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$ct%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S3)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Oligo subtypes

marker_files <- list.files(path="/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/", pattern="tableS4", full.names=T)
ct_labels <- list.files(path="/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/", pattern="tableS4")
marker_files <- marker_files[-c(2)]
ct_labels <- ct_labels[-c(2)]
ct_labels <- sapply(strsplit(ct_labels, "marker_genes_"), "[", 2)
ct_labels <- sapply(strsplit(ct_labels, "_DE"), "[", 1)
for(i in 1:length(marker_files)) {
  marker_df1 <- read.csv(marker_files[i])
  marker_df1 <- marker_df1[,!grepl("^X", colnames(marker_df1))]
  marker_df1$ct <- ct_labels[i]
  if(i==1) {
    marker_df <- marker_df1
  } else {
    marker_df <- rbind(marker_df, marker_df1)
  }
}
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$ct[!is.na(marker_df$ct)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$ct%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_OG_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to oligdendrocyte lineage clusters (table S4)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Oligo subtype markers final

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/marker_genes_final_oligo_DE_oligo_clusters_tableS4.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(colnames(marker_df))
marker_list <- lapply(cts, function(ct) {
  temp <- marker_df[,colnames(marker_df)%in%ct]
  temp[!is.na(temp)]
})
names(marker_list) <- cts
names(marker_list) <- gsub(".", "_", names(marker_list), fixed=T)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_OG_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Final set of oligdendrocyte lineage cell type marker genes differentially expressed with respect to oligdendrocyte lineage clusters (table S4)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

##################################################################################################
############################################# Kamath 2021 ########################################
##################################################################################################

pmid <- "35513515"
author <- "kamath_SN_2021"

## DA neuron cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/marker_genes_DA_neurons_DE_DA_neurons_tableS8.csv") %>% na_if("")
cts <- unique(marker_df$DA_subtype[!is.na(marker_df$DA_subtype)])
marker_list <- lapply(cts, function(ct) marker_df$primerid[marker_df$DA_subtype%in%ct])
names(marker_list) <- cts
names(marker_list) <- paste0("DAN_", names(marker_list))

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_DAN_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to substantia nigra (SN) dopaminergic neuron (DAN) clusters (table S8)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

##################################################################################################
############################################# Bakken 2019 ########################################
##################################################################################################

pmid <- "34616062"
author <- "bakken_2019"

## All cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/marker_genes_tableS4.csv") %>% na_if("")
marker_df <- marker_df[,-c(2:7)]
cts <- unique(marker_df$clusterName[!is.na(marker_df$clusterName)])
marker_list <- lapply(cts, function(ct) {
  markers <- as.character(marker_df[marker_df$clusterName%in%ct,2:ncol(marker_df)])
  markers[!is.na(markers)]
})
names(marker_list) <- cts
names(marker_list)[grep("^L[0-9]", names(marker_list))] <- paste(names(marker_list)[grep("^L[0-9]", names(marker_list))], "neuron")
gaba_list <- c("LAMP5", "PVALB", "SNCG", "SST CHODL", "SST", "VIP")
names(marker_list)[names(marker_list)%in%gaba_list] <- paste(names(marker_list)[names(marker_list)%in%gaba_list], "GABAergic")

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type marker genes (table S8)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## GABAergic cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/marker_genes_GABAergic_DE_GABAergic_tableS7.csv") %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- cts
names(marker_list) <- paste(names(marker_list), "GABAergic")

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_GABA_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect GABAergic clusters (table S7)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Glutamatergic cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/marker_genes_glutamatergic_DE_glutamatergic_tableS10.csv") %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- cts
names(marker_list) <- paste(names(marker_list), "glutamatergic")

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_GLUT_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect glutamatergic clusters (table S10)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#################################################################################################
############################################# Hodge 2018 ########################################
#################################################################################################

pmid <- "31435019"
author <- "hodge_2018"

## Single vs. all

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv") %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$single_markers_vs_all[marker_df$cluster%in%ct]
  unlist(strsplit(markers, ", ", fixed=T))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "single_vs_all", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Single vs. level 1

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv") %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$single_markers_vs_level1[marker_df$cluster%in%ct]
  unlist(strsplit(markers, ", ", fixed=T))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "single_vs_L1", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect level 1 clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Single vs. level 2

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv")
marker_df <- marker_df %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$single_markers_vs_level2[marker_df$cluster%in%ct]
  unlist(strsplit(markers, ", ", fixed=T))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "single_vs_L2", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to level 2 clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Single vs. level 3

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv")
marker_df <- marker_df %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$single_markers_vs_level3[marker_df$cluster%in%ct]
  unlist(strsplit(markers, ", ", fixed=T))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "single_vs_L3", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to level 3 clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Single vs. level 4

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv")
marker_df <- marker_df %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$single_markers_vs_level4[marker_df$cluster%in%ct]
  unlist(strsplit(markers, ", ", fixed=T))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "single_vs_L4", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to level 4 clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Level 4 vs. all

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv")
marker_df <- marker_df %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$level4_markers_vs_all[marker_df$cluster%in%ct]
  unlist(strsplit(markers, ", ", fixed=T))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "L4_vs_all", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Level 4 cell type genes differentially expressed with respect to all clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Combo vs. all

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv")
marker_df <- marker_df %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$combo_markers_vs_all[marker_df$cluster%in%ct]
  markers <- unlist(strsplit(markers, "|", fixed=T))
  unlist(lapply(markers, trimws))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "combo_vs_all", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Combinatorial cell type genes differentially expressed with respect to all clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Combo vs. level 1

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv")
marker_df <- marker_df %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$combo_markers_vs_level1[marker_df$cluster%in%ct]
  markers <- unlist(strsplit(markers, "|", fixed=T))
  unlist(lapply(markers, trimws))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "combo_vs_L1", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Combinatorial cell type genes differentially expressed with respect to level 1 clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Combo vs. level 2

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv")
marker_df <- marker_df %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$combo_markers_vs_level2[marker_df$cluster%in%ct]
  markers <- unlist(strsplit(markers, "|", fixed=T))
  unlist(lapply(markers, trimws))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "combo_vs_L2", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Combinatorial cell type genes differentially expressed with respect to level 2 clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Combo vs. level 3

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/marker_genes_tableS2.csv")
marker_df <- marker_df %>% na_if("") %>% na_if("none")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) {
  markers <- marker_df$combo_markers_vs_level3[marker_df$cluster%in%ct]
  markers <- unlist(strsplit(markers, "|", fixed=T))
  unlist(lapply(markers, trimws))
})
names(marker_list) <- cts
n_NAs <- unlist(lapply(marker_list, function(x) sum(is.na(x))))
marker_list <- marker_list[n_NAs==0]

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "combo_vs_L3", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Combinatorial cell type genes differentially expressed with respect to level 3 clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

#########################################################################################################
############################################# Ayhan 2021 #############################################
#########################################################################################################

pmid <- "34051145"
author <- "ayhan_2021"

## All cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/marker_genes_DE_all_clusters_tableS2.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to all clusters (table S2)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Neuronal cell types

clust_labels <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/marker_genes_neurons_label_mapping.csv")
marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/marker_genes_neuron_DE_neuron_clusters_tableS4.csv")
ids2label <- clust_labels$Cluster.Name[match(marker_df$cluster, clust_labels$Cluster.number)]
marker_df$cluster <- ids2label
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- cts

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_NEU_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to neuronal clusters (table S4)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

## Inhibitory cell types

marker_df <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/marker_genes_inhibitory_DE_inhibitory_clusters_tableS5.csv")
marker_df <- marker_df %>% na_if("")
cts <- unique(marker_df$cluster[!is.na(marker_df$cluster)])
marker_list <- lapply(cts, function(ct) marker_df$gene[marker_df$cluster%in%ct])
names(marker_list) <- paste0("Inh.GAD1.", cts)

j <- max(as.numeric(gsub("[A-Z]", "", new_legend$SetID)))
idx <- nrow(new_legend)
for(i in 1:length(marker_list)) {
  new_legend[idx+i,] <- c(
    paste0("MOSET", j+i),
    toupper(paste(author, names(marker_list)[i], "DE_IN_GAD1_clusters", sep="_")),
    "Cell type",
    length(marker_list[[i]]),
    "Homo sapiens",
    pmid,
    "Cell type genes differentially expressed with respect to GAD1+ inhibitory neuron clusters (table S5)"
  )
  fwrite(data.frame(marker_list[[i]]), file=paste0("MOSET", j+i, ".csv"), col.names=F)
}

############################################# Clean up & save #############################################

new_legend$SetName <- gsub(" ", "_", new_legend$SetName, fixed=T)
new_legend$SetName <- gsub("/", "_", new_legend$SetName, fixed=T)
new_legend$SetName <- gsub(".", "_", new_legend$SetName, fixed=T)
new_legend$SetName <- gsub("-", "_", new_legend$SetName, fixed=T)
new_legend$SetName <- gsub("TCELL", "T_CELL", new_legend$SetName, fixed=T)

# fwrite(new_legend, file="MyGeneSetsLEGEND.csv")
# 
# new_legend <- rbind(current_legend, new_legend)
# 
# fwrite(current_legend, file="/home/shared/genesets/OurSets/MyGeneSetsLEGENDold5.csv")
# fwrite(new_legend, file="/home/shared/genesets/OurSets/MyGeneSetsLEGEND.csv")
# 
