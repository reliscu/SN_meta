setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/6_gene_biotype")

library(biomaRt)  
library(GenomicFeatures)
library(AnnotationHub)
library(splicejam)

## Get union of all genes:

gene_list <- lapply(1:length(datinfo), function(i){
  temp <- fread(datinfo$Author_Genes_Mapped[i], data.table=F)
  return(temp$SYMBOL)
})

union_genes <- Reduce(union, gene_list)
length(union_genes)
# [1] 101273

## Get biotype of genes with annotation db:
  
biotype <- AnnotationDbi::select(EnsDb.Hsapiens.v86,  keys=union_genes, columns=c("SYMBOL", "TXBIOTYPE"), keytype="SYMBOL") 
biotype <- biotype %>% dplyr::group_by(SYMBOL) %>% dplyr::summarise(Type=paste(unique(TXBIOTYPE), collapse=" | "))
biotype_df <- data.frame(biotype)
colnames(biotype_df) <- c("Gene", "Type")

missing1 <- union_genes[!is.element(union_genes, biotype_df$Gene)]

length(missing1)
# [1] 46381

## Get biotype of genes with biomart:

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
biotype <- biomaRt::getBM(attributes=c("external_gene_name", "gene_biotype"), filters=c("external_gene_name"), values=union_genes, mart=mart)
biotype <- biotype %>% dplyr::group_by(external_gene_name) %>% dplyr::summarise(Type=paste(unique(gene_biotype), collapse=" | "))
missing2 <- union_genes[!union_genes %in% biotype$external_gene_name]

length(missing2)
# [1] 63413

colnames(biotype) <- c("Gene", "Type")
biotype_df <- rbind(biotype_df, biotype)

## How many genes are still unannotated?

missing <- intersect(missing1, missing2)
length(missing)
# [1] 42622

missing <- missing[missing!=""]

## Cross. ref remaining unannotated genes with txdb / GTF:

## How many unannoated genes are in txdb file?

hub <- AnnotationHub()
query(hub, c("homo sapiens","EnsDb"))
txdb <- hub[["AH98047"]]
txdb_df <- genes(txdb, return.type="data.frame")
txdb_df$gene_name <- toupper(txdb_df$gene_name)

sum(is.element(missing, txdb_df$gene_name))
# [1] 35

## How many unannotated genes are in GTF file?

txdb_df <- txdb_df[is.element(txdb_df$gene_name, missing),]
txdb_df <- txdb_df %>% dplyr::select(gene_name, gene_biotype)
colnames(txdb_df) <- c("Gene", "Type")

gtf_df <- makeTx2geneFromGtf("/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf")
gtf_df$gene_name <- toupper(gtf_df$gene_name)

sum(is.element(missing, gtf_df$gene_name))
# [1] 259

gtf_df <- gtf_df[is.element(gtf_df$gene_name, missing),]
gtf_df <- gtf_df %>% dplyr::select(gene_name, gene_type)
colnames(gtf_df) <- c("Gene", "Type")

## Add these annotations:

biotype_df <- rbind(biotype_df, gtf_df, txdb_df)
biotype_df$Gene <- toupper(biotype_df$Gene)
biotype_df <- biotype_df %>% dplyr::group_by(Gene) %>% dplyr::summarise(Type=paste(unique(Type), collapse=" | "))

missing <- missing[!is.element(toupper(missing), biotype_df$Gene)]
missing <- missing[!grepl("^LOC", missing)]
missing <- missing[!grepl("^RPL", missing)]
missing <- missing[!grepl("^RPS", missing)]
missing <- missing[!grepl("^SNO", missing)]
missing <- missing[!grepl("^SNAR", missing)]
missing <- missing[!grepl("P[0-9]*$", missing)]

## How genes are many still unannotated?

length(missing)
# [1] 24577

## How many genes annotated as PC?

biotype_df <- biotype_df[grep("protein", biotype_df$Type, ignore.case=T),]
colnames(biotype_df)[1] <- "SYMBOL"
dim(biotype_df)
# [1] 20051     2

fwrite(biotype_df[,-c(2)], file=paste0("protein_coding_", nrow(biotype_df), "_union_genes.csv"))
