library(ggdendro)
library(dendextend)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(propr)

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

plot_dendro <- function(dataset_class_expr, datinfo, cell_class, data_type, expr_type, sim_type, prop_metric=c("rho", "phs"), clust_method, top_n=NULL){

  datinfo <- datinfo[match(colnames(dataset_expr)[-c(1)], datinfo$Dataset),]
  
  if(!identical(datinfo$Dataset, colnames(dataset_expr)[-c(1)])){
    stop("!identical(datinfo$Dataset, colnames(dataset_expr)[-c(1)])")
  }

  n_genes <- nrow(dataset_class_expr)
  
  if(!is.null(top_n)){
    
    ## Select genes with highest CELL CLASS mean expression across all datasets:
    
    mean_dataset_class_expr <- rowMeans(dataset_class_expr[,-c(1)], na.rm=T)
    dataset_class_expr <- dataset_class_expr[order(-mean_dataset_class_expr),] %>% dplyr::slice(1:top_n)
    n_genes <- paste("top", nrow(dataset_class_expr), sep="_")
    
  }
  
  if(sim_type=="prop"){
    
    sim <- propr(dataset_class_expr[,-c(1)], metric=prop_metric, symmetrize=T, p=1)@matrix
    sim_type <- paste(sim_type, prop_metric, sep="_")
    if(grepl("ph", prop_metric)){ # phs and phi are already distance measures
      sim <- 1-sim
    }
    
  } else {
    sim <- cor(dataset_class_expr[,-c(1)], method=sim_type, use="pairwise.complete.obs")
  }
  
  dendro <- as.dendrogram(stats::hclust(as.dist(1-sim), method=clust_method))

  dendro_info <- ggdendro::dendro_data(dendro)
  datinfo <- datinfo[match(dendro_info$labels$label, datinfo$Dataset),] ## Make sure metadata matches order of dendrogram x-axis
  labels(dendro) <- datinfo$Plot_Label ## Assign plot-friendly labels to dendrogram leafs

  pdf(file=paste0("figures/", data_type, "/", expr_type, "/", cell_class, "_dataset_mean_expr_dendrogram_", data_type, "_", expr_type, "_", toupper(sim_type), "_", toupper(clust_method), "_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies_", n_genes, "_genes.pdf"), width=14, height=10)

  color_var <- "Platform"; color_lab <- "Platform"
  dendro_template(datinfo, dendro, color_var, color_lab, n_genes)
  
  color_var <- "Study"; color_lab <- "Study"
  dendro_template(datinfo, dendro, color_var, color_lab, n_genes)
  
  dev.off()
  
}

dendro_template <- function(datinfo, dendro, color_var, color_lab, n_genes){
  
  plot_sub <- paste(n_distinct(datinfo$Dataset), "datasets from", n_distinct(datinfo$Study), "studies\n", comma(n_genes), "protein coding genes")
  
  n_elements <- n_distinct(datinfo[,color_var])
  
  if(n_elements==2){
    color_vals <- c("#3492EA", "#D92F4B")
  } else if(n_elements<8){
    color_vals <- brewer.pal(n_elements, "Set2")
  } else {
    color_vals <- brewer_fxn(n_elements)
  }
  
  temp <- cbind(unique(datinfo[,color_var]), color_vals)
  idx <- match(datinfo[,color_var], temp[,1])
  dendextend::labels_colors(dendro) <- color_vals[idx]
  
  plot_title <- paste("Dataset", cell_class_full_name(cell_class), "Mean Expression Clustering")
  
  # if(combat){
  #   plot_title <- paste0(plot_title, "\nBatch Corrected")
  # }
  
  if(sim_type=="pearson"){
    
    y_lab <- upper_first(sim_type)
    
  } else if(sim_type=="prop"){
    
    y_lab <- "Proportionality"
    
  } else {
    
    y_lab <- "Cosine"
    
  }
  
  n_col <- 1; h <- T; y <- -.75
  if(color_lab=="Platform"){
    n_col <- 4; h <- F; y <- -.76
  } else if(n_distinct(dendextend::labels_colors(dendro))>=10){
    n_col <- 4; h <- F; y <- -.86
  }
  
  marg <- c(23, 5, 4, 2)
  sub_line <- 9.5
  
  par(mar=marg, xpd=T)
  plot(dendro, main=plot_title, ylab=paste0("1-", y_lab), cex.main=1.3)
  mtext(paste(clust_method, "linkage"), side=3, line=-.5, cex=1.1)
  mtext(plot_sub, side=1, line=sub_line, cex=1.2)
  legend("bottom", title=color_lab, unique(datinfo[,color_var]), col=color_vals, ncol=n_col, pch=c(20, 20, 20, 20), bty="n", pt.cex=3, cex=1.2, text.col="black", horiz=h, inset=c(.8, y))
  
  
}
