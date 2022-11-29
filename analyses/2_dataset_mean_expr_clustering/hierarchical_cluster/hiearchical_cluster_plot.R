library(ggdendro)
library(dendextend)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(propr)

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")

plot_dendro <- function(dataset_expr, datinfo, data_type, expr_type, combat=F, sim_type, prop_metric=c("rho", "phs"), clust_method, top_n=NULL){
  
  datinfo <- datinfo[match(colnames(dataset_expr)[-c(1)], datinfo$Dataset),]
  
  if(!identical(datinfo$Dataset, colnames(dataset_expr)[-c(1)])){
    stop("!identical(datinfo$Dataset, colnames(dataset_expr)[-c(1)])")
  }
  
  n_genes <- nrow(dataset_expr)
  
  if(!is.null(top_n)){
    
    ## Select genes with highest mean expression across all datasets:
    
    mean_dataset_expr <- rowMeans(dataset_expr[,-c(1)])
    dataset_expr <- dataset_expr[order(-mean_dataset_expr),] %>% dplyr::slice(1:top_n)
    n_genes <- paste("top", nrow(dataset_expr), sep="_")
    
  } 
  
  if(sim_type=="prop"){
    
    sim <- propr(dataset_expr[,-c(1)], metric=prop_metric, symmetrize=T, p=1)@matrix
    sim_type <- paste(sim_type, prop_metric, sep="_")
    if(grepl("ph", prop_metric)){ # phs and phi are already distance measures
      sim <- 1-sim
    }
    
  } else {
    sim <- cor(dataset_expr[,-c(1)], method=sim_type, use="pairwise.complete.obs")
  }
  
  dendro <- as.dendrogram(stats::hclust(as.dist(1-sim), method=clust_method))
  
  dendro_info <- ggdendro::dendro_data(dendro)
  
  ## Make sure metadata matches order of dendrogram x-axis:
  
  datinfo <- datinfo[match(dendro_info$labels$label, datinfo$Dataset),]
  
  ## Assign plot-friendly labels to dendrogram leafs:
  
  labels(dendro) <- datinfo$Plot_Label
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/dataset_mean_expr_dendrogram_", data_type, "_", expr_type, "_", toupper(sim_type), "_", toupper(clust_method), "_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies_", n_genes, "_genes.pdf")
  
  if(combat){
    file_path <- gsub(paste0(expr_type, "_"), paste0(expr_type, "_COMBAT_"), file_path)
  }
  
  pdf(file=file_path, width=10, height=8)
  
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
  
  plot_title <- paste("Dataset Mean Expression Clustering")
  
  if(combat){
    plot_title <- paste0(plot_title, "\nBatch Corrected")
  }
  
  if(sim_type=="pearson"){
    
    y_lab <- upper_first(sim_type)
    
  } else if(sim_type=="prop"){
    
    y_lab <- "Proportionality"
    
  } else {
    
    y_lab <- "Cosine"
    
  }
  
  n_col <- 1; h <- T; sub_line <- 12; y <- -1.5
  if(n_elements>5){
    h <- F; n_col <- ceiling(n_elements/2)
  }
  
  if(n_elements>20){
    n_col <- 5; sub_line <- 11; y <- -1.7
  }
  
  marg <- c(23, 8, 4, 5)
  par(mar=marg, xpd=T)
  plot(dendro, main=plot_title, ylab=paste0("1-", y_lab), cex.main=1.3)
  mtext(paste(upper_first(clust_method), "linkage"), side=3, line=0, cex=1.1)
  mtext(plot_sub, side=1, cex=1.2, line=sub_line)
  legend("bottom", title=color_lab, legend=unique(datinfo[,color_var]), col=color_vals, ncol=n_col, pch=c(20, 20, 20, 20), bty="n", pt.cex=3, cex=1.2, text.col="black", horiz=h, inset=c(0, y))
  
}


