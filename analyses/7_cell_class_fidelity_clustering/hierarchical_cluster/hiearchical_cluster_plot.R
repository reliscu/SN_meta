library(ggdendro)
library(dendextend)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(propr)

palette <- readRDS("/home/rebecca/SCSN_meta_analysis/palette_15.RDS")

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_data_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

plot_dendro <- function(class_fid, datinfo, cell_class, data_type=c("author_data"), expr_type=c("counts", "normalized_counts"), prop_scaled=F, which_genes=c("intersection", "union"), gene_list=NULL, sim_type, prop_metric=c("rho", "phs"), clust_method, top_n=NULL, class_expr=NULL){
  
  ct_data <- prep_CT_data(
    datinfo, 
    ct_df=class_fid, 
    expr_type, 
    which_genes, 
    pc_genes=T,
    top_n, 
    ct_expr=class_expr
  )
  
  datinfo <- ct_data[[1]]
  class_fid <- ct_data[[2]]
  
  if(!is.null(gene_list)){
    class_fid <- class_fid[is.element(class_fid$SYMBOL, gene_list),]
  }
  
  if(is.null(top_n)){
    
    n_genes <- nrow(class_fid)
    
  } else {
    
    n_genes <- paste("top", nrow(class_fid), sep="_")
    
  }
  
  datinfo$Cell_Type <- paste(datinfo$Study, datinfo$Cell_Type)
  
  n_cts <- n_distinct(datinfo$Cell_Type)
  
  if(sim_type=="prop"){
    
    sim <- propr(class_fid[,-c(1)], metric=prop_metric, symmetrize=T, p=1)@matrix
    sim_type <- paste(sim_type, prop_metric, sep="_")
    if(grepl("ph", prop_metric)){ # phs and phi are already distance measures
      sim <- 1-sim
    }
    
  } else {
    sim <- cor(class_fid[,-c(1)], method=sim_type, use="pairwise.complete.obs")
  }
  
  dendro <- as.dendrogram(stats::hclust(as.dist(1-sim), method=clust_method))
  
  dendro_info <- ggdendro::dendro_data(dendro)
  
  ## Make sure metadata matches order of dendrogram x-axis:
  
  datinfo$Label <- make.names(datinfo$Label)
  datinfo <- datinfo[match(dendro_info$labels$label, datinfo$Label),]
  
  ## Assign plot-friendly labels to dendrogram leafs:
  
  labels(dendro) <- datinfo$Plot_Label_CT
  
  leaf_cex <- .75
  if(nobs(dendro)>75){
    leaf_cex <- .5
  }
  
  dendro <- dendextend::set(dendro, "labels_cex", leaf_cex) ## Set dendrogram leaf font size
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/", cell_class, "_fidelity_dendrogram_", data_type, "_", expr_type, "_", toupper(sim_type), "_", toupper(clust_method), "_", n_cts, "_CTs_", n_distinct(datinfo$Dataset), "_datasets_", n_genes, "_intersection_genes.pdf")
  
  if(prop_scaled){
    file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
  }
  
  pdf(file=file_path, width=14, height=10)
  
  color_var <- "Platform"; color_lab <- "Platform"
  dendro_template(datinfo, dendro, color_var, color_lab, n_genes, n_cts)
  
  color_var <- "Study"; color_lab <- "Study"
  dendro_template(datinfo, dendro, color_var, color_lab, n_genes, n_cts)
  
  color_var <- "Cell_Type"; color_lab <- "Cell Type"
  dendro_template(datinfo, dendro, color_var, color_lab, n_genes, n_cts)
  
  dev.off()
  
}

dendro_template <- function(datinfo, dendro, color_var, color_lab, n_genes, n_cts){
  
  plot_title <- paste(cell_class_full_name(cell_class), "Fidelity Clustering")
  if(prop_scaled){
    plot_title <- gsub("Fidelity", "Proportional Fidelity", plot_title)
  }
  
  plot_sub <- paste(n_cts, tolower(cell_class_full_name(cell_class)), "cell types found in", n_distinct(datinfo$Dataset), "datasets from", n_distinct(datinfo$Study), "studies\n", comma(n_genes), "protein coding genes")
  
  n_elements <- n_distinct(datinfo[,color_var])
  
  if(n_elements==2){
    color_vals <- c("#3492EA", "#D92F4B")
  } else if(n_elements<8){
    color_vals <- brewer.pal(n_elements, "Set1")
  } else {
    color_vals <- brewer_fxn(n_elements)
  }
  
  temp <- cbind(unique(datinfo[,color_var]), color_vals)
  idx <- match(datinfo[,color_var], temp[,1])
  dendextend::labels_colors(dendro) <- color_vals[idx]
  
  if(sim_type=="pearson"){
    
    y_lab <- upper_first(sim_type)
    
  } else if(sim_type=="prop"){
    
    y_lab <- "Proportionality"
    
  } else {
    
    y_lab <- "Cosine"
    
  }
  
  if(nobs(dendro)<75){
    
    marg <- c(30, 8, 4, 5)
    
    n_col <- 2; sub_line <- 17; y <- -1.6
    if(n_elements>10){
      y <- -1.75
    }
    
  } else { ## if(nobs(dendro)<75){ 
    
    marg <- c(27, 8, 4, 5)
    
    sub_line <- 12.5; y <- -1.1
    if(n_elements>10){
      y <- -1.2
    }
    
  } ## if(nobs(dendro)<75){} else {
  
  if(color_var=="Cell_Type"){
    
    par(mar=c(21, 5, 4, 2), xpd=T)
    plot(dendro, main=plot_title, ylab=paste0("1-", y_lab), cex.main=1.3)
    mtext(paste(clust_method, "linkage"), side=3, line=-.5, cex=1.1)
    mtext(plot_sub, side=1, line=sub_line, cex=1.2)
    
  } else { ## if(color_var=="Cell Type"){
    
    par(mar=marg, xpd=T)
    plot(dendro, main=plot_title, ylab=paste0("1-", y_lab), cex.main=1.3)
    mtext(paste(clust_method, "linkage"), side=3, line=0, cex=1.1)
    mtext(plot_sub, side=1, line=sub_line, cex=1.2)
    legend("bottom", title=color_var, legend=unique(datinfo[,color_var]), col=color_vals, ncol=4, pch=c(20, 20, 20, 20), bty="n", pt.cex=3, cex=1.2, text.col="black", horiz=F, inset=c(.8, y))
    
  } ## if(color_var=="Cell Type"){} else {
  
}




