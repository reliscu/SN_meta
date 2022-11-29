library(ggdendro)
library(dendextend)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(propr)

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_data_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

plot_dendro <- function(class_expr, datinfo, cell_class, neu_subtypes, data_type=c("author_data"), expr_type=c("counts", "normalized_counts"), combat=F, sim_type, prop_metric=c("rho", "phs"), clust_method, top_n=NULL){
  
  ct_data <- prep_CT_data(datinfo, ct_df=class_expr, expr_type, top_n, ct_expr=NULL)
  
  datinfo <- ct_data[[1]]
  class_expr <- ct_data[[2]]
  
  if(is.null(top_n)){
    
    n_genes <- nrow(class_expr)
    
  } else {
    
    n_genes <- paste("top", nrow(class_expr), sep="_")
    
  }
  
  datinfo$Cell_Type <- paste(datinfo$Study, datinfo$Cell_Type)
  
  n_cts <- n_distinct(datinfo$Cell_Type)
  
  if(sim_type=="prop"){
    
    sim <- propr(class_expr[,-c(1)], metric=prop_metric, symmetrize=T, p=1)@matrix
    sim_type <- paste(sim_type, prop_metric, sep="_")
    if(grepl("ph", prop_metric)){ # phs and phi are already distance measures
      sim <- 1-sim
    }
    
  } else {
    sim <- cor(class_expr[,-c(1)], method=sim_type, use="pairwise.complete.obs")
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
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/", cell_class, "_subtypes_mean_expr_dendrogram_", data_type, "_", expr_type, "_", toupper(sim_type), "_", toupper(clust_method), "_", n_cts, "_CTs_", n_distinct(datinfo$Dataset), "_datasets_", n_genes, "_intersection_genes.pdf")
  
  if(combat){
    file_path <- gsub(paste0(expr_type, "_"), paste0(expr_type, "_COMBAT_"), file_path)
  }
  
  pdf(file=file_path, width=14, height=10)
  
  color_var <- "Class_Level1"
  color_lab <- "Subtype"
  dendro_template(datinfo, dendro, color_var, color_lab, n_genes=nrow(class_expr), n_cts)
  
  color_var <- "Class_Level2"
  color_lab <- "Subtype"
  dendro_template(datinfo, dendro, color_var, color_lab, n_genes=nrow(class_expr), n_cts)
  
  color_var <- "Hodge_2018_Annotation"
  color_lab <- "Hodge 2018 Subtype"
  dendro_template(datinfo, dendro, color_var, color_lab, n_genes=nrow(class_expr), n_cts)
  
  color_var <- "Bakken_2019_Annotation"
  color_lab <- "Bakken 2019 Subtype"
  dendro_template(datinfo, dendro, color_var, color_lab, n_genes=nrow(class_expr), n_cts)
  
  dev.off()
  
}

dendro_template <- function(datinfo, dendro, color_var, color_lab, n_genes, n_cts){
  
  plot_title <- paste(cell_class_full_name(cell_class), "Mean Expression Clustering")
  if(combat){
    plot_title <- paste0(plot_title, "\nBatch Corrected")
  }
  
  plot_sub <- paste(n_cts, tolower(cell_class_full_name(cell_class)), "cell types found in", n_distinct(datinfo$Dataset), "datasets from", n_distinct(datinfo$Study), "studies\n", comma(n_genes), "protein coding genes")
  
  subtype_labs <- datinfo[,color_var]
  neu_cols <- brewer_fxn(n_distinct(subtype_labs))
  neu_cols <- neu_cols[match(subtype_labs, unique(subtype_labs))]
  neu_cols[is.na(subtype_labs)] <- "#000000"
  subtype_labs[is.na(subtype_labs)] <- "NA"
  dendextend::labels_colors(dendro) <- neu_cols

  if(sim_type=="pearson"){
    
    y_lab <- upper_first(sim_type)
    
  } else if(sim_type=="prop"){
    
    y_lab <- "Proportionality"
    
  } else {
    
    y_lab <- "Cosine"
    
  }
  
  n_elements <- n_distinct(datinfo[,color_var])
  
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
    mtext(paste(clust_method, "linkage"), side=3, line=-.5, cex=1.1)
    mtext(plot_sub, side=1, line=sub_line, cex=1.2)
    legend("bottom",  title=color_lab, legend=unique(subtype_labs), col=unique(neu_cols), ncol=4, pch=c(20, 20, 20, 20), bty="n", pt.cex=3, cex=1.2, text.col="black", horiz=F, inset=c(.8, y))
    
  } ## if(color_var=="Cell Type"){} else {
  
}


