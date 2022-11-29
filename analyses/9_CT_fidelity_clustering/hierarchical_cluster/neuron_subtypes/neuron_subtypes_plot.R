library(ggdendro)
library(dendextend)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(propr)
library(igraph) ## compare()
library(scales) ## comma()

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_data_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
# source("/home/rebecca/SCSN_meta_analysis/code/clust_metric_fxns.R")

plot_neuron_subtypes <- function(
  ct_fid, 
  datinfo, 
  data_type=c("author_data"), 
  expr_type=c("counts", "normalized_counts"), 
  which_genes=c("intersection", "union"), 
  gene_list=NULL,
  sim_type, 
  prop_metric=c("rho", "phs"), 
  clust_method, 
  top_n=NULL,
  ct_expr,
  n_clusters,
  neu_subtypes,
  prop_scaled=F
){
  
  ct_data <- prep_CT_data(
    datinfo, 
    ct_df=ct_fid, 
    expr_type, 
    which_genes,
    pc_genes=T,
    top_n, 
    ct_expr=NULL
  )
  
  datinfo <- ct_data[[1]]
  ct_fid <- ct_data[[2]]
  
  if(!is.null(gene_list)){
    ct_fid <- ct_fid[is.element(ct_fid$SYMBOL, gene_list),]
  }
  
  neu_subtypes <- neu_subtypes %>%
    dplyr::select(
      Label, 
      Class_Level1, Class_Level2,
      Hodge_2018_Annotation, 
      Bakken_2019_Annotation
    )
  
  neu_subtypes$Label <- make.names(neu_subtypes$Label)
  datinfo <- merge(datinfo, neu_subtypes, by="Label", all.x=T)
  
  if(is.null(top_n)){
    
    n_genes <- nrow(ct_fid)
    
  } else {
    
    n_genes <- paste("top", nrow(ct_fid), sep="_")
    
  }
  
  n_cts <- n_distinct(paste(datinfo$Study, datinfo$Cell_Type))
  
  if(sim_type=="prop"){
    
    sim <- propr(ct_fid[,-c(1)], metric=prop_metric, symmetrize=T, p=1)@matrix
    sim_type <- paste(sim_type, prop_metric, sep="_")
    if(grepl("ph", prop_metric)){ # phs and phi are already distance measures
      sim <- 1-sim
    }
    
  } else {
    sim <- cor(ct_fid[,-c(1)], method=sim_type, use="pairwise.complete.obs")
  }
  
  clust <- stats::hclust(
    as.dist(1-sim), method=clust_method
  )
  
  dendro <- as.dendrogram(clust)
  
  dend_labs <-  ggdendro::dendro_data(dendro)$labels$label
  datinfo <- datinfo[match(dend_labs, datinfo$Label),] ## Make sure metadata matches order of dendrogram x-axis
  labels(dendro) <- datinfo$Plot_Label_CT ## Assign plot-friendly labels to dendrogram leafs
  dendro <- dendextend::set(dendro, "labels_cex", .4) ## Dendrogram leaf font size
  
  true_labs <- datinfo$Cell_Class
  clust_labs <- stats::cutree(clust, k=n_clusters)
  clust_labs <- clust_labs[match(datinfo$Label, make.names(names(clust_labs)))] ## Make sure cell type cluster assignment matches datinfo
  
  nmi_idx <- signif(igraph::compare(true_labs, clust_labs, method="nmi"), 2)
  rand_idx <- signif(igraph::compare(true_labs, clust_labs, method="rand"), 2)
  clust_metrics <- paste0("NMI=", nmi_idx, ", RAND=", rand_idx)
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/neuronal_subtypes_fidelity_dendrogram_", n_clusters, "_clusters_", data_type, "_", expr_type, "_", toupper(sim_type), "_", toupper(clust_method), "_", n_cts, "_CTs_", n_distinct(datinfo$Dataset), "_datasets_", n_genes, "_intersection_genes.pdf")
  
  if(prop_scaled){
    file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
  }
  
  pdf(file=file_path, width=14, height=9)
  
  color_var <- "Class_Level1"
  color_lab <- "Subtype"
  dendro_template(datinfo, dendro, clust, clust_metrics, color_var, color_lab, n_genes=nrow(ct_fid), n_cts)
  
  color_var <- "Class_Level2"
  color_lab <- "Subtype"
  dendro_template(datinfo, dendro, clust, clust_metrics, color_var, color_lab, n_genes=nrow(ct_fid), n_cts)
  
  color_var <- "Hodge_2018_Annotation"
  color_lab <- "Hodge 2018 Subtype"
  dendro_template(datinfo, dendro, clust, clust_metrics, color_var, color_lab, n_genes=nrow(ct_fid), n_cts)
  
  color_var <- "Bakken_2019_Annotation"
  color_lab <- "Bakken 2019 Subtype"
  dendro_template(datinfo, dendro, clust, clust_metrics, color_var, color_lab, n_genes=nrow(ct_fid), n_cts)
  
  dev.off()
  
}

dendro_template <- function(datinfo, dendro, clust, clust_metrics, color_var, color_lab, n_genes, n_cts){
  
  plot_title <- "Cell Type Fidelity Clustering"
  if(prop_scaled){
    plot_title <- gsub("Fidelity", "Proportional Fidelity", plot_title)
  }

  plot_sub <- paste(comma(n_cts), "cell types found in", n_distinct(datinfo$Dataset), "datasets from", n_distinct(datinfo$Study), "studies\n", comma(n_genes), "protein coding genes")
  
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
  
  n_col <- 1; h <- T; y <- -.9
  if(n_distinct(subtype_labs)>5){
    n_col <- 4; h <- F; y <- -1.15
  }
  
  par(mar=c(23, 5, 4, 2), xpd=T)
  plot(dendro, main=plot_title, ylab=paste0("1-", y_lab))
  rect.hclust(tree=clust, k=n_clusters)
  mtext(paste(n_clusters, "clusters,", paste0(clust_metrics, ","), clust_method, "linkage"), side=3, line=-.25)
  mtext(plot_sub, side=1, line=10, cex=1.2)
  legend("bottom", title=color_lab, legend=unique(subtype_labs), col=unique(neu_cols), ncol=n_col, pch=c(20, 20, 20, 20), bty="n", pt.cex=3, cex=1.2, text.col="black", horiz=h, inset=c(.8, y))
  
}

