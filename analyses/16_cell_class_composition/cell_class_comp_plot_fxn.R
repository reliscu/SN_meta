library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(cowplot) # plot_grid
library(gridExtra) ## grid.arrange

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")

cell_comp_dendro <- function(datinfo, df, sim_type, clust_method=c("complete", "ward.D2", "average")){
 
  mat <- reshape2::dcast(df, Cell_Class~Dataset, value.var="Proportion")
  mat[is.na(mat)] <- 0
  
  if(sim_type=="prop"){
    
    sim <- propr(mat[,-c(1)], metric="rho", symmetrize=T, p=1)@matrix
    
  } else {
    
    sim <- cor(mat[,-c(1)], method=sim_type, use="pairwise.complete.obs")
    
  }
  
  dendro <- as.dendrogram(
    stats::hclust(as.dist(1-sim), method=clust_method)
  )
  
  dendro_info <- ggdendro::dendro_data(dendro)
  
  ## Make sure metadata matches order of dendrogram x-axis:
  
  datinfo <- datinfo[match(dendro_info$labels$label, datinfo$Dataset),] 
  
  ## Assign plot-friendly labels to dendrogram leafs:
  
  labels(dendro) <- datinfo$Plot_Label 
  
  pdf(file=paste0("figures/", data_type, "/dataset_cell_class_proportion_dendrogram_", toupper(sim_type), "_", toupper(clust_method), "_", data_type, "_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies.pdf"), width=10, height=8)
  
  color_var <- "Platform"; color_lab <- "Platform"
  dendro_template(datinfo, dendro, color_var, color_lab)
  
  color_var <- "Unbiased_Sampling"; color_lab <- "Unbiased Sampling"
  dendro_template(datinfo, dendro, color_var, color_lab)
  
  dev.off()
  
}

dendro_template <- function(datinfo, dendro, color_var, color_lab){
  
  plot_sub <- paste(n_distinct(datinfo$Dataset), "datasets from", n_distinct(datinfo$Study), "studies")
  
  n_elements <- n_distinct(datinfo[,color_var])
  
  if(n_elements==2){
    color_vals <- c("#3492EA", "#D92F4B")
  } else if(n_elements<8){
    color_vals <- brewer.pal(n_elements, "Set1")
  } else {
    color_vals <- brewer.fxn(n_elements)
  }
  
  temp <- cbind(unique(datinfo[,color_var]), color_vals)
  idx <- match(datinfo[,color_var], temp[,1])
  dendextend::labels_colors(dendro) <- color_vals[idx]
  
  n_col <- 1; h <- T; sub_line <- 10
  if(n_elements>5){
    h <- F; n_col <- ceiling(n_elements/2); sub_line <- 9.4
  }
  
  plot_title <- paste("Distribution of Cell Class Assignments")
  
  if(sim_type=="pearson"){
    
    y_lab <- upper_first(sim_type)
    
  } else if(sim_type=="prop"){
    
    y_lab <- "Proportionality"
    
  } else {
    
    y_lab <- "Cosine"
    
  }
  
  marg <- c(20, 8, 4, 5)
  par(mar=marg, xpd=T)
  plot(dendro, main=plot_title, ylab=paste0("1-", y_lab), cex.main=1.3)
  mtext(paste(upper_first(clust_method), "linkage"), side=3, line=0, cex=1.1)
  mtext(plot_sub, side=1, cex=1.2, line=sub_line)
  legend("bottom", title=color_lab, legend=unique(datinfo[,color_var]), col=color_vals, ncol=n_col, pch=c(20, 20, 20, 20), bty="n", pt.cex=3, cex=1.2, text.col="black", horiz=h, inset=c(0, -1.035))
  
}

cell_comp_heatmap <- function(datinfo, df, sim_type){

  mat <- reshape2::dcast(df, Cell_Class~Dataset, value.var="Proportion")
  mat[is.na(mat)] <- 0

  if(sim_type=="prop"){
    
    sim <- propr(mat[,-c(1)], metric="rho", symmetrize=T, p=1)@matrix
    
  } else {
    
    sim <- cor(mat[,-c(1)], method=sim_type, use="pairwise.complete.obs")
    
  }
  
  datinfo <- datinfo[match(colnames(sim), datinfo$Dataset),]
  
  platform_color <-  brewer.pal(n_distinct(datinfo$Platform), "Set1")
  names(platform_color) <- unique(datinfo$Platform)
  facs_color <- c("#3492EA", "#D92F4B")
  names(facs_color) <- c("Y", "N")
  sampling_color <- c("#F4D823", "#8B40DC")
  names(sampling_color) <- c("Y", "N")

  class_dfs <- lapply(unique(df$Cell_Class), function(class){
    df[is.element(df$Cell_Class, class),]
  })
  
  class_df <- join_all(class_dfs, by="Dataset")
  class_df <- class_df[match(datinfo$Dataset, class_df$Dataset),]
  
  idx_fxn <- function(cell_class){
    unname(which(apply(class_df, 2, function(x) sum(grepl(cell_class, x))>0)))+1
  }
  
  meta_df <- data.frame(
    Platform=datinfo$Platform,
    FACS=datinfo$FACS_Sorted,
    `Unbiased Sampling`=datinfo$Unbiased_Sampling,
    ASC=class_df[,idx_fxn("ASC")],
    END=class_df[,idx_fxn("END")],
    EXC=class_df[,idx_fxn("EXC")],
    INH=class_df[,idx_fxn("INH")],
    MIC=class_df[,idx_fxn("MIC")],
    NEU=class_df[,idx_fxn("NEU")],
    OG=class_df[,idx_fxn("OG")],
    OPC=class_df[,idx_fxn("OPC")],
    PER=class_df[,idx_fxn("PER")],
    VSMC=class_df[,idx_fxn("VSMC")],
    check.names=F
  )
  meta_df[is.na(meta_df)] <- 0
  rownames(meta_df) <- datinfo$Dataset
  
  anno_colors <- list(
    Platform=platform_color,
    FACS=facs_color,
    `Unbiased Sampling`=sampling_color
  )
  
  break_list <- seq(-1, 1, by=.2)
  leg_color <- colorRampPalette(
    rev(brewer.pal(9, name="RdYlBu"))
  )(length(break_list)-1)
  
  
  if(sim_type=="pearson"){
    
    plot_title <- paste("Distribution of Cell Class Assignments\nPearson Similarity")
    leg_lab <- paste0(upper_first(sim_type), "\nCorrelation")
    
  } else if(sim_type=="prop"){
    
    plot_title <- paste("Distribution of Cell Class Assignments\nProp Rho Similarity")
    leg_lab <- "Proportionality"
    
  } else {
    
    plot_title <- paste("Distribution of Cell Class Assignments\nCosine Similarity")
    leg_lab <- "Cosine Similarity"
    
  }
  
  pdf(file=paste0("figures/", data_type, "/dataset_cell_class_proportion_heatmap_", toupper(sim_type), "_", data_type, "_", expr_type, "_datasets.pdf"), width=12, height=11)
  
  print(
    ComplexHeatmap::pheatmap(mat=sim, name=leg_lab, main=plot_title, scale="none", cluster_rows=T, cluster_cols=T, labels_row=datinfo$Plot_Label, labels_col=datinfo$Plot_Label, border_color="black", angle_col="45", annotation_col=meta_df, annotation_colors=anno_colors, display_numbers=F, color=leg_color, cellwidth=15, cellheight=15) # breaks=break_list, legend_breaks=break_list
  )
  
  dev.off()
  
  
}

cell_comp_barplot <- function(datinfo, df){
  
  df <- merge(df, datinfo, by="Dataset")

  ## Arrange datasets in order of EXC proportion:
  
  ds_order <- df %>% 
    dplyr::filter(Cell_Class=="EXC") %>%
    dplyr::arrange(Proportion)
  
  df$Plot_Label <- factor(df$Plot_Label, levels=c("Grubman TC EC", ds_order$Plot_Label))
  df$Cell_Class <- factor( df$Cell_Class, levels=c(sort(cell_classes), "Other"))
  df$Unbiased_Sampling[df$Unbiased_Sampling=="Y"] <- "Unbiased Sampling"
  df$Unbiased_Sampling[df$Unbiased_Sampling=="N"] <- "Biased Sampling"
 
  pdf(file=paste0("figures/", data_type, "/dataset_cell_class_proportion_barplot_", n_distinct(df$Dataset), "_datasets_", n_distinct(df$Study), "_studies.pdf"), width=9, height=7)
  
  plot_title <- paste("Cell Class Proportion")
  plot_sub <- paste(n_distinct(df$Dataset), "datasets from", n_distinct(df$Study), "studies")

  color_vec <- c(colorRampPalette(brewer.pal(11, "Paired"), interpolate="spline")(n_distinct(df$Cell_Class)-1), "#989898")
  #color_vec <- c("#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A", "#E31A1C", "#FB9A99", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#989898")
  
  p <- ggplot(df, aes(x=Plot_Label, y=Proportion, fill=Cell_Class)) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(
      plot.title=element_text(hjust=.5, face="bold", size=15),
      plot.subtitle=element_text(color="#616161", hjust=.5, size=12, lineheight=1, margin=margin(t=4, b=15)),
      legend.position="right",
      legend.box="vertical",
      legend.title=element_blank(),
      legend.text=element_text(size=10),
      axis.title.x=element_blank(), 
      axis.text.x=element_text(color="black", angle=90, size=9, hjust=1, vjust=.5, margin=margin(t=-5, r=0)),
      axis.line.x.bottom=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_text(size=13, color="black", margin=margin(r=8)),
      axis.ticks.y=element_blank(),
      axis.line.y.left=element_blank(),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    guides(fill=guide_legend(title.position="top")) +
    scale_fill_manual(values=color_vec)
  
  print(p)
  
  print(
    p + facet_grid(.~Unbiased_Sampling, scales="free_x") +
      theme(
        strip.text.x=element_text(size=12, face="bold"), 
        strip.background=element_rect(color="white")
      )
  )
  
  dev.off()
  
}

  
# p2 <- ggplot(df, aes(x=Plot_Label, y=1, fill=Platform)) +
#   geom_bar(stat="summary") +
#   theme_void() 
# 
# legend <- plot_grid(get_legend(p2), get_legend(p), ncol=1)
# 
# p <- p + theme(legend.position="none", axis.ticks=element_blank(), axis.text.y=element_blank())
# p2 <- p2 + theme(legend.position="none")
# 
# #plot <- plot_grid(p2, p, align="v", ncol=1, axis="tb", rel_heights=c(2, 15))
# #plot_grid(plot, legend, nrow=1, rel_widths=c(10, 1.5), )
# 
# p2/p
# 
# plot_grid(plot, legend, nrow=1, rel_widths=c(10, 1.5))

