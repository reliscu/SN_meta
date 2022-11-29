library(ggplot2)
library(pheatmap)
library(plyr)
library(dplyr)
library(tidyverse) ## purrr::reduce()
library(scales) ## comma()
library(APAstyler) ## snip()
library(colorRamps)

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/code/misc/overlap_matrix_fxn.R")

feature_overlap_heatmap <- function(datinfo, data_type, expr_type, umi_cut_list=c(1, 5, 10, 20), jaccard=F, pc_genes=NULL){
  
  platform_color <- brewer_fxn(9)
  names(platform_color) <- unique(datinfo$Platform)
  facs_color <- c("#3492EA", "#D92F4B")
  names(facs_color) <- c("Y", "N")
  sampling_color <- c("#F4D823", "#8B40DC")
  names(sampling_color) <- c("Y", "N")
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/feature_overlap_heatmap_", data_type, "_", expr_type, "_datasets.pdf")
  
  if(jaccard){
    file_path <- gsub("overlap", "overlap_percent",  file_path)
  }
  
  pdf(file=file_path, width=15, height=15)

  for(i in 1:length(umi_cut_list)){
    
    umi_cut <- umi_cut_list[i]
    
    gene_list <- lapply(1:length(datinfo$Author_Counts_QC_PC_Genes), function(i){
      expr <- fread(datinfo$Author_Counts_QC_PC_Genes[i], data.table=F)
      mean_expr <- rowMeans(expr[,-c(1)])
      return(expr[mean_expr>=umi_cut,1])
    })
    names(gene_list) <- datinfo$Dataset
    
    overlap_mat <- overlap_matrix(gene_list, jaccard)
  
    datinfo$No.Genes <- lengths(gene_list)
    
    meta_df <- data.frame(
      Platform=datinfo$Platform,
      FACS=datinfo$FACS_Sorted,
      `Unbiased Sampling`=datinfo$Unbiased_Sampling,
      Temp=datinfo$No.Genes,
      check.names=F
    )
    
    colnames(meta_df)[ncol(meta_df)] <- paste("# Genes UMIs >=", format(umi_cut, scientific=F))
    
    if(umi_cut==0){
      colnames(meta_df)[ncol(meta_df)] <- gsub("=", "", colnames(meta_df)[ncol(meta_df)], fixed=T)
    }
    
    rownames(meta_df) <- datinfo$Dataset
    
    anno_colors <- list(
      Platform=platform_color,
      FACS=facs_color,
      `Unbiased Sampling`=sampling_color
    )
    
    union_genes <- length(Reduce(union, gene_list))
    intersection_genes <- length(Reduce(intersect, gene_list))
    
    plot_sub <- paste(
      format(union_genes, big.mark=","), "union genes,", format(intersection_genes, big.mark=","), "genes shared among all datatsets"
    )
    
    if(jaccard){
      
      overlap_mat[is.nan(overlap_mat)] <- 0
  
      if(i==1){
        
        p <- pheatmap(mat=overlap_mat, scale="none", cluster_rows=T, cluster_cols=T, silent=T)
        row_idx <- p$tree_row$order; col_idx <- p$tree_col$order
        row_labs <- p$tree_row$labels; col_labs <- p$tree_col$labels
        
      }
      
      overlap_mat <- overlap_mat[row_idx, col_idx]
      display_mat <- apply(overlap_mat, 2, function(x) sprintf("%.0f%%", x))
      
      plot_title <- paste(
        "Feature Overlap Percent\nRestricting to Protein Coding Genes with Mean UMIs >=", format(umi_cut, scientific=F), "\n", plot_sub
      )
      
      if(umi_cut==0){
        plot_title <- gsub("=", "", plot_title, fixed=T)
      }
      
      break_list <- seq(0, 100, 10)
      leg_color <- colorRampPalette(rev(brewer.pal(9, name="RdYlBu")))(length(break_list)-1)
      
      print(
        pheatmap(
          mat=overlap_mat, 
          main=plot_title, 
          scale="none", 
          cluster_rows=F, 
          cluster_cols=F, 
          labels_row=datinfo$Plot_Label[row_idx], 
          labels_col=datinfo$Plot_Label[col_idx], 
          border_color="black", 
          angle_col=45, 
          annotation_col=meta_df, 
          annotation_colors=anno_colors, 
          display_numbers=display_mat, 
          number_color="black", 
          color=leg_color, 
          legend_breaks=break_list, 
          breaks=break_list, 
          legend_labels=break_list, 
          cellwidth=21, 
          cellheight=21
        )
      )
      
      ## treeheight_row=0, 
      
    } else { ## if(jaccard){
      
      if(i==1){
        p <- pheatmap(mat=overlap_mat, scale="none", cluster_rows=T, cluster_cols=T, silent=T)
        row_idx <- p$tree_row$order; col_idx <- p$tree_col$order
        row_labs <- p$tree_row$labels; col_labs <- p$tree_col$labels
      }
      
      overlap_mat <- overlap_mat[row_idx, col_idx]
      
      plot_title <- paste(
        "Number of Shared Genes\nRestricting to Protein Coding Genes with Mean UMIs >=", format(umi_cut, scientific=F), "\n",  plot_sub
      )
      if(umi_cut==0){
        plot_title <- gsub("=", "", plot_title, fixed=T)
      }
      
      meta_df <- meta_df[,-ncol(meta_df)]
      
      if(floor(union_genes/1e3)==0){
        
        max_val <- 1e3
        break_size <- 100
        
      } else {
        
        max_val <- 18e3
        break_size <- 2e3
        
      }
      
      break_list <- seq(0, max_val, by=break_size)
      leg_color <- colorRampPalette(rev(brewer.pal(9, name="RdYlBu")))(length(break_list)-1)

      print(
        pheatmap(
          mat=overlap_mat, 
          main=plot_title, 
          scale="none", 
          cluster_rows=F, 
          cluster_cols=F, 
          labels_row=datinfo$Plot_Label[row_idx], 
          labels_col=datinfo$Plot_Label[col_idx], 
          border_color="black", 
          angle_col=45, 
          annotation_col=meta_df, 
          annotation_colors=anno_colors, 
          display_numbers=F, 
          color=leg_color, 
          breaks=break_list, 
          legend_breaks=break_list, 
          legend_labels=comma(break_list), 
          cellwidth=20, 
          cellheight=20
        )
      )
      
    } ## if(jaccard){} else {
    
  } ## for(i in 1:length(umi_cut_list)){

  dev.off()
  
} ## feature_overlap_heatmap <- function(

feature_overlap_boxplot <- function(datinfo, data_type, expr_type, umi_cut_list=c(1, 5, 10, 20), jaccard=F, pc_genes=NULL){
  
  df <- c()
  
  for(i in 1:length(umi_cut_list)){
    
    umi_cut <- umi_cut_list[i]
    
    dataset_genes <- genes_UMI_threshold(
      datinfo, expr_type, umi_cut, pc_genes
    )
    
    overlap_mat <- overlap_matrix(list=dataset_genes, jaccard)
    diag(overlap_mat) <- NA
    
    df_umi_cut <- reshape2::melt(overlap_mat)
    colnames(df_umi_cut)[ncol(df_umi_cut)] <- "Overlap"
    df_umi_cut$UMI_Cut <- umi_cut
    
    df_umi_cut <- df_umi_cut %>% dplyr::select(UMI_Cut, Overlap) %>% tidyr::drop_na()
    
    if(is.null(df)){
      df <- df_umi_cut
    } else {
      df <- rbind(df_umi_cut, df)
    }
    
  } ## for(i in 1:length(umi_cut_list)){
  
  df$UMI_Cut[df$UMI_Cut==5e-04] <- format(df$UMI_Cut[df$UMI_Cut==5e-04], scientific=F)
  df$UMI_Cut <- factor(df$UMI_Cut, levels=sort(unique(df$UMI_Cut)))
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/feature_overlap_boxplot_", data_type, "_", expr_type, "_datasets.pdf")
  
  if(jaccard){
    file_path <- gsub("overlap", "overlap_percent",  file_path)
  }
  
  pdf(file=file_path, width=8, height=6.5)
  
  if(jaccard){
    
    plot_title <- paste("Feature Overlap")
    
    y_max <- 100
    y_breaks <- seq(0, y_max, by=20)
    y_lab <- paste("% Overlap")
    y_labs <- comma(y_breaks)
    
  } else {
    
    plot_title <- paste("Number of Shared Features")
    
    y_max <- round_any(max(df$Overlap), 1e3)
    y_breaks <- seq(0, y_max, by=2e3)
    y_lab <- ("# Shared Features")
    y_labs <- comma(y_breaks)
    
  }

  plot_sub <- paste("Pairwise comparison between", n_distinct(datinfo$Dataset), "datasets\nRestricting to protein coding genes")
  
  print(
    ggplot(df, aes(x=UMI_Cut, y=Overlap, group=UMI_Cut)) +
      geom_violin(trim=F, color="#505050", fill="lightblue", width=1.1, show.legend=F) +
      geom_boxplot(notch=T, alpha=1, width=.03, outlier.shape=NA, color="black", fill="white", show.legend=F) +
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=.5, face="bold", size=15),
        plot.subtitle=element_text(hjust=.5, size=12, lineheight=1.1, margin=margin(t=4, b=8)),
        legend.title=element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.box="vertical",
        axis.title.x=element_text(size=14, color="black", face="bold", margin=margin(t=15)), 
        axis.text.x=element_text(size=12, margin=margin(t=5), angle=0),
        axis.ticks.x=element_line(size=.4),
        axis.title.y=element_text(size=12, color="black", margin=margin(r=15)),
        axis.ticks.y=element_line(size=.4),
        plot.margin = unit(c(2, 2, 2, 2), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab(expression("Restricting to Features with Mean UMIs ">="")) + ylab(y_lab) +
      scale_y_continuous(labels=y_labs, breaks=y_breaks, limits=c(0, y_max))
  )
  
  print(
    ggplot(df, aes(x=UMI_Cut, y=Overlap, group=UMI_Cut)) +
      geom_boxplot(notch=T, width=.45, fill="lightblue", outlier.size=.1, show.legend=F) +
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=.5, face="bold", size=15),
        plot.subtitle=element_text(hjust=.5, size=12, lineheight=1.1, margin=margin(t=4, b=8)),
        legend.title=element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.box="vertical",
        axis.title.x=element_text(size=14, color="black", face="bold", margin=margin(t=15)), 
        axis.text.x=element_text(size=12, margin=margin(t=5), angle=0),
        axis.ticks.x=element_line(size=.4),
        axis.title.y=element_text(size=12, color="black", margin=margin(r=15)),
        axis.ticks.y=element_line(size=.4),
        plot.margin = unit(c(2, 2, 2, 2), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab(expression("Restricting to Features with Mean UMIs ">="")) + ylab(y_lab) +
      scale_y_continuous(labels=y_labs, breaks=y_breaks, limits=c(0, y_max)) 
      #scale_fill_manual(values=colorRampPalette(rev(brewer.pal(9, "Blues")))(n_distinct(df$UMI_Cut)))
  )
  
  dev.off()
  
}