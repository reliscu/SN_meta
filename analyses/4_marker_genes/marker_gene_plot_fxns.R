library(ggplot2)
library(plyr) ## round_any()
library(dplyr)
library(scales) ## comma()

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/continuous_plot_breaks_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

library(pheatmap)
library(data.table)
library(plyr) ## round_any
library(dplyr)
library(scales)
#library(plotscale)

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/code/misc/overlap_matrix_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

marker_overlap_heatmap <- function(datinfo, cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"), legend, gene_sets, jaccard=F, pc_genes, min_size=5, remove_subclust=T){
  
  datinfo$Study <- paste(datinfo$First_Author, datinfo$Year)
  datinfo$Study[is.element(datinfo$PubMedID, "35165441")] <- paste(datinfo$Study[is.element(datinfo$PubMedID, "35165441")], " ")
  datinfo <- datinfo %>% dplyr::select(-Region)
  
  legend <- legend[legend$Cortical==T,]
  
  legend$Cell_Type <- gsub("^.*_[0-9]{4}_", "", legend$SetName)
  legend$Cell_Type <- gsub("_L[0-9]_VS_ALL", "", legend$Cell_Type)
  legend$Cell_Type <- gsub("_COMBO_VS_ALL", "", legend$Cell_Type)
  legend$Cell_Type <- gsub("_SINGLE_VS_ALL", "", legend$Cell_Type)
  legend$Cell_Type <- gsub("_SINGLE_VS.*", "", legend$Cell_Type)
  
  for(j in 1:length(cell_classes)){
    
    cell_class <- cell_classes[j]
    print(cell_class)
    class_name <- cell_class_full_name(cell_class)
    
    file_path <- paste0("figures/", cell_class, "_number_of_marker_genes_overlap_heatmap.pdf")
    
    if(jaccard){
      file_path <- gsub("number_of", "jaccard",  file_path)
    }
    
    class_legend <- legend[is.element(legend$Cell_Class, cell_class),]
    class_legend <- class_legend[class_legend$SetSize>=min_size,]
    
    if(remove_subclust){
      class_legend <- class_legend[!grepl("_DE_", class_legend$SetName),]
      class_legend <- class_legend[!grepl("COMBO", class_legend$SetName),]
      class_legend <- class_legend[!grepl("VS_L[0-9]", class_legend$SetName),]
      class_legend <- class_legend[!grepl("L[0-9]_VS", class_legend$SetName),]
    }
    
    if(nrow(class_legend)<2){
      next
    }
    
    class_legend <- merge(class_legend, datinfo, by.x="PubMed", by.y="PubMedID")
    class_legend <- class_legend %>% dplyr::arrange(SetName)
    class_legend <- class_legend[!duplicated(class_legend$SetName),]
    
    class_sets <- gene_sets[match(class_legend$SetID, names(gene_sets))]
    
    if(!identical(class_legend$SetID, names(class_sets))){
      stop("!identical(names(class_sets, class_legend$SetID))")
    }
    
    class_legend$Plot_Label <- gsub("(L[0-9]{1})(_)([0-9])", "\\1\\-\\3", class_legend$Cell_Type)
    class_legend$Plot_Label <- gsub("_", " ", class_legend$Plot_Label)
    class_legend$Plot_Label <- paste(
      class_legend$First_Author, class_legend$Plot_Label
    )
    
    class_legend$Platform[is.element(class_legend$First_Author, "Bakken")] <- "10x Chromium V3"
    class_legend$Platform[is.element(class_legend$First_Author, "Luo")] <- "snmCAT-seq"
    
    # class_legend <- class_legend %>%
    #   dplyr::group_by(Plot_Label) %>%
    #   dplyr::mutate(ID=1:n())
    # class_legend$Plot_Label[class_legend$ID>1] <- paste(
    #   class_legend$Plot_Label[class_legend$ID>1], class_legend$ID[class_legend$ID>1]
    # )
    
    study_color <-  colorRampPalette(brewer.pal(9, "Set1"), interpolate="spline")(n_distinct(class_legend$Study))
    names(study_color) <- unique(class_legend$Study)
    platform_color <-  colorRampPalette(brewer.pal(8, "Set2"))(n_distinct(class_legend$Platform))
    names(platform_color) <- unique(class_legend$Platform)
    facs_color <- c("#3492EA", "#D92F4B")
    names(facs_color) <- c("Y", "N")
    sampling_color <- c("#F4D823", "#8B40DC")
    names(sampling_color) <- c("Y", "N")
    
    if(!is.null(pc_genes)){
      class_sets <- lapply(class_sets, function(x){
        x[is.element(x, pc_genes$SYMBOL)]
      })
    }
    
    overlap_mat <- overlap_matrix(list=class_sets, jaccard)
    
    union_genes <- n_distinct(do.call(c, class_sets))
    
    meta_df <- data.frame(
      Platform=class_legend$Platform,
      FACS=class_legend$FACS_Sorted,
      `Unbiased Sampling`=class_legend$Unbiased_Sampling,
      `Set Size`=as.numeric(class_legend$SetSize),
      Study=class_legend$Study,
      check.names=F
    )
    rownames(meta_df) <- class_legend$SetID
    
    anno_colors <- list(
      Study=study_color,
      Platform=platform_color,
      FACS=facs_color,
      `Unbiased Sampling`=sampling_color
    )
    
    class_name <- cell_class_full_name(cell_class)
    
    plot_sub <- paste(
      format(union_genes, big.mark=","), "union cell type marker genes"
    )
    
    h <- 12; w <- 14; cs <- 20
    if(nrow(overlap_mat)>=40){
      h <- 20; w <- 22; cs <- 12
    } else if(nrow(overlap_mat)>=20){
      h <- 15; w <- 18
    }
    
    pdf(file=file_path, width=w, height=h)
    
    if(jaccard){
      
      overlap_mat[is.nan(overlap_mat)] <- 0
      display_mat <- apply(overlap_mat, 2, function(x) sprintf("%.0f%%", x))
      if(nrow(overlap_mat)>30){
        display_mat <- F
      }
      
      plot_title <- paste(class_name, "Marker Overlap Percent\n", plot_sub)
      
      if(!is.null(pc_genes)){
        plot_title <- paste0("Marker Overlap Percent\nRestricting to Protein Coding Genes\n", plot_sub)
      }
      
      break_list <- seq(0, 100, 10)
      leg_color <- colorRampPalette(rev(brewer.pal(9, name="RdYlBu")))(length(break_list)-1)
      
      print(
        pheatmap::pheatmap(mat=overlap_mat, main=plot_title, scale="none", cluster_rows=T, cluster_cols=T,  labels_row=class_legend$Plot_Label, labels_col=class_legend$Plot_Label, border_color="black", angle_col=45, annotation_col=meta_df, annotation_colors=anno_colors, display_numbers=display_mat, number_color="black", color=leg_color, class_legend_breaks=break_list, breaks=break_list, class_legend_labels=break_list, cellwidth=cs+1, cellheight=cs+1)
      )
      
    } else {
      
      plot_title <- paste("Number of Shared", class_name, "Markers\n", plot_sub)
      
      if(!is.null(pc_genes)){
        plot_title <- paste("Number of Shared", class_name, "Markers\nRestricting to Protein Coding Genes\n",  plot_sub)
      }
      
      if(floor(max(overlap_mat)/1e3)==0){
        
        max_val <- 1e3
        break_size <- 100
        
      } else {
        
        max_val <- round_any(max(overlap_mat), 1000, f=ceiling)
        break_size <- 200
        
      }
      
      break_list <- seq(0, max_val, by=break_size)
      leg_color <- colorRampPalette(rev(brewer.pal(9, name="RdYlBu")))(length(break_list)-1)
      
      print(
        pheatmap::pheatmap(mat=overlap_mat, main=plot_title, scale="none", angle_col=45, labels_row=class_legend$Plot_Label, labels_col=class_legend$Plot_Label, annotation_col=meta_df, annotation_colors=anno_colors, border_color="black", color=leg_color, breaks=break_list, class_legend_breaks=break_list, class_legend_labels=comma(break_list), cellwidth=cs, cellheight=cs)
      )
      
    } ## if(jaccard){} else {
    
    dev.off()
    
  } ## for(j in 1:length(cell_classes)){
  
}

marker_overlap_barplot <- function(datinfo, cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"), legend, gene_sets, jaccard=F, pc_genes, min_size=0, remove_subclust=T){
  
  datinfo <- datinfo %>% na_if("") %>% dplyr::select(-Region)
  datinfo$Study <- paste(datinfo$First_Author, datinfo$Year)
  datinfo$Study[is.element(datinfo$PubMedID, "35165441")] <- paste(datinfo$Study[is.element(datinfo$PubMedID, "35165441")], " ")
  
  legend <- legend[is.element(legend$Cell_Class, cell_classes),]
  legend <- legend[class_legend$SetSize>=min_size,]
  
  if(remove_subclust){
    legend <- legend[!grepl("_DE_", legend$SetName),]
    legend <- legend[!grepl("COMBO", legend$SetName),]
    legend <- legend[!grepl("VS_L[0-9]", legend$SetName),]
    legend <- legend[!grepl("L[0-9]_VS", legend$SetName),]
  }
  
  df <- c()
  
  for(j in 1:length(cell_classes)){
    
    cell_class <- cell_classes[j]
    print(cell_class)
    
    class_legend <- legend[is.element(legend$Cell_Class, cell_class),]
    
    if(nrow(class_legend)==0){
      next
    }
    
    class_legend <- merge(class_legend, datinfo, by.x="PubMed", by.y="PubMedID")
    class_legend <- class_legend %>% dplyr::arrange(SetName)
    class_legend <- class_legend[!duplicated(class_legend$SetName),]
    
    class_sets <- gene_sets[match(class_legend$SetName, names(gene_sets))] ## class_legend$SetName
    
    if(!identical(class_legend$SetName, names(class_sets))){
      stop("!identical(names(class_sets, class_legend$SetName))")
    }
    
    if(!is.null(pc_genes)){
      class_sets <- lapply(class_sets, function(x){
        x[is.element(x, pc_genes$SYMBOL)]
      })
    }
    
    overlap_mat <- overlap_matrix(list=class_sets, jaccard)
    diag(overlap_mat) <- NA
    
    df_class <- reshape2::melt(overlap_mat)
    colnames(df_class)[ncol(df_class)] <- "Overlap"
    df_class$Cell_Class <- cell_class
    df_class <- df_class %>% 
      dplyr::select(Cell_Class, Overlap) %>%
      tidyr::drop_na()
    
    if(is.null(df)){
      df <- df_class
    } else {
      df <- rbind(df, df_class)
    }
    
  } ## for(j in 1:length(cell_classes)){
  
  legend <- merge(legend, datinfo, by.x="PubMed", by.y="PubMedID")
  legend$Plot_Label <- paste(legend$First_Author, legend$Region, legend$Year, legend$Original_Label)
  legend <- legend %>%
    dplyr::group_by(Plot_Label) %>%
    dplyr::mutate(ID=1:n())
  legend$Plot_Label[legend$ID>1] <- paste(legend$Plot_Label[legend$ID>1], legend$ID[legend$ID>1])
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/number_of_marker_genes_overlap_barplot.pdf")
  
  if(jaccard){
    file_path <- gsub("number_of", "jaccard",  file_path)
  }
  
  pdf(file=file_path, width=8, height=6.5)
  
  if(jaccard){
    
    plot_title <- paste("Marker Gene Overlap Percent")
    
    y_max <- 100
    y_breaks <- seq(0, y_max, by=20)
    y_lab <- paste("% Overlap")
    y_labs <- comma(y_breaks)
    
  } else {
    
    plot_title <- paste("Number of Overlapping Marker Genes")
    
    y_max <- round_any(max(df$Overlap), 1000)
    y_breaks <- seq(0, y_max, by=200)
    y_lab <- ("# Overlapping Markers")
    y_labs <- comma(y_breaks)
    
  }
  plot_sub <- paste("Pairwise comparison between cell types within each class\nRestricting to protein coding genes")
  
  temp <- df %>%
    dplyr::group_by(Cell_Class) %>%
    dplyr::summarise(Mean=mean(Overlap)) %>%
    dplyr::arrange(Mean)
  
  class_name <- cell_class_full_name(cell_class)
  df$Cell_Class <- factor(df$Cell_Class, levels=temp$Cell_Class)
  
  print(
    ggplot(df, aes(x=Cell_Class, y=Overlap, group=Cell_Class)) +
      geom_boxplot(notch=F, outlier.size=0, fill="#48AAEA", color="#505050", show.legend=F) +
      #geom_violin() +
      #geom_jitter(size=.1, alpha=.1, show.legend=F) +
      theme_classic() + 
      theme(
        plot.title=element_text(hjust=.5, face="bold", size=15),
        plot.subtitle=element_text(hjust=.5, size=12, lineheight=1.1, margin=margin(t=4, b=8)),
        legend.title=element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.box="vertical",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=12, margin=margin(t=5), angle=0),
        axis.ticks.x=element_line(size=.4),
        axis.title.y=element_text(size=14, color="black", margin=margin(r=15)),
        axis.ticks.y=element_line(size=.4),
        plot.margin = unit(c(2, 2, 2, 2), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + ylab(y_lab) 
    #scale_y_continuous(trans="log")
    # +
    #   scale_y_continuous(labels=y_labs, breaks=y_breaks)
  )
  
  dev.off()
  
}

plot_marker_distro <- function(datinfo, legend){

  legend <- legend[legend$Cortical==T,]
  legend <- legend %>% na_if("") %>% na.omit()
  legend <- legend[!is.element(legend$Cell_Class, "NEU"),]
  
  if(remove_subclust){
    legend <- legend[!grepl("_DE_", legend$SetName),]
    legend <- legend[!grepl("COMBO", legend$SetName),]
    legend <- legend[!grepl("VS_L[0-9]", legend$SetName),]
    legend <- legend[!grepl("L[0-9]_VS", legend$SetName),]
  }
  
  df <- merge(legend, datinfo, by.x="PubMed", by.y="PubMedID")
  df <- df %>% dplyr::arrange(SetName)
  df <- df[!duplicated(df$SetName),]

  df$Platform[is.element(df$First_Author, "Bakken")] <- "10x Chromium V3"
  df$Platform[is.element(df$First_Author, "Luo")] <- "snmCAT-seq"
  
  pdf(file=paste0("figures/number_of_markers_boxplot_", n_distinct(legend$PubMed), "_studies.pdf"), width=10, height=8)
  
  x_var <- "Study"; x_lab <- "Study"
  x_var_type <- "discrete"
  boxplot_template(df, x_var, x_lab, x_var_type, color_var=NULL, color_lab=NULL, color_var_type=NULL)
  
  # x_var <- "Cell_Class"; x_lab <- "Cell Class"
  # x_var_type <- "discrete"
  # boxplot_template(df, x_var, x_lab, x_var_type, color_var=NULL, color_lab=NULL, color_var_type=NULL)
  
  dev.off()
  
}

boxplot_template <- function(df, x_var, x_lab, x_var_type=c("continuous", "discrete"), color_var=NULL, color_lab=NULL, color_var_type=NULL){
  
  #r2 <- summary(lm(df$SetSize~df[,x_var]))$adj.r.squared
  
  plot_title <- paste("Number of Cell Type Markers")
  plot_sub <- paste(comma(nrow(df)), "cell types from", n_distinct(df$Study), "studies") # "adj. R2 =", signif(r2, 2), "\n", 
  
  temp <- df %>%
    dplyr::group_by(!!as.name(x_var)) %>%
    dplyr::summarise(n=median(SetSize)) %>%
    dplyr::arrange(n) %>%
    as.data.frame()
  
  color_vals <- rev(brewer.pal(n_distinct(df$Platform), "Set1"))
  
  df[,x_var] <- factor(df[,x_var], levels=temp[,x_var])
  
  p <- 
    ggplot(df, aes_string(x=x_var, y="SetSize", fill="Platform")) +
    geom_boxplot(notch=F, width=.4, outlier.shape=NA, show.legend=F) +
    geom_jitter(size=1.5, height=.4, width=.3, shape=21) + # color="black",
    theme_minimal() +
    theme(
      plot.title=element_text(hjust=0, size=15, face="bold", margin=margin(b=10)),
      plot.subtitle=element_text(hjust=0, size=12, lineheight=1.1, margin=margin(b=15)),
      legend.title=element_blank(),
      legend.text=element_text(size=11),
      legend.position="bottom",
      legend.direction="horizontal",
      legend.box="vertical",
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.title.y=element_text(size=12, margin=margin(r=20)),
      axis.text.y=element_text(color="black"),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    ylab("# Marker Genes") +
    scale_y_continuous(labels=comma) +
    scale_fill_manual(values=rev(color_vals)) +
    scale_color_manual(values=rev(color_vals))
  
  if(n_elements>6){
    
    p <- p + theme(axis.text.x=element_text(color="black", size=11, hjust=1, angle=45, margin=margin(t=-28, b=35))) + 
      guides(fill=guide_legend(override.aes=list(size=6)))
    
  } else {
    p <- p + theme(axis.text.x=element_text(face="bold", size=12, margin=margin(t=-10)))
  }
  
  print(p)
  
  print(p + scale_y_continuous(trans="log2", labels=comma))
  
} ## violin_template <- function(



