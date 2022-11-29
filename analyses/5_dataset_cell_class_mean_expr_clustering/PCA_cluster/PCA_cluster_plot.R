library(ggplot2)
library(ggrepel)
library(data.table)
library(dplyr)
library(plyr) ## round_any()
library(scales) ## comma
library(nipnTK) ## outliersMD()
library(RColorBrewer) ## brewer.pal()
library(flexiblas)

flexiblas_switch(invisible(flexiblas_load_backend("NETLIB")))

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/continuous_plot_breaks_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/PC_scatter_template_fxn.R")

plot_PCs <- function(dataset_class_expr, datinfo, cell_class, data_type=c("author_data"), expr_type=c("QC_counts", "normalized_counts"), pc_x, pc_y, gene_list=NULL, top_n=NULL){
  
  if(ncol(dataset_class_expr)-1<as.numeric(gsub("PC", "", pc_y))){
    print("pc_y is > # of columns in dataset")
  } else {
    
    datinfo <- datinfo[match(colnames(dataset_expr)[-c(1)], datinfo$Dataset),]
    
    if(!identical(datinfo$Dataset, colnames(dataset_expr)[-c(1)])){
      stop("!identical(datinfo$Dataset, colnames(dataset_expr)[-c(1)])")
    }
    
    n_genes <- nrow(dataset_expr)
    
    if(!is.null(top_n)){
      
      ## Select genes with highest DATASET CELL CLASS mean expression:
      
      mean_dataset_class_expr <- rowMeans(dataset_class_expr[,-c(1)])
      dataset_class_expr <- dataset_class_expr[order(-mean_dataset_class_expr),] %>% dplyr::slice(1:top_n)
      n_genes <- paste("top", nrow(dataset_class_expr), sep="_")
      
    } 
    
    pcs <- prcomp(t(dataset_class_expr[,-c(1)]), center=T, scale=T, rank.=4)
    var_exp <- pcs$sdev^2/sum(pcs$sdev^2)
    
    pc_df <- data.frame(Dataset=rownames(pcs$x), Variance_Explained=var_exp, pcs$x)
    
    pdf(file=paste0("figures/", data_type, "/", expr_type, "/", cell_class, "_dataset_mean_expr_", pc_x, "_vs_", pc_y, "_scatterplot_", data_type, "_", expr_type, "_", ncol(dataset_class_expr)-1, "_datasets_", n_genes, "_intersection_genes.pdf"), width=10, height=10)
    
    ## Scree plot
    print(
      ggplot(pc_df, aes(x=1:length(var_exp), y=Variance_Explained)) +
        geom_line() + geom_point() + 
        labs(title=paste(cell_class, "Subtypes\n", paste(sapply(unlist(strsplit(gsub("_", " ", expr_type), " ")), upper_first), collapse=" "), "Mean Expression Scree Plot"), subtitle=paste(nrow(pc_df), "datasets,", comma(nrow(dataset_class_expr)), "genes")) +
        theme_bw() + 
        theme(plot.title=element_text(hjust=.5),
              plot.subtitle=element_text(hjust=.5),
              plot.margin = margin(5, 3, 5, 3, "cm")) +
        ylab("Variance Explained") + ylim(0, 1) +
        xlab("PC") + scale_x_continuous(breaks=seq(1, nrow(pc_df), by=2))
    )
    
    pc_df <- merge(pc_df, datinfo, by="Dataset")
    pc_df$Year <- as.character(pc_df$Year)
    
    plot_title <- paste("Dataset", cell_class_full_name(cell_class), "Mean Expression Projections")
    plot_sub <- paste(comma(n_genes), "protein coding genes")
    
    pc_df$Outlier_Label <- pc_df$Plot_Label
    
    if(nrow(pc_df)>=15){
      
      ## Only label "outlier" cell types when plotting:
      
      pc_df$Outlier_Label[!outliersMD(pc_df[,pc_x], pc_df[,pc_y])] <- ""
      
      if(sum(pc_df$Outlier_Label!="")==0){
        
        pc_df$Outlier_Label <- pc_df$Plot_Label
        pc_df$Outlier_Label[!outliersMD(pc_df[,pc_x], pc_df[,pc_y], alpha=.05)] <- ""
        
      }
      
    } ## if(nrow(pc_df)>=15){
    
    size_range <- c(.4, 8)
    
    color_var <- "Year"; color_lab <- "Year"
    color_var_type <- "discrete"; color_leg_rows <- 1
    size_var <- "Author_No.Nuclei"; size_lab <- "# Nuclei"
    size_leg_breaks <- fivenum(pc_df[,size_var], na.rm=T)
    size_leg_labs <- comma(round_any(size_leg_breaks, 100))
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
    )
    
    color_var <- "Study"; color_lab <- "Study"
    color_var_type <- "discrete"; color_leg_rows <- 6
    size_var <- "Author_Median_Unique_Genes_PC"; size_lab <- "Median # Unique Genes"
    size_leg_breaks <- fivenum(pc_df[,size_var], na.rm=T)
    size_leg_labs <- comma(round_any(size_leg_breaks, 100))
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
    )
    
    color_var <- "Platform"; color_lab <- "Platform"
    color_var_type <- "discrete"; color_leg_rows <- 2
    size_var <- "Author_Median_Unique_Genes_PC"; size_lab <- "Median # Unique Genes"
    size_leg_breaks <- n_unique_genes_breaks(pc_df[,size_var])
    size_leg_labs <- comma(size_leg_breaks)
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
    )
    
    color_var <- "FACS_Sorted"; color_lab <- "FACS"
    color_var_type <- "discrete"; color_leg_rows <- 1
    size_var <- "Author_Median_Unique_Genes_PC"; size_lab <- "Median # Unique Genes"
    size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
    size_leg_labs <- comma(size_leg_breaks)
    pc_df$FACS_Sorted <- factor(pc_df$FACS_Sorted, levels=c("Y", "N"))
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
    )
    
    color_var <- "Unbiased_Sampling"; color_lab <- "Unbiased Sampling"
    color_var_type <- "discrete"; color_leg_rows <- 1
    size_var <- "Author_Median_Unique_Genes_PC"; size_lab <- "Median # Unique Genes"
    size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
    size_leg_labs <- comma(size_leg_breaks)
    pc_df$Unbiased_Sampling <- factor(pc_df$Unbiased_Sampling, levels=c("Y", "N"))
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
      
    )
    
    dev.off()
    
  }
  
}

