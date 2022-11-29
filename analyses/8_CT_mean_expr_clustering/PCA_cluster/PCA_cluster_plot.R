library(ggplot2)
library(ggrepel)
library(nipnTK) ## outliersMD()
library(data.table)
library(plyr) ## round_any()
library(dplyr) 
library(RColorBrewer) ## brewer.pal()
library(scales)
library(flexiblas)

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_data_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/PC_scatter_template_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/continuous_plot_breaks_fxn.R")

plot_PCs <- function(
  ct_expr, 
  datinfo, 
  data_type=c("author_data"), 
  expr_type=c("counts", "QC_counts"), 
  pc_x, pc_y, 
  gene_list=NULL,
  top_n=NULL,
  neu_subtypes
  ){
  
  ct_data <- prep_CT_data(
    datinfo, 
    ct_df=ct_expr, 
    expr_type, 
    which_genes="intersection", 
    pc_genes=T, 
    top_n, 
    ct_expr=NULL
  )
  
  datinfo <- ct_data[[1]]
  ct_expr <- ct_data[[2]]
  
  if(!is.null(gene_list)){
    ct_expr <- ct_expr[is.element(ct_expr$SYMBOL, gene_list),]
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
    
    n_genes <- nrow(ct_expr)
    
  } else {
    
    n_genes <- paste("top", nrow(ct_expr), sep="_")
    
  }
  
  n_cts <- length(unique(paste(datinfo$Study, datinfo$Cell_Type)))
  
  pcs <- prcomp(t(ct_expr[,-c(1)]), center=T, scale=F, rank.=4)
  var_exp <- pcs$sdev^2/sum(pcs$sdev^2)
  
  pc_df <- data.frame(Label=make.names(rownames(pcs$x)), Variance_Explained=var_exp, pcs$x)
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/CT_mean_expr_", pc_x, "_vs_", pc_y, "_scatterplot_", data_type, "_", expr_type, "_", n_cts, "_CTs_", n_distinct(datinfo$Dataset), "_datasets_", n_genes, "_intersection_genes.pdf"), width=10, height=10)
  
  n_pcs <- nrow(pc_df)
  
  if(nrow(pc_df)>15){
    n_pcs <- 15
  }
  
  print(
    ggplot(pc_df[1:n_pcs,], aes(x=1:n_pcs, y=Variance_Explained)) +
      geom_line() + geom_point() +
      labs(title=paste("Cell Type", paste(sapply(unlist(strsplit(gsub("_", " ", expr_type), " ")), upper_first), collapse=" "), "Mean Expression Scree Plot"), 
           subtitle=paste(n_cts, "cell types,", n_distinct(datinfo$Dataset), "datasets,", format(nrow(ct_expr), big.mark=","), "genes")) +
      theme_bw() +
      theme(plot.title=element_text(hjust=.5, face="bold"),
            plot.subtitle=element_text(hjust=.5),
            plot.margin=margin(5, 3, 5, 3, "cm")) +
      ylab("Variance Explained") + ylim(0, 1) + xlab("PC") + 
      scale_x_continuous(breaks=seq(1, n_pcs, by=2))
  )
  
  pc_df <- merge(pc_df, datinfo, by="Label")
  pc_df$Year <- factor(pc_df$Year, levels=sort(unique(pc_df$Year)))
  
  plot_title <- paste("Cell Type Mean Expression")
  plot_sub <- paste(n_cts, "cell types from", n_distinct(pc_df$Study), "studies\n", comma(nrow(ct_expr)), "protein coding genes")
  
  size_range <- c(.4, 7)
  
  color_var <- "Cell_Class"; color_lab <- "Cell Class"
  color_var_type <- "discrete"; color_leg_rows <- 3
  size_var <- "Median_Unique_Genes"; size_lab <- "Median # Unique Genes per Cell Type"
  size_leg_breaks <- n_unique_genes_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks=NULL, color_leg_labs=NULL, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "Cell_Class"; color_lab <- "Cell Class"
  color_var_type <- "discrete"; color_leg_rows <- 3
  size_var <- "Median_UMIs"; size_lab <- "Median # UMIs per Cell Type"
  size_leg_breaks <- n_UMIs_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks=NULL, color_leg_labs=NULL, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "No.Nuclei"; color_lab <- "# Nuclei"; color_var_type <- "continuous"
  color_leg_breaks <- n_nuclei_breaks(pc_df[,color_var])
  color_leg_labs <- comma(color_leg_breaks)
  size_var <- "Median_UMIs"; size_lab <- "Median # UMIs per Cell Type"
  size_leg_breaks <- n_UMIs_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows=NULL, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "Year"; color_lab <- "Year"
  color_var_type <- "discrete"; color_leg_rows <- 1
  size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
  size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "Study"; color_lab <- "Study"
  color_var_type <- "discrete"; color_leg_rows <- 6
  size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
  size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "Study"; color_lab <- "Study"
  color_var_type <- "discrete"; color_leg_rows <- 6
  size_var <- "Median_Unique_Genes"; size_lab <- "Median # Unique Genes per Cell Type"
  size_leg_breaks <- n_unique_genes_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "Median_Unique_Genes"; color_lab <- "Median # Unique Genes per Cell Type"
  color_var_type <- "continous"
  color_leg_breaks <- n_unique_genes_breaks(pc_df[,color_var])
  color_leg_labs <- comma(color_leg_breaks)
  size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
  size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
 
  color_var <- "Platform"; color_lab <- "Platform"
  color_var_type <- "discrete"; color_leg_rows <- 2
  size_var <- "Median_Unique_Genes"; size_lab <- "Median # Unique Genes per Cell Type"
  size_leg_breaks <- n_unique_genes_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "FACS_Sorted"; color_lab <- "FACS"
  color_var_type <- "discrete"; color_leg_rows <- 1
  size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
  size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  pc_df$FACS_Sorted <- factor(pc_df$FACS_Sorted, levels=c("Y", "N"))
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "Unbiased_Sampling"; color_lab <- "Unbiased Sampling"
  color_var_type <- "discrete"; color_leg_rows <- 1
  size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
  size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  pc_df$Unbiased_Sampling <- factor(pc_df$Unbiased_Sampling, levels=c("Y", "N"))
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
    )
  
  color_var <- "Class_Level1"; color_lab <- "Subtype"
  color_var_type <- "discrete"; color_leg_rows <- 1
  size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
  size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "Class_Level2"; color_lab <- "Subtype"
  color_var_type <- "discrete"; color_leg_rows <- 3
  size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
  size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "Hodge_2018_Annotation"; color_lab <- "Hodge 2018 Subtype"
  color_var_type <- "discrete"; color_leg_rows <- 1
  size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
  size_leg_breaks <- fivenum(pc_df[,size_var], na.rm=T)
  size_leg_labs <- comma(round_any(size_leg_breaks, 100))
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  color_var <- "Bakken_2019_Annotation"; color_lab <- "Bakken 2019 Subtype"
  color_var_type <- "discrete"; color_leg_rows <- 3
  size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
  size_leg_breaks <- fivenum(pc_df[,size_var], na.rm=T)
  size_leg_labs <- comma(round_any(size_leg_breaks, 100))
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range
  )
  
  dev.off()
  
}

