library(ggplot2)
library(ggrepel)
library(nipnTK) ## outliersMD()
library(data.table)
library(dplyr) 
library(plyr) ## round_any()
library(RColorBrewer) ## brewer.pal()
library(scales)

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/continuous_plot_breaks_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/PC_scatter_template_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_data_fxn.R")

plot_PCs <- function(class_expr, datinfo, cell_class, data_type=c("author_data"), expr_type=c("counts", "QC_counts"), combat=F, pc_x, pc_y, top_n=NULL, neu_subtypes){

  if(ncol(class_expr)-1<as.numeric(gsub("[A-Z]", "", pc_y))){
    print("pc_y is > # of columns in dataset")
  } else {
    
    if(!is.null(top_n)){
      if((ncol(class_expr)-1)>as.numeric(top_n)){
        top_n <- NULL
        print("# of columns in dataset > top_n genes")
      }
    }
    
    ct_data <- prep_CT_data(
      datinfo, 
      ct_df=class_expr, 
      expr_type, 
      which_genes="intersection", 
      top_n, 
      pc_genes=T,
      ct_expr=NULL
    )
    
    datinfo <- ct_data[[1]]
    class_expr <- ct_data[[2]]
    
    if(sum(is.element(cell_class, c("EXC", "INH")))>0){
      
      neu_subtypes <- neu_subtypes %>%
        dplyr::select(
          Label, 
          Class_Level1, Class_Level2,
          Hodge_2018_Annotation, 
          Bakken_2019_Annotation
        )
      
      neu_subtypes$Label <- make.names(neu_subtypes$Label)
      datinfo <- merge(datinfo, neu_subtypes, by="Label", all.x=T)
      
    }
    
    if(is.null(top_n)){
      
      n_genes <- nrow(class_expr)
      
    } else {
      
      n_genes <- paste("top", nrow(class_expr), sep="_")
      
    }
    
    n_cts <- length(unique(paste(datinfo$Study, datinfo$Cell_Type)))
    
    pcs <- prcomp(t(class_expr[,-c(1)]), center=T, scale=F, rank.=4)
    var_exp <- pcs$sdev^2/sum(pcs$sdev^2)
    
    pc_df <- data.frame(Label=make.names(rownames(pcs$x)), Variance_Explained=var_exp, pcs$x)
    
    file_path <- paste0("figures/", data_type, "/", expr_type, "/", cell_class, "_mean_expr_", pc_x, "_vs_", pc_y, "_scatterplot_", data_type, "_", expr_type, "_", n_cts, "_CTs_", n_distinct(datinfo$Dataset), "_datasets_", n_genes, "_intersection_genes.pdf")
    
    if(combat){
      file_path <- gsub(paste0(expr_type, "_"), paste0(expr_type, "_COMBAT_"), file_path)
    }
    
    pdf(file=file_path, width=12, height=12)
    
    n_pcs <- nrow(pc_df)
    
    if(nrow(pc_df)>15){
      n_pcs <- 15
    }
    
    print(
      ggplot(pc_df[1:n_pcs,], aes(x=1:n_pcs, y=Variance_Explained)) +
        geom_line() + geom_point() +
        labs(title=paste(cell_class, "Subtypes\n", paste(sapply(unlist(strsplit(gsub("_", " ", expr_type), " ")), upper_first), collapse=" "), "Mean Expression Scree Plot"), 
             subtitle=paste(n_cts, "cell types,", n_distinct(datinfo$Dataset), "datasets,", format(nrow(class_expr), big.mark=","), "genes")) +
        theme_bw() +
        theme(
          plot.title=element_text(hjust=.5, face="bold"),
          plot.subtitle=element_text(hjust=.5),
          plot.margin = margin(5, 3, 5, 3, "cm")
        ) +
        ylab("Variance Explained") + ylim(0, 1) + xlab("PC") + 
        scale_x_continuous(breaks=seq(1, nrow(pc_df[1:n_pcs,]), by=2))
    )
    
    pc_df <- merge(pc_df, datinfo, by="Label")
    pc_df$Year <- factor(pc_df$Year, levels=sort(unique(pc_df$Year)))
    pc_df <- pc_df %>% dplyr::arrange(Region) %>% 
      dplyr::mutate(
        Outlier_Label=gsub(
          "_", " ", gsub(".", "-", paste0(Plot_Label, "\n", Cell_Type), fixed=T), fixed=T)
      )
    
    if(nrow(pc_df)>=15){
      
      ## Only label "outlier" cell types when plotting:
      
      pc_df$Outlier_Label[!outliersMD(pc_df[,pc_x], pc_df[,pc_y])] <- ""
      
      if(sum(pc_df$Outlier_Label!="")==0){
        
        pc_df <- pc_df %>% dplyr::mutate(
          Outlier_Label=gsub(
            "_", " ", gsub(".", "-", paste0(Plot_Label, "\n", Cell_Type), fixed=T), fixed=T)
        )
        pc_df$Outlier_Label[!outliersMD(pc_df[,pc_x], pc_df[,pc_y], alpha=.05)] <- ""
        
      }
      
    } ## if(nrow(pc_df)>=15){
    
    plot_title <- paste(cell_class_full_name(cell_class), "Mean Expression Projections")
    plot_sub <- paste(n_cts, "cell types found in", n_distinct(datinfo$Dataset), "datasets from", n_distinct(datinfo$Study), "studies\n", comma(nrow(class_expr)), "genes")
    
    if(combat){
      plot_title <- gsub("Mean", "Batch Corrected Mean", plot_title)
    }
    
    size_range <- c(.4, 7)
    
    color_var <- "Year"; color_lab <- "Year"
    color_var_type <- "discrete"; color_leg_rows <- 1
    size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
    size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
    size_leg_labs <- comma(size_leg_breaks)
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
    )
    
    color_var <- "Study"; color_lab <- "Study"
    color_var_type <- "discrete"; color_leg_rows <- 6
    size_var <- "Median_Unique_Genes"; size_lab <- "Median # Unique Genes"
    size_leg_breaks <- fivenum(pc_df[,size_var], na.rm=T)
    size_leg_labs <- comma(round_any(size_leg_breaks, 100))
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
    )
    
    color_var <- "Platform"; color_lab <- "Platform"
    color_var_type <- "discrete"; color_leg_rows <- 2
    size_var <- "Median_Unique_Genes"; size_lab <- "Median # Unique Genes"
    size_leg_breaks <- n_unique_genes_breaks(pc_df[,size_var])
    size_leg_labs <- comma(size_leg_breaks)
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
    )
    
    color_var <- "FACS_Sorted"; color_lab <- "FACS"
    color_var_type <- "discrete"; color_leg_rows <- 1
    size_var <- "Median_Unique_Genes"; size_lab <- "Median # Unique Genes"
    size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
    size_leg_labs <- comma(size_leg_breaks)
    pc_df$FACS_Sorted <- factor(pc_df$FACS_Sorted, levels=c("Y", "N"))
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
    )
    
    color_var <- "Unbiased_Sampling"; color_lab <- "Unbiased Sampling"
    color_var_type <- "discrete"; color_leg_rows <- 1
    size_var <- "Median_Unique_Genes"; size_lab <- "Median # Unique Genes"
    size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
    size_leg_labs <- comma(size_leg_breaks)
    pc_df$Unbiased_Sampling <- factor(pc_df$Unbiased_Sampling, levels=c("Y", "N"))
    PC_scatter_template(
      df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
      
    )
    
    if(sum(is.element(cell_class, c("EXC", "INH")))>0){
      
      color_var <- "Class_Level1"; color_lab <- "Subtype"
      color_var_type <- "discrete"; color_leg_rows <- 1
      size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
      size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
      size_leg_labs <- comma(size_leg_breaks)
      PC_scatter_template(
        df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
      )
      
      color_var <- "Class_Level2"; color_lab <- "Subtype"
      color_var_type <- "discrete"; color_leg_rows <- 1
      size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
      size_leg_breaks <- n_nuclei_breaks(pc_df[,size_var])
      size_leg_labs <- comma(size_leg_breaks)
      PC_scatter_template(
        df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
      )
      
      color_var <- "Hodge_2018_Annotation"; color_lab <- "Hodge 2018 Subtype"
      color_var_type <- "discrete"; color_leg_rows <- 1
      size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
      size_leg_breaks <- fivenum(pc_df[,size_var], na.rm=T)
      size_leg_labs <- comma(round_any(size_leg_breaks, 100))
      PC_scatter_template(
        df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
      )
      
      color_var <- "Bakken_2019_Annotation"; color_lab <- "Bakken 2019 Subtype"
      color_var_type <- "discrete"; color_leg_rows <- 3
      size_var <- "No.Nuclei"; size_lab <- "# Nuclei"
      size_leg_breaks <- fivenum(pc_df[,size_var], na.rm=T)
      size_leg_labs <- comma(round_any(size_leg_breaks, 100))
      PC_scatter_template(
        df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
      )
      
    } ## if(sum(is.element(cell_class, c("EXC", "INH")))>0){
     
    dev.off()
    
  }
  
}

