library(ggplot2)
library(ggrepel)
library(data.table)
library(plyr) ## round_any()
library(dplyr)
library(scales) ## comma
library(flexiblas)

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/PC_scatter_template_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/continuous_plot_breaks_fxn.R")

plot_PCs <- function(dataset_expr, datinfo, data_type, expr_type=c("counts", "normalized_counts"), combat=F, gene_list=NULL, top_n=NULL, pc_x, pc_y){
  
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
  
  pcs <- prcomp(t(dataset_expr[,-c(1)]), center=T, scale=T, rank.=4)
  var_exp <- pcs$sdev^2/sum(pcs$sdev^2)
  
  pc_df <- data.frame(Dataset=rownames(pcs$x), Variance_Explained=var_exp, pcs$x)

  file_path <- paste0("figures/", data_type, "/", expr_type, "/dataset_mean_expr_", pc_x, "_vs_", pc_y, "_scatterplot_", data_type, "_", expr_type, "_", nrow(pc_df), "_datasets_", n_distinct(datinfo$Study), "_studies_", n_genes, "_intersection_genes.pdf")
  
  if(combat){
    file_path <- gsub(paste0(expr_type, "_"), paste0(expr_type, "_COMBAT_"), file_path)
  }
  
  pdf(file=file_path, width=10, height=10)
  
  ## Scree plot
  print(
    ggplot(pc_df, aes(x=1:length(var_exp), y=Variance_Explained)) +
      geom_line() + geom_point() + 
      labs(title=paste(paste(sapply(unlist(strsplit(gsub("_", " ", expr_type), " ")), upper_first), collapse=" "), "Mean Expression Scree Plot"), subtitle=paste(nrow(pc_df), "datasets,", format(nrow(dataset_expr), big.mark=","), "genes")) +
      theme_bw() + 
      theme(
        plot.title=element_text(hjust=.5),
        plot.subtitle=element_text(hjust=.5),
        plot.margin=margin(5, 3, 5, 3, "cm")
      ) +
      ylab("Variance Explained") + ylim(0, 1) + xlab("PC") + 
      scale_x_continuous(breaks=seq(1, nrow(pc_df), by=2))
  )
  
  pc_df <- merge(pc_df, datinfo, by="Dataset")
  pc_df$Outlier_Label <- pc_df$Plot_Label
  
  plot_title <- paste("Dataset Mean Expression Projections")
  plot_sub <- paste(comma(n_genes), "protein coding genes")
  
  size_range <- c(.4, 8)
  
  color_var <- "Study"; color_lab <- "Study"; 
  color_var_type="discrete"; color_leg_rows <- 4
  size_var <- "Author_Median_Unique_Genes_PC"; size_lab <- "Median # Unique Genes per Nucleus"
  size_leg_breaks <- n_unique_genes_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks=NULL, color_leg_labs=NULL, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
  )
  
  color_var <- "Platform"; color_lab <- "Platform"; 
  color_var_type="discrete"; color_leg_rows <- 2
  size_var <- "Author_Median_Unique_Genes_PC"; size_lab <- "Median # Unique Genes per Nucleus"
  size_leg_breaks <- n_unique_genes_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks=NULL, color_leg_labs=NULL, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
  )
  
  color_var <- "Author_Median_UMIs_PC"; color_lab <- "Median # UMIs per Nucleus"
  color_var_type="continuous"; color_leg_rows <- 1
  color_leg_breaks <- n_UMIs_breaks(pc_df[,color_var])
  color_leg_breaks[which.min(color_leg_breaks)] <- 505
  color_leg_labs <- comma(color_leg_breaks)
  color_leg_labs[1] <- "500"
  size_var <- "Author_Median_Unique_Genes_PC"; size_lab <- "Median # Unique Genes per Nucleus"
  size_leg_breaks <- n_unique_genes_breaks(pc_df[,size_var])
  size_leg_labs <- comma(size_leg_breaks)
  PC_scatter_template(
    df=pc_df, x=pc_x, y=pc_y, plot_title, plot_sub, color_var, color_lab, color_var_type, color_leg_rows, color_leg_breaks, color_leg_labs, size_var, size_lab, size_leg_breaks, size_leg_labs, size_range, label_outliers=T
  )
  
  dev.off()
  
}

PC_vs_covariates <- function(dataset_expr, datinfo, data_type=c("author_data"), expr_type=c("counts", "normalized_counts"), combat=F, gene_list=NULL, top_n=NULL, which_pc){
  
  datinfo <- datinfo %>% na_if("") %>% dplyr::filter(is.element(Dataset, colnames(dataset_expr)))
  dataset_expr <- dataset_expr[,c(1, match(datinfo$Dataset, colnames(dataset_expr)))]
  
  ## Remove genes that are NA in ANY datasets:
  
  dataset_expr <- dataset_expr[apply(dataset_expr[,-c(1)], 1, function(x) sum(is.na(x)))==0,]
  
  n_genes <- nrow(dataset_expr)
  
  if(!is.null(top_n)){
    mean_dataset_expr <- apply(dataset_expr[,-c(1)], 1, function(x) mean(x, na.rm=T))
    dataset_expr <- dataset_expr[order(-mean_dataset_expr),] %>% dplyr::slice(1:top_n)
    n_genes <- paste("top", nrow(dataset_expr), sep="_")
  } 
  
  pcs <- prcomp(t(dataset_expr[,-c(1)]), center=T, scale=T, rank.=4)
  pc_df <- data.frame(Dataset=rownames(pcs$x), pcs$x)
  pc_df <- merge(pc_df, datinfo, by="Dataset")
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/mean_expr_covariates_vs_", which_pc, "_scatterplot_", data_type, "_", expr_type, "_", ncol(dataset_expr)-1, "_datasets_", n_genes, "_intersection_genes.pdf")
  
  if(combat){
    file_path <- gsub(expr_type, paste0(expr_type, "_COMBAT"), file_path)
  }
  
  pdf(file=file_path, width=10, height=10)
  
  covar <- "Median # UMIs per Nucleus"
  model <- summary(lm(pc_df[,which_pc]~pc_df$Author_Median_UMIs_PC))
  plot(pc_df[,which_pc], pc_df$Author_Median_UMIs, main=paste(covar, "vs.", which_pc), xlab=which_pc, ylab=covar, sub=paste0("adj. R2=", signif(model$adj.r.squared, 2)))
  
  covar <- "Median # Unique Genes per Nucleus"
  model <- summary(lm(pc_df[,which_pc]~pc_df$Author_Median_Unique_Genes_PC))
  plot(pc_df[,which_pc], pc_df$Author_Median_Unique_Genes, main=paste(covar, "vs.", which_pc), xlab=which_pc, ylab=covar, sub=paste0("adj. R2=", signif(model$adj.r.squared, 2)))
  
  covar <- "# Nuclei"
  model <- summary(lm(pc_df[,which_pc]~pc_df$Author_No.Nuclei))
  plot(pc_df[,which_pc], pc_df$Author_No.Nuclei, main=paste(covar, "vs.", which_pc), xlab=which_pc, ylab=covar, sub=paste0("adj. R2=", signif(model$adj.r.squared, 2)))
  
  covar <- "# Cell Types"
  model <- summary(lm(pc_df[,which_pc]~pc_df$Author_No.Cell_Types))
  plot(pc_df[,which_pc], pc_df$Author_No.Cell_Types, main=paste(covar, "vs.", which_pc), xlab=which_pc, ylab=covar, sub=paste0("adj. R2=", signif(model$adj.r.squared, 2)))
  
  covar <- "FACS sorted"
  model <- summary(lm(pc_df[,which_pc]~as.factor(pc_df$FACS_Sorted)))
  plot(pc_df[,which_pc], as.numeric(as.factor(pc_df$FACS_Sorted)), main=paste(covar, "vs.", which_pc), xlab=which_pc, ylab=covar, sub=paste0("adj. R2=", signif(model$adj.r.squared, 2)))
  
  dev.off()
  
}
