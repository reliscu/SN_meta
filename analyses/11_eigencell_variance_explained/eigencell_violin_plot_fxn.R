library(ggplot2)
library(dplyr)
library(viridis)
library(scales)

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_stats_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/continuous_plot_breaks_fxn.R")

plot_eigencell_VE <- function(ve_df, datinfo, data_type, expr_type, n_genes){

  ve_df <- merge(ve_df, datinfo, by="Dataset", all.x=T)

  ## Match cell type stats w/ datinfo:
  
  ct_stats <- prep_CT_stats(datinfo, expr_type, pc_genes=T)
  
  ct_stats$Label <- make.names(ct_stats$Label)
  
  ct_stats <- ct_stats %>% dplyr::select(
    Cell_Class, 
    Median_UMIs, 
    Median_Unique_Genes,
    Label
  )
  
  ve_df$Label <- make.names(paste0(ve_df$Dataset, "_", gsub(" ", "_", ve_df$Cell_Type)))
  
  ve_df <- merge(ve_df, ct_stats, by="Label") ## all.x=T
  ve_df$PC1_VE <- ve_df$PC1_VE*100
  
  n_cts <- length(unique(paste(ve_df$Study, ve_df$Cell_Type)))
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/eigencell_VE_violin_plot_", data_type, "_", expr_type, "_", n_cts, "_CTs_", n_distinct(ve_df$Dataset), "_datasets_", n_genes, "_intersection_genes.pdf"), width=9, height=6)
  
  x_var <- "Cell_Class"; x_lab <- "Cell Class"
  x_var_type <- "discrete"
  violin_template(ve_df, x_var, x_lab, x_var_type, n_genes, n_cts)
  
  x_var <- "Study"; x_lab <- "Study"
  x_var_type <- "discrete"
  violin_template(ve_df, x_var, x_lab, x_var_type, n_genes, n_cts)
  
  x_var <- "Platform"; x_lab <- "Platform"
  x_var_type <- "discrete"
  violin_template(ve_df, x_var, x_lab, x_var_type, n_genes, n_cts)
  
  x_var <- "Median_Unique_Genes"; x_lab <- "Median # Unique Genes per Nucleus"
  x_var_type <- "coninuous"
  color_var <- "No.Nuclei"; color_lab <- "# Nuclei"
  color_var_type <- "continuous"
  violin_template(ve_df, x_var, x_lab, x_var_type, n_genes, n_cts, color_var, color_lab, color_var_type)
  
  x_var <- "Median_Unique_Genes"; x_lab <- "Median # Unique Genes per Nucleus"
  x_var_type <- "coninuous"
  color_var <- "Platform"; color_lab <- "Platform"
  color_var_type <- "discrete"
  violin_template(ve_df, x_var, x_lab, x_var_type, n_genes, n_cts, color_var, color_lab, color_var_type)
  
  x_var <- "Median_Unique_Genes"; x_lab <- "Median # Unique Genes per Nucleus"
  x_var_type <- "coninuous"
  color_var <- "Median_UMIs"; color_lab <- "Median # UMIs per Nucleus"
  color_var_type <- "continuous"
  violin_template(ve_df, x_var, x_lab, x_var_type, n_genes, n_cts, color_var, color_lab, color_var_type)
    
  dev.off()
  
}

violin_template <- function(df, x_var, x_lab, x_var_type=c("continuous", "discrete"), n_genes, n_cts, color_var=NULL, color_lab=NULL, color_var_type=NULL){
  
  plot_title <- paste("Eigencell % Variance Explained vs.", x_lab)
  if(nchar(x_lab)>30){
    plot_title <- paste("Eigencell % Variance Explained vs.\n", x_lab)
  }

  if(x_var_type=="discrete"){
    
    plot_sub <- paste0(n_cts, " cell types from ", n_distinct(df$Study),  " studies")
    
    temp <- df %>%
      dplyr::group_by(!!as.name(x_var)) %>%
      dplyr::summarise(n=mean(PC1_VE)) %>%
      dplyr::arrange(n) %>%
      as.data.frame()
    
    df[,x_var] <- factor(df[,x_var], levels=unique(temp[,x_var]))
  
    p <- 
      ggplot(df, aes_string(x=x_var, y="PC1_VE", fill=x_var)) + 
      geom_violin(size=.5) + 
      geom_boxplot(width=.1, color="black", fill="white", outlier.shape=NA, show.legend=F) + 
      geom_jitter(color="black", alpha=.5, size=.05, show.legend=F) + # width=.3,
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=0, face="bold", size=15),
        plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=8)),
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        legend.position="bottom",
        legend.box="horizontal",
        axis.text.x=element_blank(),
        axis.ticks.x=element_line(size=0),
        axis.title.x=element_blank(), 
        axis.title.y=element_text(size=11, color="black", margin=margin(r=15)),
        plot.margin=unit(c(2, 2, 2, 2), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + 
      ylab("% Variance Explained") + 
      scale_y_continuous(limits=c(0, max(df$PC1_VE)+.5)) +
      scale_fill_manual(
        name=x_lab, values=brewer_fxn(n_distinct(df[,x_var]))
      ) +
      guides(fill="none", color="none")
    
    # if(n_distinct(df[,x_var])>10){
    #   
    #   p <- p + aes_string(color=x_var) +
    #     scale_color_manual(
    #       name=x_lab, values=brewer_fxn(n_distinct(df[,x_var]))
    #     ) 
    #   
    # } else {
    #   
    #   p <- p + aes_string(fill=x_var) +
    #     scale_fill_manual(
    #       name=x_lab, values=brewer_fxn(n_distinct(df[,x_var]))
    #     ) 
    #   
    # }
  
    if(x_var=="Cell_Class"){
      
      p <- p + theme(axis.text.x=element_text(color="black", size=11)) #+ guides(fill="none")
      
    } else if(x_var=="Platform") {
      
      p <- p + guides(fill=guide_legend(byrow=T, nrow=2, keywidth=.85, keyheight=.85, override.aes=list(size=.5))) #+ theme(plot.title=element_text(hjust=0), plot.subtitle=element_text(hjust=0))
      
    } else {
      p <- p + theme(axis.text.x=element_text(color="black", size=10, hjust=1, angle=45, margin=margin(t=-28)))
    }
    
  } else { ## if(x_var_type=="discrete"){
    
    r2 <- summary(lm(df$PC1_VE~df[,x_var]))$adj.r.squared
    
    plot_sub <- paste0("adj. R2 = ", signif(r2, 2), "\n", n_cts, " cell types from ", n_distinct(df$Study),  " studies")
    
    if(color_var_type=="discrete"){
      
      if(grepl("Genes", x_var)){
        x_breaks <- n_unique_genes_breaks(df[,x_var])
      } else if(grepl("Nuclei", x_var)){
        x_breaks <- n_nuclei_breaks(df[,color_var])
      } else {
        x_breaks <- n_UMIs_breaks(df[,color_var])
      }
      
      x_labs <- comma(x_breaks)
      
      p <- 
        ggplot(df, aes_string(x=x_var, y="PC1_VE", color=color_var)) + 
        geom_point(size=.25, alpha=.8) + 
        theme_minimal() + 
        theme(
          plot.title=element_text(hjust=.5, face="bold", size=15),
          plot.subtitle=element_text(hjust=.5, size=11, lineheight=1, margin=margin(t=4, b=8)),
          legend.title=element_text(size=11),
          legend.text=element_text(size=9),
          legend.position="bottom",
          legend.box="vertical",
          axis.title.x=element_text(size=11.5, face="bold", color="black", margin=margin(t=10, b=10)), 
          axis.text.x=element_text(angle=0),
          axis.title.y=element_text(size=11, color="black", margin=margin(r=8)),
          plot.margin=unit(c(1, 3, 1, 3), "cm")
        ) +
        labs(title=plot_title, subtitle=plot_sub) + 
        xlab(x_lab) + ylab("% Variance Explained") +
        guides(color=guide_legend(
          title=color_lab, override.aes=list(alpha=.8, size=5)
        )) +
        scale_x_continuous(breaks=x_breaks, labels=x_labs) +
        scale_color_manual(values=brewer_fxn(n_distinct(df[,color_var]))) +
        guides(color=guide_legend(
          nrow=2, override.aes=list(size=5)
        ))
      
    } else { ## if(color_var_type=="discrete"){
      
      if(grepl("Genes", x_var)){
        x_breaks <- n_unique_genes_breaks(df[,x_var])
      } else if(grepl("Nuclei", x_var)){
        x_breaks <- n_nuclei_breaks(df[,color_var])
      } else {
        x_breaks <- n_UMIs_breaks(df[,color_var])
      }
      
      x_labs <- comma(x_breaks)
      
      if(grepl("Genes", color_var)){
        color_breaks <- n_unique_genes_breaks(df[,color_var])
      } else if(grepl("Nuclei", color_var)){
        color_breaks <- n_nuclei_breaks(df[,color_var])
      } else {
        color_breaks <- n_UMIs_breaks(df[,color_var])
      }
      
      color_labs <- comma(color_breaks)

      p <- 
        ggplot(df, aes_string(x=x_var, y="PC1_VE", color=color_var)) + 
        geom_point(size=.2, alpha=.8) +
        theme_minimal() + 
        theme(
          plot.title=element_text(hjust=.5, face="bold", size=13),
          plot.subtitle=element_text(hjust=.5, size=11, lineheight=1, margin=margin(t=4, b=8)),
          legend.title=element_text(size=11),
          legend.text=element_text(size=10),
          legend.position="bottom",
          legend.box="vertical",
          axis.title.x=element_text(size=11.5, face="bold", color="black", margin=margin(t=10, b=10)), 
          axis.text.x=element_text(angle=0),
          axis.title.y=element_text(size=11, color="black", margin=margin(r=8)),
          plot.margin=unit(c(.75, 3.5, .75, 3.5), "cm")
        ) +
        labs(title=plot_title, subtitle=plot_sub) + 
        xlab(x_lab) + ylab("% Variance Explained") +
        scale_x_continuous(breaks=x_breaks, labels=x_labs) +
        #scale_x_continuous(labels=comma) +
        scale_y_continuous(breaks=seq(0, 100, 25), limits=c(0, 100)) +
        scale_color_viridis(
          guide=guide_colourbar(
            title=color_lab, 
            barwidth=20, barheight=1,
            frame.colour="black",
            label.theme=element_text(
              size=8.5, angle=40, hjust=.8
            ),
            title.position="top",
            title.hjust=.54,
          ),
          trans="log",
          limits=c(min(df[,color_var]), max(df[,color_var])), 
          labels=color_labs, 
          breaks=color_breaks
        )
      
    } ## if(color_var_type=="discrete"){} else {
    
  } ## if(x_var_type=="discrete"){} else {
  
  print(p)
 
} ## violin_template <- function(
