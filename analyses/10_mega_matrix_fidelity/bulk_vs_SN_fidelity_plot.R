library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)

source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

plot_bulk_kME_vs_SN_fidelity <- function(sn_fid, bulk_kme, data_type, expr_type, cell_classes, prop_scaled=F){
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/cell_class_bulk_vs_SN_fidelity_scatterplot_", data_type, "_", expr_type, ".pdf")
  
  if(prop_scaled){
    file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
  }
  
  pdf(file=file_path, width=10, height=10)
  
  for(j in 1:length(cell_classes)){
    
    cell_class <- cell_classes[j]
    cell_class_name <- cell_class_full_name(cell_class)
    
    class_sn_fid <- sn_fid %>% dplyr::select(
      SYMBOL, 
      paste(cell_class, "Fidelity", sep="_"), 
      paste(cell_class, "Mean_Expr_Percentile", sep="_"), 
      paste(cell_class, "Mean_UMIs", sep="_"), 
      paste(cell_class, "No.Nuclei", sep="_")
    )
    
    colnames(class_sn_fid)[-c(1)] <- paste("SN", colnames(class_sn_fid)[-c(1)], sep="_")
    
    class_bulk_fid <- bulk_fid %>% 
      na.omit() %>%
      dplyr::select(
        SYMBOL, "kME", paste(cell_class, "No.samples", sep="_")
      ) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::slice_max(
        order_by=kME, with_ties=F
      )
    
    colnames(class_bulk_fid)[2] <- paste("Bulk", colnames(class_bulk_fid)[2], sep="_")
    
    class_sn_fid$SYMBOL <- toupper(class_sn_fid$SYMBOL)
    class_bulk_fid$SYMBOL <- toupper(class_bulk_fid$SYMBOL)
    
    df <- merge(class_sn_fid, class_bulk_fid, by=1)
    
    sn_col <- intersect(grep("Fidelity", colnames(df)), grep("SN", colnames(df)))
    bulk_col <- intersect(grep("Fidelity", colnames(df)), grep("Bulk", colnames(df)))
    mean_perc_col <- grep("Mean_Expr_Percentile", colnames(df))
    umi_col <- grep("Mean_UMIs", colnames(df))
    
    df$Mean_UMIs_Status <- df[,umi_col]>=1
    df$Mean_UMIs_Status <- factor(df$Mean_UMIs_Status, levels=c(T, F))
    
    r2 <- signif(summary(lm(df[,sn_col]~df[,bulk_col]))$adj.r.squared, 2)
    
    color_breaks <- seq(0, 100, by=20)
    color_breaks[which.min(color_breaks)] <- ceiling(min(df[,mean_perc_col]))
    color_breaks[which.max(color_breaks)] <- floor(max(df[,mean_perc_col]))
    color_labs <- round_any(color_breaks, 5)
    
    plot_title <- paste(cell_class_name, "Bulk vs. Single-Nucleus Fidelity")
    plot_sub <- paste0("adj. R2 = ", r2, "\n", format(nrow(df), big.mark=","), " genes") 
    leg_title <- paste(cell_class_name, "Single-Nucleus Mean Expression %tile")
    x_lab <- paste0(cell_class_name, " Bulk Fidelity\n(", format(max(df[,grep("No.samples", colnames(df))]), big.mark=","), " samples)")
    y_lab <- paste0(cell_class_name, " Single-Nucleus Fidelity\n(", format(max(df[,grep("No.Nuclei", colnames(df))]), big.mark=","), " nuclei)")
    if(prop_scaled){
      y_lab <- gsub("Fidelity", "Proportional Fidelity", y_lab)
    }
    
    p <- ggplot(df, aes_string(
      x=colnames(df)[bulk_col], y=colnames(df)[sn_col], color=colnames(df)[mean_perc_col]
    )) + 
      geom_point(alpha=.6, size=1) +
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=.5, face="bold", size=16),
        plot.subtitle=element_text(hjust=.5, size=13, lineheight=1.1, margin=margin(t=4, b=8)),
        axis.title.x=element_text(face="bold", size=13, margin=margin(t=8, b=8)),
        axis.title.y=element_text(face="bold", size=13, margin=margin(r=12)),
        legend.title=element_text(size=11, lineheight=1, margin=margin(b=3)),
        legend.position="bottom",
        legend.box="vertical",
        plot.margin=unit(c(2, 2, 2, 2), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab(x_lab) + ylab(y_lab) + 
      scale_color_viridis(
        guide=guide_colourbar(
          title=leg_title, 
          frame.colour="black",
          barwidth=20, barheight=1.5,
          label.theme=element_text(
            size=10, angle=0, hjust=.8
          ),
          title.position="top",
          title.hjust=.53,
        ),
        limits=c(0, 100), 
        labels=color_labs, 
        breaks=color_breaks
      )
    
    if(prop_scaled){
      print(p)
    } else {
      
      print(
        p + scale_y_continuous(
          labels=signif_fxn,
          limits=c(0, 1),
          breaks=seq(0, 1, by=.25)
        )
      )
      
    } ## if(prop_scaled){} else {
    
    p <- ggplot(df, aes_string(
      x=colnames(df)[bulk_col], y=colnames(df)[sn_col], color="Mean_UMIs_Status"
    )) + 
      geom_point(alpha=.3, size=1) +
      geom_point(data=df[df$Mean_UMIs_Status==F,], alpha=.3, size=1) +
      theme_minimal() + 
      theme(plot.title=element_text(hjust=.5, face="bold", size=16),
            plot.subtitle=element_text(hjust=.5, size=13, lineheight=1.1, margin=margin(t=4, b=8)),
            axis.title.x=element_text(face="bold", size=13, margin=margin(t=8, b=8)),
            axis.title.y=element_text(face="bold", size=13, margin=margin(r=12)),
            legend.title=element_text(size=12, margin=margin(t=3, b=3)),
            legend.text=element_text(size=11),
            legend.position="bottom",
            legend.box="vertical",
            plot.margin = unit(c(2, 2, 2, 2), "cm")) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab(x_lab) + ylab(y_lab) + 
      guides(color=guide_legend(
        title=expression("Mean UMIs ">="1"),
        title.position="top",
        title.hjust=.5,
        override.aes=list(size=5, alpha=1)
      )) +
      scale_color_manual(values=c("#3492EA", "#D92F4B"))
    
    if(prop_scaled){
      print(p)
    } else {
      
      print(
        p + scale_y_continuous(
          labels=signif_fxn,
          limits=c(0, 1),
          breaks=seq(0, 1, by=.25)
        )
      )
      
    } ## if(prop_scaled){} else {
    
  } ##  for(j in 1:length(cell_classes)){
  
  dev.off()
  
} ## plot_bulk_vs_SN_fidelity <- function(

plot_bulk_vs_SN_fidelity <- function(sn_fid, bulk_fid, data_type, expr_type, cell_classes, prop_scaled=F){
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/cell_class_bulk_vs_SN_fidelity_scatterplot_", data_type, "_", expr_type, ".pdf")

  if(prop_scaled){
    file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
  }
  
  pdf(file=file_path, width=10, height=10)

  for(j in 1:length(cell_classes)){
    
    cell_class <- cell_classes[j]
    cell_class_name <- cell_class_full_name(cell_class)
    
    class_sn_fid <- sn_fid %>% dplyr::select(
      SYMBOL, 
      paste(cell_class, "Fidelity", sep="_"), 
      paste(cell_class, "Mean_Expr_Percentile", sep="_"), 
      paste(cell_class, "Mean_UMIs", sep="_"), 
      paste(cell_class, "No.Nuclei", sep="_")
    )
    
    colnames(class_sn_fid)[-c(1)] <- paste("SN", colnames(class_sn_fid)[-c(1)], sep="_")
    
    class_bulk_fid <- bulk_fid %>% 
      na.omit() %>%
      dplyr::select(
        SYMBOL, paste(cell_class, "Fidelity", sep="_"), paste(cell_class, "No.samples", sep="_")
      ) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::slice_max(
        order_by=!!as.name(
          colnames(class_bulk_fid)[grep("Fidelity", colnames(class_bulk_fid))]
        ),
        with_ties=F
      )
    
    colnames(class_bulk_fid)[2] <- paste("Bulk", colnames(class_bulk_fid)[2], sep="_")
    
    df <- merge(class_sn_fid, class_bulk_fid, by=1)
    
    sn_col <- intersect(grep("Fidelity", colnames(df)), grep("SN", colnames(df)))
    bulk_col <- intersect(grep("Fidelity", colnames(df)), grep("Bulk", colnames(df)))
    mean_perc_col <- grep("Mean_Expr_Percentile", colnames(df))
    umi_col <- grep("Mean_UMIs", colnames(df))
  
    df$Mean_UMIs_Status <- df[,umi_col]>=1
    df$Mean_UMIs_Status <- factor(df$Mean_UMIs_Status, levels=c(T, F))
    
    r2 <- signif(summary(lm(df[,sn_col]~df[,bulk_col]))$adj.r.squared, 2)
    
    color_breaks <- seq(0, 100, by=20)
    color_breaks[which.min(color_breaks)] <- ceiling(min(df[,mean_perc_col]))
    color_breaks[which.max(color_breaks)] <- floor(max(df[,mean_perc_col]))
    color_labs <- round_any(color_breaks, 5)
    
    plot_title <- paste(cell_class_name, "Bulk vs. Single-Nucleus Fidelity")
    plot_sub <- paste0("adj. R2 = ", r2, "\n", format(nrow(df), big.mark=","), " genes") 
    leg_title <- paste(cell_class_name, "Single-Nucleus Mean Expression %tile")
    x_lab <- paste0(cell_class_name, " Bulk Fidelity\n(", format(max(df[,grep("No.samples", colnames(df))]), big.mark=","), " samples)")
    y_lab <- paste0(cell_class_name, " Single-Nucleus Fidelity\n(", format(max(df[,grep("No.Nuclei", colnames(df))]), big.mark=","), " nuclei)")
    if(prop_scaled){
      y_lab <- gsub("Fidelity", "Proportional Fidelity", y_lab)
    }
    
    p <- ggplot(df, aes_string(
      x=colnames(df)[bulk_col], y=colnames(df)[sn_col], color=colnames(df)[mean_perc_col]
    )) + 
      geom_point(alpha=.6, size=1) +
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=.5, face="bold", size=16),
        plot.subtitle=element_text(hjust=.5, size=13, lineheight=1.1, margin=margin(t=4, b=8)),
        axis.title.x=element_text(face="bold", size=13, margin=margin(t=8, b=8)),
        axis.title.y=element_text(face="bold", size=13, margin=margin(r=12)),
        legend.title=element_text(size=11, lineheight=1, margin=margin(b=3)),
        legend.position="bottom",
        legend.box="vertical",
        plot.margin=unit(c(2, 2, 2, 2), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab(x_lab) + ylab(y_lab) + 
      scale_color_viridis(
        guide=guide_colourbar(
          title=leg_title, 
          barwidth=20, barheight=1.5,
          label.theme=element_text(
            size=10, angle=0, hjust=.8
          ),
          title.position="top",
          title.hjust=.53,
        ),
        limits=c(0, 100), 
        labels=color_labs, 
        breaks=color_breaks
      )
    
    if(prop_scaled){
      print(p)
    } else {
      
      print(
        p + scale_y_continuous(
          labels=signif_fxn,
          limits=c(0, 1),
          breaks=seq(0, 1, by=.25)
        )
      )
      
    } ## if(prop_scaled){} else {
    
    p <- ggplot(df, aes_string(
      x=colnames(df)[bulk_col], y=colnames(df)[sn_col], color="Mean_UMIs_Status"
    )) + 
      geom_point(alpha=.3, size=1) +
      geom_point(data=df[df$Mean_UMIs_Status==F,], alpha=.3, size=1) +
      theme_minimal() + 
      theme(plot.title=element_text(hjust=.5, face="bold", size=16),
            plot.subtitle=element_text(hjust=.5, size=13, lineheight=1.1, margin=margin(t=4, b=8)),
            axis.title.x=element_text(face="bold", size=13, margin=margin(t=8, b=8)),
            axis.title.y=element_text(face="bold", size=13, margin=margin(r=12)),
            legend.title=element_text(size=12, margin=margin(t=3, b=3)),
            legend.text=element_text(size=11),
            legend.position="bottom",
            legend.box="vertical",
            plot.margin = unit(c(2, 2, 2, 2), "cm")) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab(x_lab) + ylab(y_lab) + 
      guides(color=guide_legend(
        title=expression("Mean UMIs ">="1"),
        title.position="top",
        title.hjust=.5,
        override.aes=list(size=5, alpha=1)
      )) +
      scale_color_manual(values=c("#3492EA", "#D92F4B"))
    
    if(prop_scaled){
      print(p)
    } else {
      
      print(
        p + scale_y_continuous(
          labels=signif_fxn,
          limits=c(0, 1),
          breaks=seq(0, 1, by=.25)
        )
      )
      
    } ## if(prop_scaled){} else {
  
  } ##  for(j in 1:length(cell_classes)){
  
  dev.off()
  
} ## plot_bulk_vs_SN_fidelity <- function(

signif_fxn <- function(x) sprintf("%.2f", x)
