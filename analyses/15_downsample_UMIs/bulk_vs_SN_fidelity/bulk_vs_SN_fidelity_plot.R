library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggtext)
library(scales)
library(viridis)

source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

plot_bulk_vs_SN_fidelity <- function(datinfo, class_sn_fid, bulk_fid, cell_class=c("NEU", "ASC", "MIC", "OG", "OPC"), data_type, expr_type, prop_scaled=F){

  cell_class_name <- cell_class_full_name(cell_class)
  
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

  class_sn_fid$Threshold <- as.numeric(class_sn_fid$Threshold)
  class_sn_fid <- class_sn_fid[!is.element(class_sn_fid$Threshold, 1e6),]

  idx <- class_sn_fid$Threshold%%100>0
  pre_umis <- mean(class_sn_fid$Threshold[idx])
  class_sn_fid$Threshold[idx] <- 1e6
  class_sn_fid <- class_sn_fid[class_sn_fid$Threshold%in%c(500, 1e3, 1e4, 1e5, 1e6),]

  idx <- order(unique(class_sn_fid$Threshold))
  class_sn_fid$Threshold <- sapply(class_sn_fid$Threshold, comma)
  class_sn_fid$Threshold[class_sn_fid$Threshold=="1,000,000"] <- "Full coverage"
  class_sn_fid$Threshold[class_sn_fid$Threshold!="Full coverage"] <- paste(class_sn_fid$Threshold[class_sn_fid$Threshold!="Full coverage"], "UMIs")
  class_sn_fid$Threshold <- factor(
    class_sn_fid$Threshold, levels=unique(class_sn_fid$Threshold)[idx]
  )
  
  class_sn_fid$SYMBOL <- toupper(class_sn_fid$SYMBOL)
  class_bulk_fid$SYMBOL <- toupper(class_bulk_fid$SYMBOL)
  
  thresholds <- unique(class_sn_fid$Threshold)

  df_list <- lapply(1:length(thresholds), function(k){
    
    temp <- class_sn_fid[is.element(class_sn_fid$Threshold, thresholds[k]),]
    colnames(temp)[-c(1)] <- paste("SN", colnames(temp)[-c(1)], sep="_")
    df <- merge(temp, class_bulk_fid, by=1)
    
    return(df)
    
  })
  df <- do.call(rbind, df_list)
  
  sn_col <- intersect(grep("Fidelity", colnames(df)), grep("SN", colnames(df)))
  bulk_col <- intersect(grep("Fidelity", colnames(df)), grep("Bulk", colnames(df)))
  mean_perc_col <- grep("Mean_Expr_Percentile", colnames(df))
  #umi_col <- grep("Mean_UMIs", colnames(df))
  
  color_breaks <- seq(0, 100, by=20)
  color_breaks[which.min(color_breaks)] <- ceiling(min(df[,mean_perc_col]))
  color_breaks[which.max(color_breaks)] <- floor(max(df[,mean_perc_col]))
  color_labs <- round_any(color_breaks, 5)
  
  plot_title <- paste(cell_class_name, "Single-Nucleus Fidelity vs. Read Depth per Nucleus")
  plot_sub <- paste0(
    "Average median # UMIs per nucleus prior to downsampling: ", comma(pre_umis), "\n", comma(n_distinct(df$SYMBOL)), " protein coding genes"
  ) 
   
  leg_title <- paste("Mean Expression %tile") # cell_class_name, "Single-Nucleus 
  x_lab <- paste0("Bulk Fidelity\n(", comma(max(df[,grep("No.samples", colnames(df))])), " samples)")
  y_lab <- paste0("Single-Nucleus Fidelity") #\n(", comma(max(df[,grep("No.Nuclei", colnames(df))])), " nuclei)")
  
  if(prop_scaled){
    y_lab <- gsub("Fidelity", "Proportional Fidelity", y_lab)
  }
  
  datasets <- unique(df$SN_Dataset)
  
  r2_list <- lapply(1:length(datasets), function(i){
    
    r2 <- unlist(lapply(1:length(thresholds), function(k){
      idx <- which(
        is.element(df$SN_Threshold, thresholds[k])&is.element(df$SN_Dataset, datasets[i])
      )
      return(signif(
        summary(lm(df[idx, sn_col]~df[idx ,bulk_col]))$adj.r.squared, 2
      ))
    }))
    
    return(data.frame(r2, SN_Threshold=thresholds, SN_Dataset=datasets[i]))
    
  })
  r2_df <- do.call(rbind, r2_list)
  r2_df$Label <- paste("r^2 ==", r2_df$r2)

  df_temp <- merge(df, r2_df, by=c("SN_Threshold", "SN_Dataset"))
  df_temp <- merge(df_temp, datinfo, by.x="SN_Dataset", by.y="Dataset")
  df_temp$Plot_Label <- gsub("SSv4", "", df_temp$Plot_Label)
  df_temp$Plot_Label <- gsub("FC", "", df_temp$Plot_Label)
  df_temp$Plot_Label <- gsub("TC", "", df_temp$Plot_Label)
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/", cell_class, "_bulk_vs_SN_fidelity_scatterplot_UMIs_downsampled_", data_type, "_", expr_type, ".pdf")
  
  if(prop_scaled){
    file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
  }
  
  pdf(file=file_path, width=8, height=9)
  
  p <- ggplot(df_temp, aes_string(
    x=colnames(df)[bulk_col], y=colnames(df)[sn_col], 
    color=colnames(df)[mean_perc_col]
  )) + 
    geom_point(alpha=.75, size=.1) +
    theme_classic() + 
    theme(
      plot.title=element_text(hjust=.5, face="bold", size=16),
      plot.subtitle=element_text(
        hjust=.5, size=13, lineheight=1.1, margin=margin(t=4, b=8)
      ),
      legend.title=element_text(size=11, lineheight=1, margin=margin(b=3)),
      legend.position="bottom",
      legend.box="vertical",
      axis.title.x=element_text(size=13, margin=margin(t=15, b=6)),
      axis.title.y=element_text(size=13, margin=margin(r=12)),
      panel.border=element_blank(),
      panel.background=element_blank(),
      axis.ticks.x=element_line(color="black"),
      panel.grid.major=element_line(colour="grey", size=.15), # colour="grey", size=.15
      panel.grid.minor=element_line(colour="grey", size=.15),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) + 
    xlab(x_lab) + ylab(y_lab) + 
    scale_color_viridis(
      guide=guide_colourbar(
        title=leg_title, 
        barwidth=20, barheight=1.5,
        frame.colour="black",
        label.theme=element_text(size=10, angle=0, hjust=.8),
        title.position="top",
        title.hjust=.53,
      ),
      limits=c(0, 100), 
      labels=color_labs, 
      breaks=color_breaks
    ) +
    facet_grid(Plot_Label~SN_Threshold) +
    theme(
      strip.background=element_blank(),
      strip.text=element_text(face="bold", size=10),
      strip.placement="inside",
      panel.spacing.x=unit(.5, "lines")
    ) + geom_text(
      aes(label=Label, x=42, y=.875),
      color="black", size=2.75,
      parse=T, check_overlap=T
    )

  print(p)
  
  dev.off()
  
} ## plot_bulk_vs_SN_fidelity <- function(

plot_bulk_kme_vs_SN_fidelity <- function(datinfo, class_sn_fid, bulk_kme, cell_class=c("NEU", "ASC", "MIC", "OG", "OPC"), data_type, expr_type, prop_scaled=F){
  
  cell_class_name <- cell_class_full_name(cell_class)
  
  class_bulk_kme <- bulk_kme %>% 
    na.omit() %>%
    dplyr::select(
      SYMBOL, 
      No.Samples,
      colnames(bulk_kme)[grep("kME", colnames(bulk_kme))]
    )

  class_sn_fid$Threshold <- as.numeric(class_sn_fid$Threshold)
  class_sn_fid <- class_sn_fid[!is.element(class_sn_fid$Threshold, 1e6),]
  
  idx <- class_sn_fid$Threshold%%100>0
  pre_umis <- mean(class_sn_fid$Threshold[idx])
  class_sn_fid$Threshold[idx] <- 1e6
  class_sn_fid <- class_sn_fid[class_sn_fid$Threshold%in%c(500, 1e3, 1e4, 1e5, 1e6),]
  
  idx <- order(unique(class_sn_fid$Threshold))
  class_sn_fid$Threshold <- sapply(class_sn_fid$Threshold, comma)
  class_sn_fid$Threshold[class_sn_fid$Threshold=="1,000,000"] <- "Full coverage"
  class_sn_fid$Threshold[class_sn_fid$Threshold!="Full coverage"] <- paste(class_sn_fid$Threshold[class_sn_fid$Threshold!="Full coverage"], "UMIs")
  class_sn_fid$Threshold <- factor(
    class_sn_fid$Threshold, levels=unique(class_sn_fid$Threshold)[idx]
  )
  
  thresholds <- unique(class_sn_fid$Threshold)
  
  df_list <- lapply(1:length(thresholds), function(k){
    
    temp <- class_sn_fid[is.element(class_sn_fid$Threshold, thresholds[k]),]
    df <- merge(temp, class_bulk_kme, by=1)
    
    return(df)
    
  })
  df <- do.call(rbind, df_list)
  
  sn_col <- grep("Fidelity", colnames(df))
  bulk_col <- grep("kME", colnames(df))
  mean_perc_col <- grep("Mean_Expr_Percentile", colnames(df))
  #umi_col <- grep("Mean_UMIs", colnames(df))
  
  color_breaks <- seq(0, 100, by=20)
  color_breaks[which.min(color_breaks)] <- ceiling(min(df[,mean_perc_col]))
  color_breaks[which.max(color_breaks)] <- floor(max(df[,mean_perc_col]))
  color_labs <- round_any(color_breaks, 5)
  
  plot_title <- paste(cell_class_name, "Single-Nucleus Fidelity vs. Read Depth per Nucleus")
  plot_sub <- paste0(
    "Average median # UMIs per nucleus prior to downsampling: ", comma(pre_umis), "\n", comma(n_distinct(df$SYMBOL)), " protein coding genes"
  ) 
  
  leg_title <- paste("Mean Expression %tile") # cell_class_name, "Single-Nucleus 
  x_lab <- paste0("Bulk Turquoise Module kME\n(", comma(max(df[,grep("No.Samples", colnames(df))])), " samples)")
  y_lab <- paste0("Single-Nucleus Fidelity") #\n(", comma(max(df[,grep("No.Nuclei", colnames(df))])), " nuclei)")
  
  if(prop_scaled){
    y_lab <- gsub("Fidelity", "Proportional Fidelity", y_lab)
  }
  
  datasets <- unique(df$Dataset)
  
  r2_list <- lapply(1:length(datasets), function(i){
    
    r2 <- unlist(lapply(1:length(thresholds), function(k){
      idx <- which(
        is.element(df$Threshold, thresholds[k])&is.element(df$Dataset, datasets[i])
      )
      return(signif(
        summary(lm(df[idx, sn_col]~df[idx ,bulk_col]))$adj.r.squared, 2
      ))
    }))
    
    return(data.frame(r2, Threshold=thresholds, Dataset=datasets[i]))
    
  })
  r2_df <- do.call(rbind, r2_list)
  r2_df$Label <- paste("r^2 ==", r2_df$r2)
  
  df_temp <- merge(df, r2_df, by=c("Threshold", "Dataset"))
  df_temp <- merge(df_temp, datinfo, by="Dataset")
  df_temp$Plot_Label <- gsub("SSv4", "", df_temp$Plot_Label)
  df_temp$Plot_Label <- gsub("FC", "", df_temp$Plot_Label)
  df_temp$Plot_Label <- gsub("TC", "", df_temp$Plot_Label)
  
  file_path <- paste0("figures/", data_type, "/", expr_type, "/", cell_class, "_bulk_kME_vs_SN_fidelity_scatterplot_UMIs_downsampled_", data_type, "_", expr_type, ".pdf")
  
  if(prop_scaled){
    file_path <- gsub("fidelity", "fidelity_proportion_scaled", file_path)
  }
  
  pdf(file=file_path, width=10, height=11)
  
  p <- ggplot(df_temp, aes_string(
    x=colnames(df)[bulk_col], y=colnames(df)[sn_col], 
    color=colnames(df)[mean_perc_col]
  )) + 
    geom_point(alpha=.75, size=.1) +
    theme_classic() + 
    theme(
      plot.title=element_text(hjust=.5, face="bold", size=16),
      plot.subtitle=element_text(hjust=.5, size=13, lineheight=1.1, margin=margin(t=4, b=8)),
      legend.title=element_text(size=11, lineheight=1, margin=margin(b=3)),
      legend.position="bottom",
      legend.box="vertical",
      axis.title.x=element_text(face="bold", size=13, margin=margin(t=15, b=6)),
      axis.title.y=element_text(face="bold", size=13, margin=margin(r=12)),
      panel.border=element_blank(),
      panel.background=element_blank(),
      axis.ticks.x=element_line(color="black"),
      panel.grid.major=element_line(colour="grey", size=.15), # colour="grey", size=.15
      panel.grid.minor=element_line(colour="grey", size=.15),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) + 
    xlab(x_lab) + ylab(y_lab) + 
    scale_color_viridis(
      guide=guide_colourbar(
        title=leg_title, 
        barwidth=20, barheight=1.5,
        frame.colour="black",
        label.theme=element_text(size=10, angle=0, hjust=.8),
        title.position="top",
        title.hjust=.53,
      ),
      limits=c(0, 100), 
      labels=color_labs, 
      breaks=color_breaks
    ) +
    facet_grid(Plot_Label~Threshold) +
    theme(
      strip.background=element_blank(),
      strip.text=element_text(face="bold", size=10),
      strip.placement="inside",
      panel.spacing.x=unit(.5, "lines")
    ) + geom_text(
      aes(label=Label, x=.5, y=.75),
      color="black", size=2.75,
      parse=T, check_overlap=T
    )
  
  print(p)
  
  dev.off()
  
} ## plot_bulk_vs_SN_fidelity <- function(
