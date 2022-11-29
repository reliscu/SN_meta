library(ggplot2)
library(data.table)
library(dplyr)
library(viridis)
library(scales)
library(ggrepel)
library(RColorBrewer)

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_stats_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/continuous_plot_breaks_fxn.R")

palette <- readRDS("/home/rebecca/SCSN_meta_analysis/palette_15.RDS")

summary_plots <- function(datinfo){
  
  datinfo$Year <- as.factor(datinfo$Year)
  datinfo$FACS_Sorted <- factor(datinfo$FACS_Sorted, levels=c("Y", "N"))
  datinfo$Unbiased_Sampling <- factor(datinfo$Unbiased_Sampling, levels=c("Y", "N"))

  plot_sub <- paste(n_distinct(datinfo$Study), "studies,", n_distinct(datinfo$Dataset), "regional datasets")

  pdf("figures/platform_barplot.pdf", width=15, height=11)
  
  df <- datinfo %>% 
    dplyr::group_by(Platform) %>% 
    dplyr::summarise(n=n()) %>%
    dplyr::arrange(n) %>%
    dplyr::mutate(
      Platform=factor(Platform, levels=Platform)
    )
  
  x_var <- "Platform"; x_lab <- "Platforms"
  
  ggplot(df, aes(x=!!as.name(x_var), y=n, fill=!!as.name(x_var))) +
    geom_bar(stat="identity", width=.75, show.legend=F) +
    theme_minimal() +
    theme(
      plot.title=element_text(face="bold", hjust=.5, size=20),
      plot.subtitle=element_blank(), # element_text(hjust=.5, size=14, lineheight=1, margin=margin(t=2, b=10))
      axis.title.x=element_blank(),
      axis.text.x=element_text(color="black", vjust=.5, size=8),
      axis.title.y=element_text(color="black", size=12, margin=margin(r=30)),
      axis.text.y=element_text(color="black", size=12),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    labs(title=x_lab, subtitle=plot_sub) +
    ylab("# Datasets") +
    scale_fill_brewer(palette="Dark2", direction=-1) +
    scale_y_continuous(limits=c(0, max(df$n)), breaks=seq(0, max(df$n), by=2)) 
  
  dev.off()
  
  pdf("figures/region_barplot.pdf", width=12, height=9)
  
  df <- datinfo %>% 
    dplyr::group_by(Region_Code) %>% 
    dplyr::summarise(n=n()) %>%
    dplyr::arrange(n) %>%
    dplyr::mutate(
      Region_Code=factor(Region_Code, levels=Region_Code)
    )
  
  x_var <- "Region_Code"; x_lab <- "Regions"
  
  ggplot(df, aes(x=!!as.name(x_var), y=n, fill=!!as.name(x_var))) +
    geom_bar(stat="identity", width=.75, show.legend=F) +
    theme_minimal() +
    theme(
      plot.title=element_text(face="bold", hjust=.5, size=20),
      plot.subtitle=element_blank(), # element_text(hjust=.5, size=14, lineheight=1, margin=margin(t=2, b=10))
      axis.title.x=element_blank(),
      axis.text.x=element_text(color="black", size=14),
      axis.title.y=element_text(color="black", size=13, margin=margin(r=30)),
      axis.text.y=element_text(color="black", size=12),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    labs(title=x_lab, subtitle=plot_sub) +
    ylab("# Datasets") +
    scale_fill_brewer(palette="Set1", direction=-1) +
    scale_y_continuous(limits=c(0, max(df$n)), breaks=seq(0, max(df$n), by=5))
  
  dev.off()
  
  pdf("figures/FACS_pie_chart.pdf", width=8, height=8)
  
  df <- datinfo %>% 
    dplyr::group_by(FACS_Sorted) %>% 
    dplyr::summarise(n=n()) %>%
    dplyr::arrange(n) %>%
    dplyr::mutate(
      FACS_Sorted=factor(FACS_Sorted, levels=c("Y", "N"))
    )
  
  df$Labs <- NA
  df$Labs[df$FACS_Sorted=="Y"] <- paste0("Y\n(", df$n[df$FACS_Sorted=="Y"], ")")
  df$Labs[df$FACS_Sorted=="N"] <- paste0("N\n(", df$n[df$FACS_Sorted=="N"], ")")
  
  x_var <- "FACS_Sorted"; x_lab <- "FACS"
  
  ggplot(df, aes(x="", y=n, fill=!!as.name(x_var))) +
    geom_bar(stat="identity", show.legend=F) +
    theme_minimal() +
    theme(
      plot.title=element_text(size=27, hjust=.5, face="bold", margin=margin(b=-20)),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text=element_blank(),
      panel.border=element_blank(),
      panel.grid=element_blank(),
      axis.ticks=element_blank(),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    coord_polar("y", start=1) +
    labs(title=x_lab) +
    scale_fill_brewer(palette="Set1", direction=-1, labels=labs) +
    geom_text(aes(label=Labs), color="white", size=10, position=position_stack(vjust=.5))

  dev.off()
  
  pdf("figures/sampling_pie_chart.pdf", width=8, height=8)
  
  df <- datinfo %>% 
    dplyr::group_by(Unbiased_Sampling) %>% 
    dplyr::summarise(n=n()) %>%
    dplyr::arrange(n) %>%
    dplyr::mutate(
      Unbiased_Sampling=factor(Unbiased_Sampling, levels=c("Y", "N"))
    )
  
  df$Labs <- NA
  df$Labs[df$Unbiased_Sampling=="Y"] <- paste0("Y\n(", df$n[df$Unbiased_Sampling=="Y"], ")")
  df$Labs[df$Unbiased_Sampling=="N"] <- paste0("N\n(", df$n[df$Unbiased_Sampling=="N"], ")")

  x_var <- "Unbiased_Sampling"; x_lab <- "Unbiased Sampling"
  
  print(
    ggplot(df, aes(x="", y=n, fill=!!as.name(x_var))) +
      geom_bar(stat="identity", show.legend=F) +
      theme_minimal() +
      theme(
        plot.title=element_text(size=27, hjust=.5, face="bold", margin=margin(b=-20)),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        plot.margin=unit(c(2, 2, 2, 2), "cm")
      ) +
      coord_polar("y", start=1) +
      labs(title=x_lab) +
      scale_fill_brewer(palette="Set1", direction=-1) +
      geom_text(aes(label=Labs), color="white", size=10, position=position_stack(vjust=.5))
  )
  
  dev.off()

  
} ## summary_plots <- function(

trend_plots <- function(datinfo){
  
  datinfo$Year <- as.factor(datinfo$Year)
  datinfo$FACS_Sorted <- factor(datinfo$FACS_Sorted, levels=c("Y", "N"))
  datinfo$Unbiased_Sampling <- factor(datinfo$Unbiased_Sampling, levels=c("Y", "N"))

  pdf(file=paste0("figures/summary_plots_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies.pdf"), width=12, height=10)
  
  x_var <- "Plot_Label"; x_lab <- "Dataset"
  x_var_type <- "discrete"
  y_var <- "Author_Counts_QC_No.Genes_Mean_UMIs_Greaterthan_1_Protein_Coding"; y_lab <- "Number of Genes with Mean UMIs >=1"
  y_var_type <- "continuous"
  plot_template(df=datinfo, x_var, x_lab, x_var_type, y_var, y_lab, y_var_type, color_var=NULL, color_lab=NULL, color_var_type=NULL)
  
  x_var <- "Plot_Label"; x_lab <- "Dataset"
  x_var_type <- "discrete"
  y_var <- "Author_Median_Unique_Genes_QC_Protein_Coding"; y_lab <- "Median # Unique Genes per Nucleus"
  y_var_type <- "continuous"
  plot_template(df=datinfo, x_var, x_lab, x_var_type, y_var, y_lab, y_var_type, color_var=NULL, color_lab=NULL, color_var_type=NULL)
  
  x_var <- "Author_No.Nuclei_QC"; x_lab <- "# Nuclei"
  y_var <- "Author_No.Cell_Types"; y_lab <- "# Cell Types"
  color_var <- "Platform"; color_lab <- "Platform"
  color_var_type <- "discrete"
  text_template(df=datinfo, x_var, x_lab, y_var, y_lab, color_var, color_lab, color_var_type)
  
  x_var <- "Author_No.Nuclei_QC"; x_lab <- "# Nuclei"
  y_var <- "Author_Median_Unique_Genes_QC_Protein_Coding"; y_lab <- "Median # Unique Genes per Nucleus"
  color_var <- "Year"; color_lab <- "Year"
  color_var_type <- "discrete"
  text_template(df=datinfo, x_var, x_lab, y_var, y_lab, color_var, color_lab, color_var_type)
  
  x_var <- "Author_No.Nuclei_QC"; x_lab <- "# Nuclei"
  y_var <- "Author_Median_Unique_Genes_QC_Protein_Coding"; y_lab <- "Median # Unique Genes per Nucleus"
  color_var <- "Platform"; color_lab <- "Platform"
  color_var_type <- "discrete"
  text_template(df=datinfo, x_var, x_lab, y_var, y_lab, color_var, color_lab, color_var_type)

  x_var <- "Author_No.Nuclei_QC"; x_lab <- "# Nuclei"
  x_var_type <- "coninuous"
  y_var <- "Author_Median_Unique_Genes_QC_Protein_Coding"; y_lab <- "Median # Unique Genes per Nucleus"
  y_var_type <- "continuous"
  color_var <- "Platform"; color_lab <- "Platform"
  color_var_type <- "discrete"
  plot_template(df=datinfo, x_var, x_lab, x_var_type, y_var, y_lab, y_var_type, color_var, color_lab, color_var_type)
  
  x_var <- "Author_No.Nuclei_QC"; x_lab <- "# Nuclei"
  x_var_type <- "coninuous"
  y_var <- "Author_No.Cell_Types"; y_lab <- "# Cell Types"
  y_var_type <- "continuous"
  color_var <- "Author_Median_Unique_Genes_QC_Protein_Coding"; color_lab <- "Median # Unique Genes per Nucleus"
  color_var_type <- "continuous"
  plot_template(df=datinfo, x_var, x_lab, x_var_type, y_var, y_lab, y_var_type, color_var, color_lab, color_var_type)
  
  x_var <- "Author_No.Nuclei_QC"; x_lab <- "# Nuclei"
  x_var_type <- "coninuous"
  y_var <- "Author_Median_Unique_Genes_QC_Protein_Coding"; y_lab <- "Median # Unique Genes per Nucleus"
  y_var_type <- "continuous"
  color_var <- "Study"; color_lab <- "Study"
  color_var_type <- "discrete"
  
  x_var <- "Author_Median_Unique_Genes_QC_Protein_Coding"; x_lab <- "Median # Unique Genes per Nucleus"
  x_var_type <- "coninuous"
  y_var <- "Author_No.Cell_Types"; y_lab <- "# Cell Types Reported"
  y_var_type <- "continuous"
  color_var <- "Platform"; color_lab <- "Platform"
  color_var_type <- "discrete"
  
  dev.off()

} ## trend_plots <- function(

plot_template <- function(df, x_var, x_lab, x_var_type=c("continuous", "discrete"), y_var, y_lab, y_var_type=c("continuous", "discrete"), color_var=NULL, color_lab=NULL, color_var_type=NULL){
  
  plot_sub <- paste(n_distinct(df$Dataset), "regional datasets from", n_distinct(df$Study), "studies")
  
  if(grepl("Genes", y_var)){
  
    plot_sub <- paste0(plot_sub, "\nRestricted to protein coding genes")
    
  }
  
  if(x_var_type=="discrete"){
    
    n_elements <- n_distinct(df[,x_var])
    
  } else {
    
    n_elements <- n_distinct(df[,color_var])
    
  }
  
  if(n_elements>15){
    
    color_vals <- brewer_fxn(n_elements)
    
  } else if(n_elements==2){
    
    color_vals <- c("#3492EA", "#D92F4B")
    
  } else {
    
    color_vals <- palette[1:n_elements]
    
  }
  
  if(x_var_type=="discrete"&y_var_type=="continuous"){

    plot_title <- paste(y_lab)
    
    if(grepl("Genes", y_var)){
      y_lab <- "# Genes"
    }
    
    temp <- df %>%
      dplyr::group_by(!!as.name(x_var)) %>%
      dplyr::summarise(n=mean(!!as.name(y_var))) %>%
      dplyr::arrange(n) %>%
      as.data.frame()
    
    df[,x_var] <- factor(df[,x_var], levels=unique(temp[,x_var]))

    if(x_var=="Platform"){
      
      p <- ggplot(df, aes_string(x=x_var, y=y_var, fill=x_var)) +
        geom_boxplot(width=.5, outlier.shape=NA) +
        geom_jitter(size=.3, width=.1, color="black", show.legend=F) +
        theme_minimal() +
        theme(
          plot.title=element_text(hjust=.5, size=14, face="bold", margin=margin(b=10)),
          plot.subtitle=element_text(hjust=.5, size=12, lineheight=1.1, margin=margin(b=15)),
          legend.title=element_text(face="bold", size=14, margin=margin(b=10)),
          legend.text=element_text(size=11),
          legend.position="bottom",
          legend.direction="horizontal",
          legend.box="vertical",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_text(margin=margin(r=15)),
          axis.text.y=element_text(color="black"),
          plot.margin=unit(c(2, 2, 2, 2), "cm")
        ) +
        labs(title=plot_title, subtitle=plot_sub) +
        ylab(y_lab) +
        scale_y_continuous(labels=comma) +
        scale_color_manual(values=color_vals) + scale_fill_manual(values=color_vals)
      
      if(n_elements>10){
        
        p <- p + theme(legend.title=element_blank(), axis.text.x=element_text(color="black", size=11, hjust=1, angle=45, margin=margin(t=-28))) + guides(fill="none", color="none")
        
      } else {

        p <- p + guides(fill=guide_legend(byrow=F, keywidth=1, keyheight=1, title=x_lab, title.position="top", title.hjust=.5, override.aes=list(size=.5)))
        
      }
      
    } else if(grepl("Genes", y_var)){ ## if(x_var=="Platform"){
    
      y_breaks <- seq(0, 15e3, by=2e3)
      if(grepl("Median", y_var)){
        y_breaks <- seq(0, round_any(max(df[,y_var]), 1e3, f=ceiling), by=500)
      }
      
      p <- ggplot(df, aes_string(x=x_var, y=y_var, fill=x_var)) +
        geom_bar(stat="identity", width=.75, color="black") +
        theme_minimal() +
        theme(
          plot.title=element_text(hjust=.5, size=14, face="bold", margin=margin(b=10)),
          plot.subtitle=element_text(hjust=.5, size=12, lineheight=1.1, margin=margin(b=15)),
          legend.title=element_text(face="bold", size=14, margin=margin(b=10)),
          legend.text=element_text(size=11),
          legend.position="bottom",
          legend.direction="horizontal",
          legend.box="vertical",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_text(margin=margin(r=10)),
          axis.text.y=element_text(color="black"),
          plot.margin=unit(c(2, 2, 2, 2), "cm")
        ) +
        labs(title=plot_title, subtitle=plot_sub) +
        ylab(y_lab) +
        scale_y_continuous(breaks=y_breaks, labels=comma) +
        scale_color_manual(values=color_vals) + scale_fill_manual(values=color_vals)
      
      if(n_elements>10){
        
        p <- p + theme(legend.title=element_blank(), axis.text.x=element_text(color="black", size=9, hjust=1, angle=45, margin=margin(t=-34))) + guides(fill="none", color="none")
        
      } else {
        
        p <- p + guides(color=guide_legend(byrow=F, keywidth=.85, keyheight=.85, title=x_lab, title.position="top", title.hjust=.5, override.aes=list(size=.5)))
        
      }
      
    } ## else if(grepl("Genes", y_var)){
   
  } else { ## else: x_var_type=="continuous"&y_var_type="continuous"
    
    plot_title <- paste(x_lab, "vs.", y_lab)
    if(nchar(y_lab)>30){
      plot_title <- paste(x_lab, "vs.\n", y_lab)
    }
    
    if(color_var=="Year"){
      
      color_vals <- colorRampPalette(brewer.pal(9, "YlOrRd"))(n_elements)
      
    }  
      
    if(color_var_type=="discrete"){
      
      n_row <- 2
      if(color_var=="Year"){
        n_row <- 1
      }
      
      p <- 
        ggplot(df, aes_string(x=x_var, y=y_var, fill=color_var)) +
        geom_point(size=3, alpha=1, shape=21) +
        theme_classic() + 
        theme(
          plot.title=element_text(hjust=.5, face="bold", size=14),
          plot.subtitle=element_text(hjust=.5, size=12, lineheight=1, margin=margin(t=4, b=8)),
          legend.title=element_text(size=12, face="bold"),
          legend.position="bottom",
          legend.box="vertical",
          legend.text=element_text(size=12),
          axis.title.x=element_text(size=12, color="black", margin=margin(t=12)), 
          axis.text.x=element_text(angle=0, margin=margin(t=4)),
          axis.title.y=element_text(size=12, color="black", margin=margin(r=12)),
          plot.margin = unit(c(2, 2, 2, 2), "cm")
        ) +
        labs(title=plot_title, subtitle=plot_sub) + 
        xlab(x_lab) + ylab(y_lab) +
        guides(fill=guide_legend(
          nrow=n_row,
          title.position="left",
          title.hjust=.5,
          title=color_lab,
          override.aes=list(size=5)
        )) +
        scale_x_continuous(label=comma) + scale_y_continuous(label=comma) +
        scale_fill_manual(values=color_vals) 
      
    } else { ## if(color_var_type=="discrete"){
      
      ## color_var_type=="continuous
      
      if(grepl("Genes", color_var)){
        color_breaks <- n_unique_genes_breaks(df[,color_var])
      } else if(grepl("UMIs", color_var)) {
        color_breaks <- n_UMIs_breaks(df[,color_var])
      } else {
        color_breaks <- n_nuclei_breaks(df[,color_var])
      }
      
      color_labs <- comma(color_breaks)
      
      p <- 
        ggplot(df, aes_string(x=x_var, y=y_var, fill=color_var)) + 
        geom_point(size=3, alpha=1, shape=21) +
        theme_classic() + 
        theme(
          plot.title=element_text(hjust=.5, face="bold", size=14),
          plot.subtitle=element_text(hjust=.5, size=12, lineheight=1, margin=margin(t=4, b=8)),
          legend.title=element_text(size=12, face="bold", margin=margin(b=5)),
          legend.position="bottom",
          legend.box="vertical",
          legend.text=element_text(size=12),
          axis.title.x=element_text(size=12, color="black", margin=margin(t=12, b=10)), 
          axis.text.x=element_text(angle=0, margin=margin(t=4)),
          axis.title.y=element_text(size=12, color="black", margin=margin(r=12)),
          plot.margin = unit(c(2, 2, 2, 2), "cm")
        ) +
        labs(title=plot_title, subtitle=plot_sub) + 
        xlab(x_lab) + ylab(y_lab) +
        scale_x_continuous(label=comma) + scale_y_continuous(label=comma) +
        scale_fill_viridis(
          guide=guide_colourbar(
            title=color_lab, 
            barwidth=20, barheight=1,
            frame.colour="black",
            ticks.colour="black",
            label.theme=element_text(
              size=10
            ),
            title.position="top",
            title.hjust=.53,
          ),
          labels=color_labs,
          breaks=color_breaks
        )
      
    } ## if(color_var_type=="discrete"){} else {
    
  } ## if(x_var_type=="discrete"){} else {
  
  print(p)
  
} ## plot_template <- function(

text_template <- function(df, x_var, x_lab, y_var, y_lab, color_var, color_lab, color_var_type=c("continuous", "discrete")){
  
  plot_title <- paste(x_lab, "vs.", y_lab)
  
  temp1 <- df %>%
    dplyr::group_by(Plot_Label) %>%
    dplyr::summarise(Y=mean(!!as.name(y_var)))
  
  temp2 <- df %>%
    dplyr::group_by(Plot_Label) %>%
    dplyr::summarise(X=mean(!!as.name(x_var))) 
  
  temp <- merge(temp1, temp2, by="Plot_Label")
  temp <- merge(temp, df, by="Plot_Label")
  
  plot_sub <- paste(n_distinct(temp$Dataset), "regional datasets from", n_distinct(temp$Study), "studies")
  
  if(grepl("Genes", y_var)){
    
    plot_sub <- paste0(plot_sub, "\nRestricted to protein coding genes")
    
  }
  
  n_elements <- n_distinct(df[,color_var])
  
  if(color_var=="Year"){
    
    color_vals <- colorRampPalette(brewer.pal(9, "YlOrRd"))(n_elements)
    
  } else if(n_elements==2){
    
    color_vals <- c("#3492EA", "#D92F4B")
    
  } else {
    
    color_vals <- brewer_fxn(n_elements)
    
  } 
  
  n_row <- 2
  if(color_var=="Year"){
    n_row <- 1
  }
  
  if(color_var_type=="discrete"){

    print(
      ggplot(temp, aes_string(x="X", y="Y", color=color_var)) +
        geom_point(size=1, alpha=1, show.legend=T) +
        geom_text_repel(
          aes(label=Plot_Label), size=5, fontface="bold", show.legend=F, max.overlaps=1000, seed=666, box.padding=.3
        ) +
        theme_classic() + 
        theme(
          plot.title=element_text(hjust=.5, face="bold", size=16),
          plot.subtitle=element_text(hjust=.5, size=13, lineheight=1, margin=margin(t=4, b=8)),
          legend.title=element_blank(),
          legend.position="bottom",
          legend.box="vertical",
          legend.text=element_text(size=12),
          axis.title.x=element_text(size=15, color="black", margin=margin(t=12)), 
          axis.text.x=element_text(size=11, color="black", margin=margin(t=4)),
          axis.title.y=element_text(size=15, color="black", margin=margin(r=18)),
          axis.text.y=element_text(size=11, color="black"),
          plot.margin = unit(c(2, 2, 2, 2), "cm")
        ) +
        labs(title=plot_title, subtitle=plot_sub) + 
        xlab(x_lab) + ylab(y_lab) +
        scale_x_continuous(label=comma) + scale_y_continuous(label=comma) +
        scale_color_manual(values=color_vals) +
        guides(color=guide_legend(
          nrow=n_row, title=color_var, title.position="left",
          override.aes=list(size=8)
        ))
    )

    
  } else {
    
    if(grepl("Genes", color_var)){
      color_breaks <- n_unique_genes_breaks(df[,color_var])
    } else if(grepl("UMIs", color_var)) {
      color_breaks <- n_UMIs_breaks(df[,color_var])
    } else {
      color_breaks <- n_nuclei_breaks(df[,color_var])
    }
    
    color_labs <- comma(color_breaks)
    
    print(
      ggplot(temp, aes_string(x="X", y="Y", color=color_var)) +
        geom_point(size=1, alpha=1, show.legend=T) +
        geom_text_repel(
          aes(label=Plot_Label), size=5, fontface="bold", show.legend=F, max.overlaps=1000, seed=666, box.padding=.3
        ) +
        theme_classic() + 
        theme(
          plot.title=element_text(hjust=.5, face="bold", size=16),
          plot.subtitle=element_text(hjust=.5, size=13, lineheight=1, margin=margin(t=4, b=8)),
          legend.title=element_text(face="bold", size=14, margin=margin(b=5)),
          legend.position="bottom",
          legend.box="vertical",
          legend.text=element_text(size=12),
          axis.title.x=element_text(size=15, color="black", margin=margin(t=12, b=12)), 
          axis.text.x=element_text(size=11, color="black", margin=margin(t=4)),
          axis.title.y=element_text(size=15, color="black", margin=margin(r=18)),
          axis.text.y=element_text(size=11, color="black"),
          plot.margin = unit(c(2, 2, 2, 2), "cm")
        ) +
        labs(title=plot_title, subtitle=plot_sub) + 
        xlab(x_lab) + ylab(y_lab) +
        scale_x_continuous(label=comma) + scale_y_continuous(label=comma) +
        scale_color_viridis(
          guide=guide_colourbar(
            title=color_lab, 
            barwidth=20, barheight=1,
            frame.colour="black",
            ticks.colour="black",
            label.theme=element_text(
              size=10, angle=45, vjust=.5
            ),
            title.position="top",
            title.hjust=.53,
          ),
          labels=color_labs,
          breaks=color_breaks
        )
    )
    
  }
  
}

# cell_class_prop <- function(
#   datinfo, 
#   cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"), 
#   na_cts=c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class")
# ){
#   
#   datinfo <- datinfo %>% na_if("") %>% filter(!is.na(Author_Barcode_Annotations))
#   
#   cell_prop_df <- c()
#   
#   for(i in 1:nrow(datinfo)){
#     
#     cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F) %>% na_if("")
#     
#     if(is.element("MO_Cell_Class", colnames(cellinfo))){
# 
#       ## Don't include NA/ambiguous cell types in total CT counts:
#       
#       cellinfo <- cellinfo[!is.na(cellinfo$Cell_Type),]
#       cellinfo <- cellinfo[!grepl(paste(na_cts, collapse="|"), cellinfo$Cell_Type),]
#       cellinfo$MO_Cell_Class[!is.element(cellinfo$MO_Cell_Class, cell_classes)] <- "Other" 
#         
#       cell_prop <- data.frame(
#         table(cellinfo$MO_Cell_Class)/nrow(cellinfo)
#       )
#       cell_prop$Dataset <- datinfo$Dataset[i]
#        
#       if(is.null(cell_prop_df)){
#         cell_prop_df <- cell_prop
#       } else {
#         cell_prop_df <- rbind(cell_prop_df, cell_prop)
#       }
#       
#     } ## if(is.element("MO_Cell_Class", colnames(cellinfo))){
#     
#   } ## for(i in 1:nrow(datinfo)){
#   
#   colnames(cell_prop_df) <- c("Cell_Class", "Proportion", "Dataset")
#   
#   df <- merge(cell_prop_df, datinfo, by="Dataset")
#   
#   ## Arrange datasets in order of EXC proportion:
#   
#   ds_order <- df %>% 
#     dplyr::filter(Cell_Class=="EXC") %>%
#     dplyr::arrange(Proportion)
#   
#   df$Plot_Label <- factor(df$Plot_Label, levels=c("Grubman TC EC", ds_order$Plot_Label))
#   df$Cell_Class <- factor(df$Cell_Class, levels=c(sort(cell_classes), "Other"))
#   df$Unbiased_Sampling[df$Unbiased_Sampling=="Y"] <- "Unbiased Sampling"
#   df$Unbiased_Sampling[df$Unbiased_Sampling=="N"] <- "Biased Sampling"
#   df$Study <- paste(df$First_Author, df$Year)
#   df$Study[is.element(df$PubMedID, "35165441")] <- paste(df$Study[is.element(df$PubMedID, "35165441")], " ")
#   
#   pdf(file=paste0("figures/cell_proportion_barplot_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(df$Study), "_studies.pdf"), width=9, height=7)
#   
#   plot_title <- paste("Cell Class Proportion")
#   plot_sub <- paste(n_distinct(df$Dataset), "datasets from", n_distinct(df$Study), "studies")
#   
#   color_vec <- c("#254EC4", "#33B2FF", "#E53545", "#FF9F33", "#25C48D", "#F78DD7", "#F5DB1D", "#C425C4", "#FF33C4", "#804AD8", "#989898")
#   
#   p <- ggplot(df, aes(x=Plot_Label, y=Proportion, fill=Cell_Class)) +
#     geom_bar(stat="identity") +
#     theme_classic() +
#     theme(
#       plot.title=element_text(hjust=.5, face="bold", size=15),
#       plot.subtitle=element_text(
#         color="#616161", hjust=.5, size=12, 
#         lineheight=1, margin=margin(t=4, b=15)
#       ),
#       legend.position="right",
#       legend.box="vertical",
#       legend.title=element_blank(),
#       legend.text=element_text(size=10),
#       axis.title.x=element_blank(), 
#       axis.text.x=element_text(angle=90, size=9, hjust=1, vjust=.5, margin=margin(t=-5, r=0)), # element_text(angle=45, size=9, hjust=1, margin=margin(t=5, r=11))
#       axis.line.x.bottom=element_blank(),
#       axis.ticks.x=element_blank(),
#       axis.title.y=element_text(size=12, color="black", margin=margin(r=8)),
#       axis.text.y=element_blank(),
#       axis.ticks.y=element_blank(),
#       axis.line.y.left=element_blank(),
#       plot.margin=unit(c(2, 2, 2, 2), "cm")
#     ) +
#     labs(title=plot_title, subtitle=plot_sub) +
#     guides(fill=guide_legend(
#       title.position="top"
#     )) +
#     scale_fill_manual(values=color_vec)
#   
#   print(p)
#   
#   print(
#     p + facet_grid(.~Unbiased_Sampling, scales="free_x") +
#       theme(strip.text.x=element_text(size=12, face="bold"), 
#             strip.background=element_rect(color="white"))
#   )
#   
#   dev.off()
#   
# } ## cell_class_proportion <- function(