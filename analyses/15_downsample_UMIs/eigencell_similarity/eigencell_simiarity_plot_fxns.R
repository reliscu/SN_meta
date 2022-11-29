library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=3)

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/code/misc/upper_first.R")

plot_eigencell_similarity_distro <- function(datinfo, df, sim_type){
 
  df <- df[df$Threshold!=1e6,]
  idx <- which(df$Threshold%%100>0)
  pre_umis <- mean(df$Threshold[idx])
  df$Threshold[idx] <- 1e6
  idx <- order(
    unique(df$Threshold)
  )
  df$Threshold <- sapply(df$Threshold, comma)
  df$Threshold[df$Threshold=="1,000,000"] <- "Full coverage"
  df$Threshold <- factor(
    df$Threshold, levels=unique(df$Threshold)[idx]
  )
  
  df <- merge(df, datinfo, by="Dataset")
  
  df$Plot_Label <- gsub("SSv4", "", df$Plot_Label)
  df$Plot_Label <- gsub("SSv4", "", df$Plot_Label)
  df$Plot_Label <- gsub("FC", "", df$Plot_Label)
  df$Plot_Label <- gsub("TC", "", df$Plot_Label)
  
  mean_nuclei <- comma(floor(mean(df$Author_No.Nuclei_QC)))
  
  plot_title <- "Cell vs. Eigencell Similarity"
  plot_sub <- paste("Average # nuclei:", mean_nuclei, "\nAverage median # UMIs per nucleus prior to downsampling:", comma(floor(pre_umis)))
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/cell_pairwise_", toupper(sim_type), "_similarity_UMIs_downsampled_", data_type, "_", expr_type, "_", n_distinct(df$Dataset), "_datasets.pdf"), width=8.5, height=8)
  
  p <- ggplot(df, aes(x=Threshold, y=True_Sim, fill=Plot_Label)) + 
    geom_violin(size=.5, position=position_dodge(0.9), show.legend=F) +
    geom_boxplot(
      notch=T, width=.08, size=.3, fill="white", color="black", 
      position=position_dodge(0.9), outlier.shape=NA, show.legend=F
    ) +
    theme_classic() + 
    theme(
      plot.title=element_text(hjust=.5, face="bold", size=15, margin=margin(b=10)),
      plot.subtitle=element_text(hjust=.5, size=11, lineheight=1.1, margin=margin(t=4, b=8)),
      legend.title=element_blank(),
      axis.title.x=element_text(size=13, color="black",  margin=margin(t=10, b=10)),
      axis.text.x=element_text(size=10, angle=45, hjust=1, color="black"), 
      axis.ticks.x=element_blank(),
      axis.title.y=element_text(size=12, color="black", margin=margin(r=15)),
      axis.text.y=element_text(color="black"),
      axis.ticks.y=element_line(linewidth=.4),
      panel.grid.major.y=element_line(colour="darkgrey", size=.2),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    xlab("# UMIs per nucleus") +
    scale_y_continuous(name=paste(upper_first(sim_type), "Similarity"), breaks=seq(0, 1, by=.2)) +
    scale_fill_manual(values=brewer.pal(5, "Set2")) +
    facet_grid(Plot_Label~.) +
    theme(strip.background=element_blank())
  
  print(p)
  
  dev.off()
  
}
