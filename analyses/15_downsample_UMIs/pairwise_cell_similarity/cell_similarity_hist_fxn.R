library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(future.apply)
library(flexiblas)

options(future.globals.maxSize=Inf)
plan(multicore, workers=3)

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/code/misc/upper_first.R")

plot_cell_similarity <- function(datinfo, cell_sim, sim_type){
  
  cell_sim_list <- future_lapply(1:length(cell_sim_paths), FUN=function(i){
    
    cell_sim <- readRDS(cell_sim_paths[i])
    cell_sim <- reshape2::melt(cell_sim)[,c(2,3)]
    colnames(cell_sim) <- c("Threshold", "Sim")
    
    cell_sim$Dataset <- strsplit(
      strsplit(cell_sim_paths[i], "/")[[1]][4], "_UMIs"
    )[[1]][1]
    cell_sim$Threshold <- as.numeric(
      gsub("Threshold_", "", cell_sim$Threshold)
    )
    cell_sim <- cell_sim[!is.element(cell_sim$Threshold, 1e6),]
    
    idx <- which(
      cell_sim$Threshold%%100>0
    )
    
    pre_umis <- mean(cell_sim$Threshold[idx])
    print(pre_umis)
    cell_sim$Threshold[idx] <- 1e6
    
    cell_sim <- cell_sim %>% 
      dplyr::group_by(Threshold) %>%
      dplyr::slice_sample(n=1e4) %>%
      as.data.frame()
    
    return(cell_sim)
    
  }) # future.seed=T
  
  df <- do.call(rbind, cell_sim_list)
  
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

  mean_umis <- comma(floor(mean(c(1010673, 883430, 974742, 885429, 1651414))))
  mean_nuclei <- comma(floor(mean(df$Author_No.Nuclei_QC)))
  
  plot_title <- "Pairwise Cell Similarity vs. Read Depth"
  plot_sub <- paste("Average # nuclei:", mean_nuclei, "\nAverage median # UMIs per nucleus prior to downsampling:", mean_umis)
  
  if(sim_type=="pearson"){
    
    y_lab <- paste(upper_first(sim_type), "Correlation")
    
  } else if(sim_type=="prop"){
    
    y_lab <- "Proportionality"
    
  } else {
    
    y_lab <- "Cosine Similarity"
    
  }
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/cell_pairwise_", toupper(sim_type), "_similarity_UMIs_downsampled_", data_type, "_", expr_type, "_", n_distinct(df$Dataset), "_datasets.pdf"), width=8.5, height=8)
  
  p <- ggplot(df, aes(x=Threshold, y=Sim, fill=Plot_Label)) + 
    geom_violin(trim=T, size=.5, position=position_dodge(0.9), show.legend=F) +
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
      axis.ticks.y=element_line(size=.4),
      panel.grid.major.y=element_line(colour="darkgrey", size=.2),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    xlab("# UMIs per nucleus") +
    scale_y_continuous(name=y_lab, breaks=seq(0, 1, by=.2)) +
    scale_fill_manual(values=brewer.pal(5, "Set2")) +
    facet_grid(Plot_Label~.) +
    theme(strip.background=element_blank())
  
  print(p)
  
  dev.off()
  
}

