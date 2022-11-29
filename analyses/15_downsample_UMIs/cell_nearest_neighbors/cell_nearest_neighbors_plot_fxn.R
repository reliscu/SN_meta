library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/code/misc/upper_first.R")

plot_percent_matching_neighbors <- function(datinfo, nn_paths, sim_type){
  
  df <- lapply(1:length(nn_paths), function(i){
    
    nn_df <- read.csv(nn_paths[i])
    #nn_df <- nn_df[!is.element(nn_df$Threshold, 1e6),]
    
    max_depth <- max(unique(nn_df$Threshold))
    
    ## Identify nearest neighbor at full coverage:
    
    idx <- which(
      is.element(nn_df$Threshold, max_depth)
    )
    nn_df_max <- nn_df[,]
    nn_df_max <- nn_df_max %>% dplyr::select(
      Cell_ID, Neighbor_ID
    )
    
    nn_df <- nn_df[-idx,]
    
    thresholds <- unique(nn_df$Threshold)
    
    nn_percent <- lapply(1:length(thresholds), function(j){
      
      nn_df1 <- nn_df[is.element(nn_df$Threshold, thresholds[j]),]
      
      return(
        sum(nn_df_max$Neighbor_ID==nn_df1$Neighbor_ID)/nrow(nn_df_max)*100
      )
      
    }) ## nn_fraction <- lapply(
    
    percent_df <- data.frame(
      Percent=c(unlist(nn_percent), 100),
      Threshold=c(thresholds, max_depth),
      Dataset=unique(nn_df$Dataset)
    )
    
    return(percent_df)
    
  })
  df <- do.call(rbind, df)
  
  idx <- df$Threshold%%100>0
  pre_umis <- mean(df$Threshold[idx])
  #df$Threshold[idx] <- 1e6
  df <- df[!idx,]

  idx <- order(unique(df$Threshold))
  df$Threshold <- sapply(df$Threshold, comma)
  #df$Threshold[df$Threshold=="1,000,000"] <- "Full coverage"
  df$Threshold <- factor(df$Threshold, levels=unique(df$Threshold)[idx])
  
  df <- merge(df, datinfo, by="Dataset")
  
  df$Plot_Label <- gsub(" SSv4", "", df$Plot_Label)
  df$Plot_Label <- gsub("FC ", "", df$Plot_Label)
  df$Plot_Label <- gsub("TC ", "", df$Plot_Label)
  
  df <- df[is.element(df$Plot_Label, c(
    "ABI ACC",
    "Bakken M1",
    "Hodge MTG"
  )),]
  
  plot_title <- "% Cells with the Same Nearest Neighbor After Downsampling"
  plot_sub <- paste("Average # UMIs per nucleus prior to downsampling:", comma(pre_umis))

  pdf(file=paste0("figures/", data_type, "/", expr_type, "/cell_percent_nearest_neighbor_", toupper(sim_type), "_similarity_UMIs_downsampled_", data_type, "_", expr_type, "_", n_distinct(df$Dataset), "_datasets.pdf"), width=8, height=7)
  
  ggplot(df, aes(x=Threshold, y=Percent)) +
    geom_segment(
      aes(x=Threshold, xend=Threshold, y=0, yend=Percent)
    ) +
    geom_point(size=2) +
    theme_minimal() +
    theme(
      plot.title=element_text(hjust=0, face="bold", size=14, margin=margin(b=5)),
      plot.subtitle=element_text(hjust=0, size=11, lineheight=1.1, margin=margin(t=0, b=7)),
      axis.title.x=element_text(size=12, color="black", margin=margin(t=12)),
      axis.text.x=element_text(angle=45, hjust=1, color="black"),
      axis.title.y=element_text(size=12, color="black", margin=margin(r=15)),
      panel.border=element_rect(colour="black", fill=NA),
      plot.margin=unit(c(1, 1, 1, 1), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    xlab("# UMIs per nucleus") + ylab("%") +
    scale_y_continuous(breaks=seq(0, 36, by=10), limits=c(0, 36)) +
    facet_grid(Plot_Label~.)
  
  dev.off()
  
}

# plot_nn_sim <- function(datinfo, df, sim_type){
# 
#   idx <- df$Threshold%%100>0
#   pre_umis <- mean(df$Threshold[idx])
#   df$Threshold[idx] <- 1e6
#   
#   idx <- order(unique(df$Threshold))
#   df$Threshold <- comma(df$Threshold)
#   df$Threshold[df$Threshold=="1,000,000"] <- "Full coverage"
#   df$Threshold <- factor(
#     df$Threshold, levels=unique(df$Threshold)[idx]
#   )
#   
#   df <- merge(df, datinfo, by="Dataset")
#   
#   df$Plot_Label <- gsub("SSv4", "", df$Plot_Label)
#   df$Plot_Label <- gsub("FC", "", df$Plot_Label)
#   df$Plot_Label <- gsub("TC", "", df$Plot_Label)
# 
#   plot_title <- "Cell Nearest Neighbor Similarity vs. Read Depth"
#   plot_sub <- paste("Average median # UMIs per nucleus prior to downsampling:", comma(pre_umis))
#   
#   pdf(file=paste0("figures/", data_type, "/", expr_type, "/cell_nearest_neighbor_", toupper(sim_type), "_similarity_UMIs_downsampled_", data_type, "_", expr_type, "_", n_distinct(df$Dataset), "_datasets.pdf"), width=8.5, height=8)
#   
#   if(sim_type=="pearson"){
#     
#     y_lab <- paste(upper_first(sim_type), "Correlation")
#     
#   } else if(sim_type=="prop"){
#     
#     y_lab <- "Proportionality"
#     
#   } else {
#     
#     y_lab <- "Cosine Similarity"
#     
#   }
#     
#   print(
#     ggplot(df, aes(x=Threshold, y=Sim, fill=Plot_Label)) +
#       geom_violin(trim=T, size=.5, position=position_dodge(0.9), show.legend=F) +
#       geom_boxplot(
#         notch=T, width=.1, size=.3, fill="white", color="black", 
#         position=position_dodge(0.9), outlier.shape=NA, show.legend=F
#       ) +
#       theme_classic() + 
#       theme(
#         plot.title=element_text(hjust=.5, face="bold", size=15, margin=margin(b=10)),
#         plot.subtitle=element_text(hjust=.5, size=11, lineheight=1.1, margin=margin(t=4, b=8)),
#         legend.title=element_blank(),
#         axis.title.x=element_text(size=13, color="black",  margin=margin(b=10)),
#         axis.text.x=element_text(size=10, angle=45, hjust=1, color="black"), 
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_text(size=12, color="black", margin=margin(r=15)),
#         axis.text.y=element_text(color="black"),
#         axis.ticks.y=element_line(size=.4),
#         panel.grid.major.y=element_line(colour="darkgrey", size=.2),
#         plot.margin=unit(c(2, 2, 2, 2), "cm")
#       ) +
#       labs(title=plot_title, subtitle=plot_sub) +
#       xlab("# UMIs per nucleus") +
#       scale_y_continuous(name=y_lab, breaks=seq(0, 1, by=.2)) +
#       scale_fill_manual(values=brewer.pal(5, "Set2")) +
#       facet_grid(Plot_Label~.) +
#       theme(strip.background=element_blank())
#   )
#   
#   dev.off()
#   
# }
