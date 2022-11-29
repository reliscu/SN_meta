library(ggplot2)
library(dplyr)
#library(viridis)
library(scales)

plot_eigencell_VE <- function(datinfo, df){
  
  pre_umis <- mean(df$Threshold[which(df$Threshold%%100>0)])
  df$Threshold[df$Threshold%%100>0] <- 1e6
  idx <- order(
    unique(df$Threshold)
  )

  df$Threshold <- comma(df$Threshold)
  df$Threshold[df$Threshold=="1,000,000"] <- "Full coverage"
  df$Threshold <- factor(
    df$Threshold, levels=unique(df$Threshold)[idx]
  )
  
  df <- merge(df, datinfo, by="Dataset")
  
  df$Plot_Label <- gsub("SSv4", "", df$Plot_Label)
  df$Plot_Label <- gsub("SSv4", "", df$Plot_Label)
  df$Plot_Label <- gsub("FC", "", df$Plot_Label)
  df$Plot_Label <- gsub("TC", "", df$Plot_Label)

  plot_sub <- paste("Average median # UMIs per nucleus prior to downsampling:", comma(pre_umis))
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/eigencell_VE_violin_plot_UMIs_downsampled_", data_type, "_", expr_type, "_", n_distinct(df$Dataset), "_datasets.pdf"), width=8.5, height=8)

  print(
    ggplot(df, aes(x=Threshold, y=PC1_VE, fill=Plot_Label)) +
      geom_violin(trim=T, size=.5, position=position_dodge(0.9), show.legend=F) +
      geom_boxplot(
        notch=T, width=.2, size=.3, fill="white", color="black", 
        position=position_dodge(0.9), outlier.shape=NA, show.legend=F
      ) +
      theme_classic() + 
      theme(
        plot.title=element_text(hjust=.5, face="bold", size=15, margin=margin(b=10)),
        plot.subtitle=element_text(hjust=.5, size=11, lineheight=1.1, margin=margin(t=4, b=8)),
        legend.title=element_blank(),
        axis.title.x=element_text(size=13, color="black",  margin=margin(b=10)),
        axis.text.x=element_text(size=10, angle=45, hjust=1, color="black"), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=12, color="black", margin=margin(r=15)),
        axis.text.y=element_text(color="black"),
        axis.ticks.y=element_line(size=.4),
        panel.grid.major.y=element_line(colour="darkgrey", size=.2),
        plot.margin=unit(c(2, 2, 2, 2), "cm")
      ) +
      labs(
        title="Eigencell % Variance Explained vs. Read Depth", subtitle=plot_sub
      ) +
      xlab("# UMIs per Nucleus") +
      ylab("% Variance Explained") +
      scale_y_continuous(
        breaks=seq(0, 1, by=.2), 
        limits=c(0, 1), 
        labels=seq(0, 100, by=20),
        minor_breaks=
      ) +
      scale_fill_manual(values=brewer.pal(5, "Set2")) +
      facet_grid(Plot_Label~.) +
      theme(strip.background=element_blank())
  )
  
  dev.off()
  
}
