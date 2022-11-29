library(ggplot2)
library(dplyr)
library(viridis)
library(scales)
#library(ggrastr)
library(stringr)
library(ggpubr) ## ggarrange()

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

plot_real_vs_rand_r2 <- function(df, data_type, expr_type, summary_type, n_pred, n_genes){
  
  n_nuclei <- n_distinct(paste(df$Dataset, df$Cell_ID))
  mean_r2 <- signif(mean(df$R2[df$Distro=="Real_R2"]), 2)
  
  df$Distro <- as.character(df$Distro)
  df$Distro[grep("Random", df$Distro)] <- "Random Predictors"
  df$Distro[grep("Real", df$Distro)] <- "Cell Class Predictors"

  plot_title <- paste0(mean_r2*100, "% of the Variance in ", comma(n_nuclei), " Cells\nis Explained by ", n_pred, " Cell Class Predictors")
  plot_sub <- paste("Multiple linear regression performed with", summary_type, "expression vectors of cell class cells\nor randomly selected cells over", comma(n_genes), "protein coding genes")
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/modeling_R2_", n_pred, "_real_vs_random_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", n_nuclei, "_nuclei_", n_distinct(df$Dataset), "_datasets_", n_genes, "_genes.pdf"), width=8.5, height=7)
  
  print(
    ggplot(df, aes(x=Distro, y=R2, fill=Distro)) + 
      geom_violin(linewidth=1.1, show.legend=F) + 
      geom_boxplot(notch=T, width=.05, color="black", fill="white", outlier.shape=NA, show.legend=F) +
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=0, face="bold", size=14),
        plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=8)),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        legend.position="bottom",
        legend.box="horizontal",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=11, color="black", margin=margin(t=15)),
        axis.ticks.x=element_line(size=0),
        axis.title.y=element_text(size=13, color="black", margin=margin(r=15)),
        panel.border=element_rect(colour="black", fill=NA),
        plot.margin=unit(c(1, 2, 1, 1), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) +
      ylab(expression("Adjusted"~R^2)) +
      scale_fill_manual(values=brewer.pal(3, "Set2"))
  )
  
  dev.off()
  
}

plot_real_vs_rand_r2_seq <- function(df, data_type, expr_type, summary_type, n_genes){
  
  n_nuclei <- n_distinct(paste(df$Dataset, df$Cell_ID))
  
  df$Distro <- "Cell Class"
  df$Distro[grep("Random", df$Predictor)] <- "Random Predictors"
  df$Predictors <- as.character(df$Predictors)
  df$Predictors <- gsub("_", " | ", gsub("Random_", "", df$Predictors), fixed=T)
  temp <- unique(df$Predictors)
  df$Predictors <- factor(df$Predictors, levels=temp[order(nchar(temp))])
  
  plot_title <- paste0("Variance Explained by Additional Cell Class Predictors")
  plot_sub <- paste("Multiple linear regression performed with cell class", summary_type, "expression predictors for", comma(n_nuclei), "cells from", n_distinct(df$Dataset), "datasets") 
  
  df_temp <- df %>% 
    dplyr::filter(Distro=="Cell Class") %>%
    dplyr::summarise(
      R2_Pred1=mean(R2[Predictors==levels(df$Predictors)[1]]), 
      R2_Pred2=mean(R2[Predictors==levels(df$Predictors)[2]]),
      R2_Pred3=mean(R2[Predictors==levels(df$Predictors)[3]]),
      R2_Pred4=mean(R2[Predictors==levels(df$Predictors)[4]]), 
      R2_Pred5=mean(R2[Predictors==levels(df$Predictors)[5]]),
      R2_Pred6=mean(R2[Predictors==levels(df$Predictors)[6]]),
      R2_Pred7=mean(R2[Predictors==levels(df$Predictors)[7]]),
      R2_Pred8=mean(R2[Predictors==levels(df$Predictors)[8]])
    ) %>%
    dplyr::mutate(
      R2_Diff1=R2_Pred1,
      R2_Diff1_2=R2_Pred2-R2_Pred1,
      R2_Diff2_3=R2_Pred3-R2_Pred2,
      R2_Diff3_4=R2_Pred4-R2_Pred3,
      R2_Diff4_5=R2_Pred5-R2_Pred4,
      R2_Diff5_6=R2_Pred6-R2_Pred5,
      R2_Diff6_7=R2_Pred7-R2_Pred6,
      R2_Diff7_8=R2_Pred8-R2_Pred7
    ) %>%
    dplyr::select(contains("Diff")) %>%
    reshape2::melt(
      value.name="R2", 
      variable.name="Predictors"
    )
  
  df_temp$Predictors <- as.character(df_temp$Predictors)
  df_temp$Predictors <- factor(levels(df$Predictors), levels=rev(levels(df$Predictors)))
  color_df <- data.frame(Predictors=df_temp$Predictors, Color=brewer.pal(8, "Set1"))
  df <- merge(df, color_df, by="Predictors")
  df$Color[df$Distro=="Random Predictors"] <- "#A9A9A9"
  df$Color <- factor(df$Color, levels=c(color_df$Color, "#A9A9A9"))

  p1 <- ggplot(df, aes(x=Predictors, y=R2, group=interaction(Predictors, Distro), fill=Color)) + 
    geom_violin(position=position_dodge(.7), linewidth=.5, fill="white", color="grey", show.legend=F) + 
    geom_boxplot(notch=T, linewidth=.5, alpha=1, width=.3, outlier.shape=NA, position=position_dodge(.7), key_glyph=draw_key_rect) + 
    scale_fill_identity() +
    theme_minimal() + 
    theme(
      plot.title=element_text(hjust=0, face="bold", size=17),
      plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
      legend.title=element_text(size=12),
      legend.text=element_text(size=10),
      legend.position="bottom",
      legend.direction="vertical",
      legend.box.just="left",
      legend.key=element_rect(color="black", linewidth=1), 
      axis.title.x=element_blank(), 
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_text(size=13, margin=margin(r=13)),
      axis.text.y=element_text(color="black", size=10),
      panel.grid.major.x=element_blank(),
      plot.margin=unit(c(1, 0, 0, 1), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    ylab(expression("Adjusted"~R^2)) +
    scale_y_continuous(breaks=seq(0, 1, by=.1)) +
    scale_fill_identity(guide="legend", labels=c(as.character(df_temp$Predictors), "Random Predictors")) +
    guides(fill=guide_legend(title="Predictors", keywidth=.7, keyheight=.7))
  
  p2 <- ggplot(df_temp, aes(x=1, y=R2, fill=Predictors)) + 
    geom_bar(stat="identity", linewidth=.4, color="black", position="stack", show.legend=F) +
    theme_minimal() + 
    theme(
      plot.title=element_blank(),
      plot.subtitle=element_blank(), 
      legend.title=element_text(size=12),
      legend.text=element_text(size=10),
      legend.position="bottom",
      legend.direction="vertical",
      legend.box.just="left",
      legend.key=element_rect(color="black", linewidth=.1),
      axis.title.x=element_blank(), 
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y.right=element_text(size=13, margin=margin(l=13)),
      axis.text.y=element_blank(), 
      panel.border=element_blank(), 
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.margin=unit(c(2.2, 0, 0, 0), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    ylab("Proportion of Variance Explained") +
    scale_fill_manual(values=brewer.pal(8, "Set1"), limits=levels(df$Predictors)) +
    scale_color_manual(values=brewer.pal(8, "Set1"), limits=levels(df$Predictors)) +
    scale_y_continuous(breaks=seq(0, 1, by=.1))
  
  p3 <- ggplot(subset(df_temp, grepl("MIC", Predictors)), aes(x=1, y=R2, fill=Predictors)) + 
    geom_bar(stat="identity",  linewidth=.4, color="black", position="stack", show.legend=F) +
    theme_minimal() + 
    theme(
      plot.title=element_blank(),
      plot.subtitle=element_blank(), 
      legend.title=element_text(size=12),
      legend.text=element_text(size=10),
      legend.position="bottom",
      legend.direction="vertical",
      legend.box.just="left",
      legend.key=element_rect(color="black", linewidth=.1),
      axis.title.x=element_blank(), 
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(), 
      panel.border=element_blank(), 
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.margin=unit(c(2.2, 0, 0, 0), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    ylab("Proportion of Variance Explained") +
    scale_fill_manual(values=brewer.pal(8, "Set1"), limits=levels(df$Predictors)) +
    scale_color_manual(values=brewer.pal(8, "Set1"), limits=levels(df$Predictors)) +
    scale_y_continuous(breaks=seq(0, 1, by=.1))

  pdf(file=paste0("figures/", data_type, "/", expr_type, "/modeling_R2_boxplot_real_vs_random_", summary_type, "cell_sequential_predictors_", data_type, "_", expr_type, "_", n_nuclei, "_nuclei_", n_distinct(df$Dataset), "_datasets_", n_genes, "_genes.pdf"), width=14.5, height=9)
  
  print(ggarrange(p1, p2, p3, ncol=3, widths=c(2, .2, .165), common.legend=T, legend="bottom") + theme(plot.margin=margin(1, 1, 1, 1, "cm")))
  
  dev.off()
  
}

plot_r2_predictor_groups <- function(df, data_type, expr_type, summary_type, n_genes){
  
  df <- df %>% dplyr::filter(Distro=="Real_R2")
  n_nuclei <- n_distinct(paste(df$Dataset, df$Cell_ID))

  ## Plot distribution of R2 by predictor group per dataset:
  
  plot_title <- paste0("Variance Explained by Major Cell Classes")
  plot_sub <- paste("Multiple linear regression performed with cell class", summary_type, "expression predictors for", comma(n_nuclei), "cells from", n_distinct(df$Dataset), "datasets") 
  
  temp <- df %>%
    dplyr::group_by(Predictors) %>%
    dplyr::summarise(n=median(R2)) %>%
    arrange(n)
  
  df$Predictors <- factor(df$Predictors, levels=temp$Predictors)

  ## Barplot of proportion R2 per predictor group:
  
  df_temp <- df %>% 
    dplyr::summarise(
      R2_Pred1=mean(R2[Predictors==levels(df$Predictors)[1]]), 
      R2_Pred2=mean(R2[Predictors==levels(df$Predictors)[2]]),
      R2_Pred3=mean(R2[Predictors==levels(df$Predictors)[3]])
    ) %>%
    dplyr::mutate(
      R2_Diff1_2=R2_Pred2-R2_Pred1,
      R2_Diff2_3=R2_Pred3-R2_Pred2
    ) %>%
    dplyr::select(-c(R2_Pred2, R2_Pred3)) %>%
    reshape2::melt(
      value.name="R2", 
      variable.name="Predictors"
    )
  
  df_temp$Predictors <- as.character(df_temp$Predictors)
  df_temp$Predictors <- factor(levels(df$Predictors), levels=rev(levels(df$Predictors)))
  
  p1 <- ggplot(df, aes(x=Predictors, y=R2, fill=Predictors)) + 
    geom_violin(position=position_dodge(.7), linewidth=.5, fill="white", color="grey", show.legend=F) + 
    geom_boxplot(notch=T, linewidth=.5, alpha=1, width=.3, outlier.shape=NA, position=position_dodge(.7), key_glyph=draw_key_rect) + 
    theme_minimal() + 
    theme(
      plot.title=element_text(hjust=0, face="bold", size=17),
      plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
      legend.title=element_text(size=12),
      legend.text=element_text(size=10),
      legend.position="bottom",
      legend.direction="vertical",
      legend.box.just="left",
      legend.key=element_rect(color="black", linewidth=1),
      axis.title.x=element_blank(), 
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_text(size=13, margin=margin(r=13)),
      axis.text.y=element_text(color="black", size=10),
      panel.grid.major.x=element_blank(),
      plot.margin=unit(c(1, 0, 1, 1), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    ylab(expression("Adjusted"~R^2)) +
    scale_y_continuous(breaks=seq(0, 1, by=.1)) +
    scale_fill_manual(values=rev(brewer.pal(3, "Set2"))) +
    guides(fill=guide_legend(title="Predictors", keywidth=.7, keyheight=.7))

  p2 <- ggplot(df_temp, aes(x=1, y=R2, fill=Predictors)) + 
    geom_bar(stat="identity", position="stack", linewidth=.4, color="black", show.legend=F) +
    theme_minimal() + 
    theme(
      plot.title=element_blank(),
      plot.subtitle=element_blank(),
      legend.title=element_text(size=12),
      legend.text=element_text(size=10),
      legend.position="bottom",
      legend.direction="vertical",
      legend.box.just="left",
      legend.key=element_rect(color="black", linewidth=.1),
      axis.title.x=element_blank(), 
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y.right=element_text(size=13, margin=margin(l=13)),
      axis.text.y=element_blank(), #element_text(color="black", size=10)
      panel.border=element_blank(), #element_rect(colour="black", fill=NA),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.margin=unit(c(2.2, 2, .5, 0), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    ylab("Proportion of Variance Explained") +
    scale_fill_manual(values=rev(brewer.pal(3, "Set2")), limits=levels(df$Predictors)) +
    scale_y_continuous(breaks=seq(0, 1, by=.1)) +
    guides(fill=guide_legend(title="Predictors", keywidth=.8, keyheight=.8))
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/modeling_R2_boxplot_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", n_nuclei, "_nuclei_", n_distinct(df$Dataset), "_datasets_", n_genes, "_genes.pdf"), width=10, height=7)
  
  print(ggarrange(p1, p2, ncol=2, widths=c(2, .5), common.legend=T, legend="bottom") + theme(plot.margin=margin(1, 1, 1, 1, "cm")))

  dev.off()
  
}

plot_real_vs_rand_coeff <- function(df, data_type, expr_type, summary_type, n_pred, n_genes){
  
  n_nuclei <- n_distinct(paste(df$Dataset, df$Cell_ID))

  df_real <- df[,!grepl("Random", colnames(df))]
  df_real <- reshape2::melt(df_real, value.name="Coef", variable.name="Coef_Cell_Class")
  df_real <- df_real %>% dplyr::group_by(Coef_Cell_Class) %>% slice_sample(n=1e5)
  df_real$Coef_Cell_Class <- gsub("Real_", "", df_real$Coef_Cell_Class)
  df_real$Pred <- paste0(sapply(df_real$Coef_Cell_Class, cell_class_full_name), "\nPredictor")
  
  df_rand <- df[,!grepl("Real", colnames(df))]
  df_rand <- reshape2::melt(df_rand, value.name="Coef", variable.name="Coef_Cell_Class")
  df_rand$Coef_Cell_Class <- gsub("Random_", "", gsub(".[0-9]{1}$", "", df_rand$Coef_Cell_Class))
  df_rand <- df_rand %>% dplyr::group_by(Coef_Cell_Class) %>% slice_sample(n=1e5)
  df_rand$Pred <- paste0(sapply(df_rand$Coef_Cell_Class, cell_class_full_name), "\nRandom Predictor")
  
  plot_title <- paste("Cell Class Regression Coefficients")
  plot_sub <- paste("Multiple linear regression performed with", n_pred, summary_type, "expression vectors of cell class cells or randomly selected cells for", comma(n_nuclei), "nuclei from", n_distinct(df$Dataset), "datasets") 
  
  p_real <- ggplot(df_real, aes(x=MO_Cell_Class, y=Coef, color=ifelse(
    MO_Cell_Class==Coef_Cell_Class, 'magenta', 'grey')
  )) + 
    geom_violin(linewidth=.5, show.legend=F) + 
    geom_boxplot(width=.03, fill="white", outlier.shape=NA, show.legend=F) +
    theme_minimal() + 
    theme(
      plot.title=element_text(hjust=0, face="bold", size=15),
      plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=2.5, b=12)),
      legend.title=element_blank(),
      legend.text=element_text(size=10),
      legend.position="bottom",
      legend.box="horizontal",
      axis.title.x=element_blank(),
      axis.text.x=element_text(size=8, angle=45, color="black", margin=margin(t=10)),
      axis.ticks.x=element_line(size=0),
      axis.title.y=element_blank(),
      panel.border=element_rect(colour="black", fill=NA),
      plot.margin=unit(c(1, 2, 0, .5), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    xlab("Cell Class") + ylab("Regression Coefficient") +
    facet_wrap(.~Pred, nrow=1, scales="free_x") +
    theme(strip.text.x=element_text(size=10)) +
    scale_color_identity()

  p_rand <- ggplot(df_rand, aes(x=MO_Cell_Class, y=Coef, color=ifelse(
    MO_Cell_Class==Coef_Cell_Class, 'magenta', 'grey')
  )) + 
    geom_violin(linewidth=.5, show.legend=F) + 
    geom_boxplot(width=.03, fill="white", outlier.shape=NA, show.legend=F) +
    theme_minimal() + 
    theme(
      plot.title=element_blank(), 
      plot.subtitle=element_blank(), 
      axis.title.x=element_text(size=13, color="black", margin=margin(t=10)), 
      axis.text.x=element_text(size=8, angle=45, color="black", margin=margin(t=10)),
      axis.ticks.x=element_line(size=0),
      axis.title.y=element_blank(),
      panel.border=element_rect(colour="black", fill=NA),
      plot.margin=unit(c(0, 2, 1, .5), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    xlab("Cell Class") + ylab("Regression Coefficient") +
    facet_wrap(.~Pred, nrow=1, scales="free_x") +
    theme(strip.text.x=element_text(size=10)) +
    scale_color_identity()
  
  fig <- ggarrange(p_real, p_rand, nrow=2)
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/modeling_coefficients_", n_pred, "_real_vs_random_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", n_nuclei, "_nuclei_", n_distinct(df$Dataset), "_datasets_", n_genes, "_genes.pdf"), width=13, height=9)
  
  print(annotate_figure(fig, left=text_grob("Regression Coefficient", rot=90, vjust=.5)))
  
  dev.off()
  
}

plot_r2_left_out_dataset <- function(datinfo, df, data_type, expr_type, summary_type, n_pred, n_genes){

  datinfo <- datinfo %>% dplyr::select(Dataset, Plot_Label)
  df <- merge(df[,c(1:5)], datinfo, by="Dataset")
  
  temp <- df %>%
    dplyr::group_by(Plot_Label) %>%
    dplyr::summarise(n=median(R2)) %>%
    arrange(n)
  
  df$Plot_Label <- factor(df$Plot_Label, levels=temp$Plot_Label)
  
  plot_title <- paste("Variance Explained by", n_pred, "Major Cell Classes")
  plot_sub <- paste("Multiple linear regression performed on nuclei from left out dataset using", summary_type, "expression predictors built on remaining", n_distinct(df$Dataset)-1, "datasets") 
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/modeling_R2_left_out_dataset_", n_pred, "_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", n_distinct(df$Dataset), "_datasets_", n_genes, "_genes.pdf"), width=11.5, height=8)
  
  print(
    ggplot(df, aes(x=Plot_Label, y=R2)) + 
      geom_violin(position=position_dodge(.7), linewidth=.5, fill="white", color="grey", show.legend=F) + 
      geom_boxplot(
        notch=T, linewidth=.5, alpha=1, width=.3, outlier.shape=NA, 
        position=position_dodge(.7), key_glyph=draw_key_rect
      ) + 
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=0, face="bold", size=14),
        plot.subtitle=element_text(hjust=0, size=10, lineheight=1, margin=margin(t=4, b=8)),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        legend.position="bottom",
        legend.box="horizontal",
        axis.title.x=element_text(size=13, color="black", margin=margin(t=10)),
        axis.text.x=element_text(size=11, angle=45, hjust=1, color="black", margin=margin(b=15)),
        axis.ticks.x=element_line(size=0),
        axis.title.y=element_text(size=13, color="black", margin=margin(r=10)),
        panel.border=element_rect(colour="black", fill=NA),
        plot.margin=unit(c(1, 2, 1, 1), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) +
      xlab("Left Out Dataset") +
      ylab(expression("Adjusted"~R^2)) +
      scale_fill_manual(values=brewer.pal(3, "Set2"))
  )
  
  dev.off()
  
}
