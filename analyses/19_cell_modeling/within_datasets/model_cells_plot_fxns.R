library(ggplot2)
library(dplyr)
library(scales)
library(ggrastr)
library(stringr)
library(ggrepel)
library(ggpubr) ## ggarrange()

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")

r2_scatterplot <- function(datinfo, df, data_type, expr_type, summary_type, n_pred, n_genes){
  
  n_nuclei <- n_distinct(paste(df$Dataset, df$Cell_ID))
  datinfo$Plot_Label <- gsub("FC ", "", datinfo$Plot_Label)
  datinfo <- datinfo %>% dplyr::select(
    Dataset, Plot_Label, Study, 
    Author_Median_Unique_Genes_PC, 
    Author_No.Cell_Types
  )
  
  df$Distro <- as.character(df$Distro)
  df$Distro[grep("Rand", df$Distro)] <- "Random"
  df$Distro[grep("Real", df$Distro)] <- "Real"
  df$Distro <- factor(df$Distro, levels=c("Real", "Random"))
  
  df <- merge(df, datinfo, by="Dataset")

  df_temp <- df %>%
    dplyr::group_by(Plot_Label) %>%
    dplyr::summarise(
      R2_Real=mean(R2[Distro=="Real"]), 
      R2_Random=mean(R2[Distro=="Random"]),
      Median_Genes=unique(Author_Median_Unique_Genes_PC),
      No.CTs=unique(Author_No.Cell_Types),
      Study=unique(Study)
    ) %>%
    dplyr::mutate(
      R2_Diff=R2_Real-R2_Random
    )
  
  plot_sub <- paste("Multiple linear regression performed with", n_pred, "cell class", summary_type, "expression vectors for", comma(n_nuclei), "cells from", n_distinct(df$Dataset), "datasets") 
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/modeling_R2_scatterplot_", n_pred, "_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", n_nuclei, "_nuclei_", n_distinct(df$Dataset), "_datasets_", n_distinct(df$Study), "_studies_", n_genes, "_genes.pdf"), width=11, height=8)
  
  plot_title <- paste("Median Unique Genes Expressed vs. Variance Explained by Major Cell Classes")
  
  r2 <- signif(summary(lm(df_temp$R2_Real~df_temp$Median_Genes))$adj.r.squared, 2)
  
  print(
    ggplot(df_temp, aes(Median_Genes, y=R2_Real, color=Study)) +
      geom_point(show.legend=F) +
      geom_smooth(method='lm', formula= y~x, color="black", alpha=.2, show.legend=F) +
      geom_text_repel(aes(label=Plot_Label), box.padding=.4, show.legend=F, max.overlaps=100) +
      theme_minimal() +
      theme(
        plot.title=element_text(hjust=0, face="bold", size=15, lineheight=1.1, margin=margin(b=7)),
        plot.subtitle=element_text(hjust=0, size=11, lineheight=1.1, margin=margin(b=12)),
        legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=13, color="black", margin=margin(t=18)),
        axis.text.x=element_text(size=12, color="black"),
        axis.title.y=element_text(size=13, color="black", margin=margin(r=15)),
        axis.text.y=element_text(size=9, color="black"),
        panel.border=element_rect(colour="black", fill=NA),
        plot.margin=unit(c(1, 2, 1, 1), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab("Median # Unique Genes per Nucleus") +
      ylab(expression("Adjusted"~R^2)) +
      scale_x_continuous(labels=comma) +
      scale_color_manual(values=brewer_fxn(n_distinct(df_temp$Study))) +
      coord_cartesian(expand=F, clip="off") +
      geom_label(
        aes(label=paste("R^2 ==", r2), x=950, y=.651),
        label.padding=unit(.3, "lines"),
        color="black", size=5, parse=T, label.size=.5, 
      ) 
  )
  
  ## Plot delta between real and random vs. read depth
  
  plot_title <- paste("Median Unique Genes Expressed vs. Variance Explained by Real and Random Predictors")
  
  r2 <- signif(summary(lm(df_temp$R2_Diff~df_temp$Median_Genes))$adj.r.squared, 2)
  
  print(
    ggplot(df_temp, aes(Median_Genes, y=R2_Diff, color=Study)) +
      geom_point(show.legend=F) +
      geom_smooth(method='lm', formula= y~x, color="black", alpha=.2, show.legend=F) +
      geom_text_repel(aes(label=Plot_Label), box.padding=.4, show.legend=F, max.overlaps=100) +
      theme_minimal() +
      theme(
        plot.title=element_text(hjust=0, face="bold", size=13, lineheight=1.1, margin=margin(b=7)),
        plot.subtitle=element_text(hjust=0, size=11, lineheight=1.1, margin=margin(b=12)),
        legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=13, color="black", margin=margin(t=18)),
        axis.text.x=element_text(size=12, color="black"),
        axis.title.y=element_text(size=13, color="black", margin=margin(r=15)),
        axis.text.y=element_text(size=9, color="black"),
        panel.border=element_rect(colour="black", fill=NA),
        plot.margin=unit(c(1, 2, 1, 1), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab("Median # Unique Genes per Nucleus") +
      ylab(expression("Real"~-~"Random Adjusted"~R^2)) +
      scale_x_continuous(labels=comma) +
      scale_color_manual(values=brewer_fxn(n_distinct(df_temp$Study))) +
      coord_cartesian(expand=F, clip="off") +
      geom_label(
        aes(label=paste("R^2 ==", r2), x=950, y=.23, fontface="italic"),
        label.padding=unit(.3, "lines"),
        color="black", size=5, parse=T, label.size=.5,
      )
  )
  
   dev.off()
  
}

plot_r2_real_vs_rand <- function(datinfo, df, data_type, expr_type, summary_type, n_pred, n_genes, order_by=c("real_r2", "delta")){
  
  n_nuclei <- n_distinct(paste(df$Dataset, df$Cell_ID))
  datinfo$Plot_Label <- gsub("FC ", "", datinfo$Plot_Label)
  datinfo <- datinfo %>% dplyr::select(Dataset, Plot_Label, Study)
  
  df$Distro <- as.character(df$Distro)
  df$Distro[grep("Rand", df$Distro)] <- "Random"
  df$Distro[grep("Real", df$Distro)] <- "Real"
  df$Distro <- factor(df$Distro, levels=c("Random", "Real"))
  
  df <- merge(df, datinfo, by="Dataset")
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/modeling_real_vs_random_R2_boxplot_", n_pred, "_", summary_type, "cell_predictors_order_by_", order_by, "_", data_type, "_", expr_type, "_", n_nuclei, "_nuclei_", n_distinct(df$Dataset), "_datasets_", n_distinct(df$Study), "_studies_", n_genes, "_genes.pdf"), width=13, height=7.5)
  
  plot_title <- paste("Variance Explained by", n_pred, "Major Cell Classes")
  plot_sub <- paste("Multiple linear regression performed with cell class", summary_type, "expression predictors for", comma(n_nuclei), "cells from", n_distinct(df$Dataset), "datasets") 
  
  if(order_by=="real_r2"){
    
    temp <- df %>%
      dplyr::filter(Distro=="Real") %>%
      dplyr::group_by(Plot_Label) %>%
      dplyr::summarise(n=median(R2)) %>%
      arrange(n)
    
  } else {
    
    temp <- df %>%
      dplyr::group_by(Plot_Label) %>%
      dplyr::summarise(
        n=median(R2[Distro=="Real"])-median(R2[Distro=="Random"])
      ) %>%
      arrange(n)
    
  }
  
  df$Plot_Label <- factor(df$Plot_Label, levels=temp$Plot_Label)
  
  print(
    ggplot(df, aes(x=Plot_Label, y=R2, group=interaction(Plot_Label, Distro), fill=Distro)) + 
      geom_violin(position=position_dodge(.7), linewidth=.5, fill="white", color="grey", show.legend=F) + 
      geom_boxplot(
        notch=T, linewidth=.5, alpha=1, width=.3, outlier.shape=NA, 
        position=position_dodge(.7), key_glyph=draw_key_rect
      ) + 
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=0, face="bold", size=17),
        plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.box.just="left",
        legend.key=element_rect(color="black", linewidth=1),
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=11, angle=45, color="black", hjust=1, margin=margin(b=10)),
        axis.ticks.x=element_line(linewidth=.5),
        axis.title.y=element_text(size=13, margin=margin(r=13)),
        axis.text.y=element_text(color="black", size=10),
        panel.border=element_rect(colour="black", fill=NA),
        panel.grid.major.x=element_blank(),
        plot.margin=unit(c(1, 1, 1, 1), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) +
      ylab(expression("Adjusted"~R^2)) +
      scale_y_continuous(breaks=seq(0, 1, by=.1)) +
      scale_fill_manual(values=rev(brewer.pal(3, "Set2")[c(1,2)])) +
      guides(fill=guide_legend(title="Predictors", keywidth=.7, keyheight=.7))
  )
  
  dev.off()
  
}

plot_r2_predictor_groups <- function(datinfo, df, data_type, expr_type, summary_type, n_genes){
  
  df <- df %>% dplyr::filter(Distro=="Real_R2")
  n_nuclei <- n_distinct(paste(df$Dataset, df$Cell_ID))
  datinfo$Plot_Label <- gsub("FC ", "", datinfo$Plot_Label)
  datinfo <- datinfo %>% dplyr::select(
    Dataset, Plot_Label, Study, 
    Author_Median_Unique_Genes_PC
  )
  
  df <- merge(df, datinfo, by="Dataset")
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/modeling_R2_boxplot_", summary_type, "cell_predictors_", data_type, "_", expr_type, "_", n_nuclei, "_nuclei_", n_distinct(df$Dataset), "_datasets_", n_distinct(df$Study), "_studies_", n_genes, "_genes.pdf"), width=13, height=7.5)
  
  ## Plot distribution of R2 by predictor group per dataset:
  
  plot_title <- paste0("Variance Explained by Major Cell Classes")
  plot_sub <- paste("Multiple linear regression performed with cell class", summary_type, "expression predictors for", comma(n_nuclei), "cells from", n_distinct(df$Dataset), "datasets") 
  
  temp <- df %>%
    dplyr::group_by(Plot_Label) %>%
    dplyr::summarise(n=median(R2)) %>%
    arrange(n)
  
  df$Plot_Label <- factor(df$Plot_Label, levels=temp$Plot_Label)
  
  temp <- unique(df$Predictors)
  
  df$Predictors <- factor(df$Predictors, levels=temp[order(nchar(temp))])
  
  pred1 <- paste(c("ASC", "MIC", "OG", "OPC", "NEU"), collapse=" | ")
  pred2 <- paste(c("ASC", "MIC", "OG", "OPC", "EXC", "INH"), collapse=" | ")
  pred3 <- paste(c("ASC", "MIC", "OG", "OPC", "EXC", "CGE", "MGE"), collapse=" | ")
  
  if(n_distinct(df$Predictors)>2){
    lims <- c(pred1, pred2, pred3)
  } else {
    lims <- c(pred1, pred2)
  }
  
  print(
    ggplot(df, aes(x=Plot_Label, y=R2, group=interaction(Plot_Label, Predictors), fill=Predictors)) + 
      geom_violin(position=position_dodge(.7), linewidth=.5, fill="white", color="grey", show.legend=F) + 
      geom_boxplot(
        notch=T, linewidth=.5, alpha=1, width=.3, outlier.shape=NA, 
        position=position_dodge(.7), key_glyph=draw_key_rect
      ) + 
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=0, face="bold", size=17),
        plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.position="bottom",
        legend.direction="vertical",
        legend.box.just="left",
        legend.key=element_rect(color="black", linewidth=1), #legend.spacing.y=unit(1, "cm"),
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=11, angle=45, color="black", hjust=1),
        axis.ticks.x=element_line(linewidth=.5),
        axis.title.y=element_text(size=13, margin=margin(r=13)),
        axis.text.y=element_text(color="black", size=10),
        panel.border=element_rect(colour="black", fill=NA),
        panel.grid.major.x=element_blank(),
        plot.margin=unit(c(1, 1, 1, 1), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) +
      ylab(expression("Adjusted"~R^2)) +
      scale_y_continuous(breaks=seq(0, 1, by=.1)) +
      scale_fill_manual(values=rev(brewer.pal(3, "Set2")), limits=lims) +
      guides(fill=guide_legend(title="Predictors", keywidth=.7, keyheight=.7))
  )
  
  ## Plot barplots of R2 per dataset:
  
  plot_title <- paste0("Variance Explained by Major Cell Class Predictors")
  
  if(n_distinct(df$Predictors)>2){
    
    df_temp <- df %>% 
      dplyr::group_by(Plot_Label) %>%
      dplyr::summarise(
        R2_Pred1=mean(R2[Predictors==pred1]), 
        R2_Pred2=mean(R2[Predictors==pred2]),
        R2_Pred3=mean(R2[Predictors==pred3])
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
    
  } else { ## if(n_distinct(df$Predictors)>2){
    
    df_temp <- df %>% 
      dplyr::group_by(Plot_Label) %>%
      dplyr::summarise(
        R2_Pred1=mean(R2[Predictors==pred1]), 
        R2_Pred2=mean(R2[Predictors==pred2])
      ) %>%
      dplyr::mutate(
        R2_Diff1_2=R2_Pred2-R2_Pred1
      ) %>%
      dplyr::select(-c(R2_Pred2)) %>%
      reshape2::melt(
        value.name="R2", 
        variable.name="Predictors"
      )
    
  } ## if(n_distinct(df$Predictors)>2){} else {
  
  df_temp$Predictors <- as.character(df_temp$Predictors)
  df_temp$Predictors[grep("Pred1", df_temp$Predictors)] <- pred1
  df_temp$Predictors[grep("Diff1_2", df_temp$Predictors)] <- pred2
  df_temp$Predictors[grep("Diff2_3", df_temp$Predictors)] <- pred3
  
  temp <- df_temp %>%
    dplyr::group_by(Plot_Label) %>%
    dplyr::summarise(
      n=sum(R2)
    ) %>%
    arrange(n)
  
  df_temp$Plot_Label <- factor(df_temp$Plot_Label, levels=temp$Plot_Label)
  
  if(n_distinct(df$Predictors)>2){
    
    df_temp$Predictors <- factor(df_temp$Predictors, levels=c(pred3, pred2, pred1))
    
  } else {
    
    df_temp$Predictors <- factor(df_temp$Predictors, levels=c(pred2, pred1))
    
  }
  
  print(
    ggplot(df_temp, aes(x=Plot_Label, y=R2, fill=Predictors)) + 
      geom_bar(stat="identity", position="stack", linewidth=.4, color="black") +
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=0, face="bold", size=17),
        plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.position="bottom",
        legend.direction="vertical",
        legend.box.just="left",
        legend.key=element_rect(color="black", linewidth=.1),
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=11, angle=45, color="black", hjust=1),
        axis.ticks.x=element_line(linewidth=.5),
        axis.title.y=element_text(size=13, margin=margin(r=13)),
        axis.text.y=element_text(color="black", size=10),
        panel.border=element_rect(colour="black", fill=NA),
        panel.grid.major.x=element_blank(),
        plot.margin=unit(c(1, 4, 1, 4), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) +
      ylab(expression("Adjusted"~R^2)) +
      scale_fill_manual(values=rev(brewer.pal(3, "Set2")), limits=lims) +
      scale_y_continuous(breaks=seq(0, 1, by=.05)) +
      guides(fill=guide_legend(title="Predictors", keywidth=.8, keyheight=.8))
  )
  
  dev.off()
  
}

