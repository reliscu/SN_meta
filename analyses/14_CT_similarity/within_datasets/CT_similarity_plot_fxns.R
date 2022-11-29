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

CT_similarity_boxplot <- function(datinfo){
  
  datinfo$Plot_Label <- gsub("FC ", "", datinfo$Plot_Label)
  datinfo <- datinfo %>% dplyr::select(Dataset, Plot_Label, Study)
  
  df <- do.call(rbind, lapply(sim_list, function(x) do.call(rbind, x)))
  rownames(df) <- 1:nrow(df)

  df <- df %>% 
    # dplyr::filter(
    #   Cell_Type==Distro|MO_Cell_Class==Distro
    # ) %>%
    dplyr::mutate(
      Distro=if_else(Cell_Type==Distro, "Cell Type", "Cell Class")
    )

  df <- merge(df, datinfo, by="Dataset")
  
  temp <- df %>%
    dplyr::group_by(Plot_Label) %>%
    dplyr::summarise(
      n=median(Similarity[Distro=="Cell Type"])
    ) %>%
    arrange(n)

  df$Plot_Label <- factor(df$Plot_Label, levels=temp$Plot_Label)

  pdf(file=paste0("figures/", data_type, "/", expr_type, "/CT_vs_cell_class_", summary_type, "cell_similarity_", data_type, "_", expr_type, "_", n_distinct(df$Dataset), "_datasets_", n_genes, "_genes.pdf"), width=12, height=9)
  
  plot_title <- "Cell Type vs. Cell Class Similarity"
  plot_sub <- paste0("Pearson similarity calculated between all cell class nuclei and cell type vs. cell class mean expression vectors")
  
  print(
    ggplot(df, aes(x=Plot_Label, y=Similarity, group=interaction(Plot_Label, Distro), fill=Distro)) + 
      geom_violin(position=position_dodge(.7), linewidth=.5, fill="white", color="grey", show.legend=F) + 
      geom_boxplot(linewidth=.5, alpha=1, width=.3, outlier.shape=NA, position=position_dodge(.7), key_glyph=draw_key_rect) + 
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=0, face="bold", size=17),
        plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
        legend.title=element_blank(),
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
      ylab("Pearson Similarity") +
      scale_y_continuous(breaks=seq(0, 1, by=.1)) +
      scale_fill_manual(values=rev(brewer.pal(3, "Set2")[c(1,2)])) 
  )
  
  temp <- df %>%
    dplyr::group_by(MO_Cell_Class) %>%
    dplyr::summarise(
      n=median(Similarity[Distro=="Cell Type"])
    ) %>%
    arrange(n)
  
  df$MO_Cell_Class <- factor(df$MO_Cell_Class, levels=temp$MO_Cell_Class)
  
  print(
    ggplot(df, aes(x=MO_Cell_Class, y=Similarity, group=interaction(MO_Cell_Class, Distro), fill=Distro)) + 
      geom_violin(position=position_dodge(.7), linewidth=.5, fill="white", color="grey", show.legend=F) + 
      geom_boxplot(linewidth=.5, alpha=1, width=.3, outlier.shape=NA, position=position_dodge(.7), key_glyph=draw_key_rect) + 
      theme_minimal() + 
      theme(
        plot.title=element_text(hjust=0, face="bold", size=17),
        plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.box.just="left",
        legend.key=element_rect(color="black", linewidth=1),
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=12, color="black", margin=margin(t=5, b=10)),
        axis.ticks.x=element_line(linewidth=.5),
        axis.title.y=element_text(size=13, margin=margin(r=13)),
        axis.text.y=element_text(color="black", size=10),
        panel.border=element_rect(colour="black", fill=NA),
        panel.grid.major.x=element_blank(),
        plot.margin=unit(c(1, 1, 1, 1), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) +
      ylab("Pearson Similarity") +
      scale_y_continuous(breaks=seq(0, 1, by=.1)) +
      scale_fill_manual(values=rev(brewer.pal(3, "Set2")[c(1,2)])) +
      guides(fill=guide_legend(keywidth=.7, keyheight=.7))
  )

  df_temp <- df %>%
    dplyr::group_by(Plot_Label, MO_Cell_Class) %>%
    dplyr::summarise(
      Percent=sum(
        Similarity[Distro=="Cell Type"]>max(Similarity[Distro=="Cell Class"])
      )/n()*100,
      No.CTs=n_distinct(Cell_Type)
    ) %>% 
    arrange(Percent)
  
  temp <- df_temp %>% 
    dplyr::group_by(Plot_Label) %>%
    dplyr::summarise(n=sum(Percent)) %>%
    arrange(n)
  
  df_temp$Plot_Label <- factor(df_temp$Plot_Label, levels=temp$Plot_Label)
  
  plot_title <- paste0("% Nuclei with Greater Similarity to Cell Type than Cell Class")
  
  p <- ggplot(df_temp, aes(x=Plot_Label, y=Percent, fill=MO_Cell_Class)) +
    geom_bar(stat="identity", position="stack", linewidth=.4, color="black") +
    theme_minimal() +
    theme(
      plot.title=element_text(hjust=0, face="bold", size=14),
      plot.subtitle=element_text(hjust=0, size=10, lineheight=1, margin=margin(t=4, b=11)),
      legend.title=element_blank(),
      legend.text=element_text(size=10),
      legend.position="bottom",
      legend.direction="horizontal",
      legend.box.just="left",
      legend.key=element_rect(color="black", linewidth=.1),
      axis.title.x=element_blank(),
      axis.text.x=element_text(size=11, angle=45, color="black", hjust=1),
      axis.ticks.x=element_line(linewidth=.5),
      axis.title.y=element_text(size=13, margin=margin(r=15)),
      axis.text.y=element_text(color="black", size=10),
      panel.border=element_rect(colour="black", fill=NA),
      panel.grid.major.x=element_blank(),
      plot.margin=unit(c(1, 4, 1, 4), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    ylab("%") +
    scale_fill_manual(values=rev(brewer.pal(10, "Paired"))) +
    guides(fill=guide_legend(keywidth=.8, keyheight=.8, nrow=1))
  
  print(p)
  print(p + scale_y_continuous(trans="sqrt"))
    
  dev.off()
  
}



