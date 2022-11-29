library(ggplot2)
library(ggrastr)
library(dplyr)
library(RColorBrewer)
library(scales)

source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_data_fxn.R")

cell_class_expr_distribution <- function(datinfo, na_cts=c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class"), gene_list=NULL, top_n=NULL){
  
  ct_stats <- prep_CT_stats(
    datinfo, 
    expr_type,
    pc_genes=T
  )
  
  datinfo <- merge(datinfo, ct_stats, by="Dataset")

  datinfo$Platform_Class <- "Droplet-based"
  datinfo$Platform_Class[grep("SMART", datinfo$Platform)] <- "Plate-based"
  datinfo$Platform_Class <- factor(
    datinfo$Platform_Class, levels=c("Droplet-based", "Plate-based")
  )
  
  datinfo <- datinfo[!is.element(datinfo$Cell_Class, "VSMC"),]
  datinfo <- datinfo[!grepl(paste(na_cts, collapse="|"), datinfo$Cell_Type),]
  
  n_cts <- n_distinct(paste(datinfo$Study, datinfo$Cell_Type))
  
  plot_title <- "Median # Unique Genes Expressed"
  plot_sub <- paste(n_cts, "cell types from", n_distinct(datinfo$Study), "studies")
  
  pdf(file=paste0("figures/cell_class_mean_expr_", n_cts, "_CTs_", n_distinct(datinfo$Dataset), "_datasets_", n_distinct(datinfo$Study), "_studies.pdf"), width=12, height=9)
  
  idx <- which(
    datinfo$Platform_Class=="Droplet-based"
  )
  
  temp <- datinfo[idx,] %>%
    dplyr::group_by(Cell_Class) %>%
    dplyr::summarise(
      n=median(Median_Unique_Genes)
    ) %>%
    arrange(n)
  
  datinfo$Cell_Class <- factor(datinfo$Cell_Class, levels=temp$Cell_Class)
  
  print(
    ggplot(datinfo, aes(x=Cell_Class, y=Median_Unique_Genes, fill=Cell_Class)) +
      geom_violin(trim=F, size=.7, draw_quantiles=c(.5)) +
      ggrastr::geom_jitter_rast(size=.3, alpha=.25, show.legend=F) +
      theme_minimal() +
      theme(
        plot.title=element_text(hjust=0, face="bold", size=15, lineheight=1.1, margin=margin(b=7)),
        plot.subtitle=element_text(hjust=0, size=12, lineheight=1.1, margin=margin(b=12)),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10, color="black"),
        axis.title.x=element_text(size=12, face="bold", color="black", margin=margin(t=10, b=22)),
        axis.text.y=element_blank(), #element_text(size=12, color="black"),
        panel.border=element_rect(colour="black", fill=NA),
        plot.margin=unit(c(1, 2, 1, 2), "cm"),
        panel.grid.major.y=element_blank(),
      ) +
      ylab("# Unique Genes") +
      labs(title=plot_title, subtitle=plot_sub) +
      scale_y_continuous(labels=comma) +
      scale_color_manual(values=rev(brewer.pal(8, "Set2"))) +
      scale_fill_manual(values=rev(brewer.pal(8, "Set2"))) +
      guides(fill=guide_legend(
        byrow=T, nrow=1, reverse=T
      )) +
      coord_flip() +
      facet_grid(.~Platform_Class, scales="free_x") +
      theme(strip.text=element_text(size=13, face="bold")) 
  )
  
  plot_title <- "Median # UMIs"

  temp <- datinfo[idx,] %>%
    dplyr::group_by(Cell_Class) %>%
    dplyr::summarise(
      n=median(Median_UMIs)
    ) %>%
    arrange(n)

  datinfo$Cell_Class <- factor(datinfo$Cell_Class, levels=temp$Cell_Class)

  print(
    ggplot(datinfo, aes(x=Cell_Class, y=Median_UMIs, fill=Cell_Class)) +
      geom_violin(trim=F, size=.7, draw_quantiles=c(.5)) +
      ggrastr::geom_jitter_rast(size=.35, alpha=.25, show.legend=F) +
      theme_minimal() +
      theme(
        plot.title=element_text(hjust=0, face="bold", size=15, lineheight=1.1, margin=margin(b=7)),
        plot.subtitle=element_text(hjust=0, size=12, lineheight=1.1, margin=margin(b=12)),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10, color="black"),
        axis.title.x=element_text(size=12, face="bold", color="black", margin=margin(t=10, b=22)),
        axis.text.y=element_blank(), #element_text(size=12, color="black"),
        panel.border=element_rect(colour="black", fill=NA),
        plot.margin=unit(c(1, 2, 1, 2), "cm"),
        panel.grid.major.y=element_blank(),
      ) +
      ylab("# UMIs") +
      labs(title=plot_title, subtitle=plot_sub) +
      scale_y_continuous(labels=comma) +
      # scale_color_manual(values=rev(brewer.pal(8, "Set2"))) +
      # scale_fill_manual(values=rev(brewer.pal(8, "Set2"))) +
      guides(fill=guide_legend(
        byrow=T, nrow=1, reverse=T
      )) +
      coord_flip() +
      facet_grid(.~Platform_Class, scales="free_x") +
      theme(strip.text=element_text(size=13, face="bold"))
  )

  print(p + scale_y_continuous(trans="log2", labels=comma))
  
  dev.off()
  
}