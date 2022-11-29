library(reshape2)
library(dplyr)
library(tidygraph)
library(ggraph)
library(viridis)
library(scales)

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_stats_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

MNN_dataset_network <- function(datinfo, df, data_type, expr_type, sim_type){
  
  ## Restrict to comparisons in which there are fewer CTs in dataset1 than dataset2 (not all cells from dataset1 will have the "chance" to be a nearest neighbor in dataset2 if there are fewer cell types in dataset2):
  
  df <- df[df$Dataset1_No.CTs<=df$Dataset2_No.CTs,]
  
  ## Summarize % of MNN per pairwise dataset comparison:
  
  summary_df <- df %>%
    dplyr::group_by(Dataset1, Dataset2) %>%
    dplyr::summarise(
      n=sum(Dataset1_CT==Dataset2_NN)/n()*100
    ) %>%
    arrange(desc(n))
  
  df <- reshape2::dcast(summary_df, Dataset1~Dataset2)
  
}

MNN_network <- function(datinfo, df, data_type, expr_type, sim_type, neu_subtypes){
  
  neu_subtypes <- neu_subtypes %>% na_if("") %>%
    dplyr::select(-c(Cell_Type, Cell_Class, Dataset)) 
  
  ## Get CT stats:
  
  ct_stats <- prep_CT_stats(datinfo, expr_type, pc_genes=T)
  ct_stats <- merge(ct_stats, neu_subtypes, by="Label", all.x=T)

  ## Restrict to cell types with MNN:
  
  df <- df[df$Dataset1_CT==df$Dataset2_NN,]
  
  ## Restrict to comparisons in which there are fewer CTs in dataset1 than dataset2 (not all cells from dataset1 will have the "chance" to be a nearest neighbor in dataset2 if there are fewer cell types in dataset2):
  
  df <- df[df$Dataset1_No.CTs<=df$Dataset2_No.CTs,]
  
  ## Don't include comparisons between identical cell types / cell types from the same study:
  
  temp <- datinfo %>% dplyr::select(Dataset, First_Author)
  
  df <- merge(df, datinfo, by.x="Dataset1", by.y="Dataset")
  df <- merge(df, temp, by.x="Dataset2", by.y="Dataset")
  df <- df[df$First_Author.x!=df$First_Author.y,]
  df <- df[df$Dataset1_CT!=df$Dataset1_NN,]
  
  ## Create unique identifiers:
  
  df$Plot_Label1 <- datinfo$Dataset[match(df$Dataset1, datinfo$Dataset)]
  df$Plot_Label2 <- datinfo$Dataset[match(df$Dataset2, datinfo$Dataset)]
  df$Edge1 <- paste(df$Plot_Label1, df$Dataset1_CT)
  df$Edge2 <- paste(df$Plot_Label2, df$Dataset1_NN)

  ### All cells:
  
  ## Get data into the right format:
  
  edges <- data.frame(from=df$Edge1, to=df$Edge2, Sim=df$Sim)
  net <- tidygraph::as_tbl_graph(edges)
  nodes <- net %>% activate(nodes) %>% data.frame()
  edges <- net %>% activate(edges) %>% data.frame()
  
  nodes$Dataset <- sapply(strsplit(nodes$name, " "), function(x) paste(x[-length(x)], collapse=" "))
  nodes$Study <- datinfo$Study[match(nodes$Dataset, datinfo$Dataset)]
  
  ## Add additional metadata to nodes df:
  
  ids <- make.names(gsub(" ", "_", nodes$name))
  ct_stats$Label <-  make.names(ct_stats$Label)
  ct_stats <- ct_stats[is.element(ct_stats$Label, ids),]
  
  idx <- match(
    ids, ct_stats$Label
  )
  nodes$No.Nuclei <- ct_stats$No.Nuclei[idx]
  nodes$Cell_Type <- ct_stats$Cell_Type[idx]
  nodes$Cell_Class <- ct_stats$Cell_Class[idx]
  nodes$Median_Unique_Genes <- ct_stats$Median_Unique_Genes[idx]
  
  ## Remake graph network:
  
  net <- tbl_graph(nodes=nodes, edges=edges)
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/eigencell_MNN_graph_network_", toupper(sim_type), "_similarity_", data_type, "_", expr_type, "_", n_distinct(c(df$Dataset2, df$Dataset1)), "_datasets.pdf"), width=14, height=11)
  
  plot_title <- paste("Mutual Nearest Neighbors")
  plot_sub <- paste(
    n_distinct(paste(nodes$Cell_Type, nodes$Study)), "cell types from", n_distinct(nodes$Study), "studies"
  )
  
  color_var <- "Study"; color_lab <- NA
  color_var_type <- "discrete"
  network_plot_template(net, plot_title, plot_sub, color_var, color_lab, color_var_type, layout="fr")
  
  color_var <- "Median_Unique_Genes"; color_lab <- "Median # Unique Genes"
  color_var_type <- "continuous"
  network_plot_template(net, plot_title, plot_sub, color_var, color_lab, color_var_type, layout="fr")
  
  
  ### By cell class:
  
  cell_classes <- unique(df$Dataset1_Class)
  
  for(j in 1:length(cell_classes)){
    
    cell_class <- cell_classes[j]
    cell_class_name <- cell_class_full_name(cell_class)
    
    cat("\n")
    print(cell_class)
    
    df1 <- df[is.element(df$Dataset1_Class, cell_class),]
    
    ## Get data into the right format:
    
    edges <- data.frame(from=df1$Edge1, to=df1$Edge2, Sim=df1$Sim)
    net <- tidygraph::as_tbl_graph(edges)
    nodes <- net %>% activate(nodes) %>% data.frame()
    edges <- net %>% activate(edges) %>% data.frame()
    
    nodes$Dataset <- sapply(strsplit(nodes$name, " "), function(x) paste(x[-length(x)], collapse=" "))
    nodes$Study <- datinfo$Study[match(nodes$Dataset, datinfo$Dataset)]
    
    ## Add additional metadata to nodes df:
    
    ids <- make.names(gsub(" ", "_", nodes$name))
    ct_stats1 <- ct_stats[is.element(ct_stats$Label, ids),]
    
    idx <- match(
      ids, ct_stats1$Label
    )
    nodes$No.Nuclei <- ct_stats1$No.Nuclei[idx]
    nodes$Cell_Type <- ct_stats1$Cell_Type[idx]
    nodes$Cell_Class <- ct_stats1$Cell_Class[idx]
    nodes$Class_Level1 <- ct_stats1$Class_Level1[idx]
    nodes$Class_Level2 <- ct_stats1$Class_Level2[idx]
    nodes$Median_Unique_Genes <- ct_stats1$Median_Unique_Genes[idx]
    
    ## Remake graph network:
    
    net <- tbl_graph(nodes=nodes, edges=edges)
    
    plot_title <- paste(cell_class_name, "Mutual Nearest Neighbors")
    plot_sub <- paste(
      n_distinct(paste(nodes$Cell_Type, nodes$Study)), "cell types from", n_distinct(nodes$Study), "studies"
    )    
    
    color_var <- "Study"; color_lab <- NA
    color_var_type <- "discrete"
    network_plot_template(net, plot_title, plot_sub, color_var, color_lab, color_var_type, layout="fr")
    
    color_var <- "Median_Unique_Genes"; color_lab <- "Median # Unique Genes"
    color_var_type <- "continuous"
    network_plot_template(net, plot_title, plot_sub, color_var, color_lab, color_var_type, layout="fr")
    
    if(cell_class=="INH"){
      
      color_var <- "Class_Level2"; color_lab <- NA
      color_var_type <- "discrete"
      network_plot_template(net, plot_title, plot_sub, color_var, color_lab, color_var_type, layout="fr")
      
    }
    
    if(cell_class=="EXC"){
    
      color_var <- "Class_Level1"; color_lab <- NA
      color_var_type <- "discrete"
      network_plot_template(net, plot_title, plot_sub, color_var, color_lab, color_var_type, layout="fr")
      
      color_var <- "Class_Level2"; color_lab <- NA
      color_var_type <- "discrete"
      network_plot_template(net, plot_title, plot_sub, color_var, color_lab, color_var_type, layout="fr")
      
    }
  
  } ## for(j in 1:length(cell_classes)){

  dev.off()
  
}

network_plot_template <- function(net, plot_title, plot_sub, color_var, color_lab, color_var_type=c("discrete", "continuous"), layout=c("fr", "stress")){
  
  p <- ggraph(net, layout=layout) + 
    geom_edge_link(aes(alpha=Sim), position=position_jitter(), show.legend=F) + 
    scale_edge_alpha(range=c(0, .4)) +
    geom_node_point(aes(fill=!!sym(color_var)), position=position_jitter(), size=3.5, shape=21) +
    theme_classic() +
    theme(
      plot.title=element_text(hjust=0, face="bold", size=15, margin=margin(b=5)),
      plot.subtitle=element_text(hjust=0, size=12, lineheight=1.2, margin=margin(b=5)),
      legend.text=element_text(size=12),
      legend.position="bottom",
      legend.direction="horizontal",
      legend.box="vertical",
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      panel.border=element_rect(colour="black", fill=NA),
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    labs(title=plot_title, subtitle=plot_sub) +
    guides(fill=guide_legend(override.aes=list(size=5.5)))
  
  if(color_var_type=="discrete"){
    
    p <- p + scale_fill_manual(values=brewer_fxn(n_distinct(nodes[,color_var]))) +
      theme(legend.title=element_blank()) 
    
  } else { ## if(color_var_type=="discrete"){
    
    p <- p + theme(
      legend.title=element_text(size=12, margin=margin(b=5, t=5))
    ) +
      guides(fill=guide_colourbar(
        barwidth=20, title.position="top", title.hjust=.5, 
        frame.colour="black", frame.linewidth=.3, 
        label.theme=element_text(size=9, hjust=.8, vjust=.9)
      )) +
      scale_fill_continuous(name=color_lab, type="viridis", labels=comma) 
    
  }
  
  print(p)
  
} ## network_plot_template
