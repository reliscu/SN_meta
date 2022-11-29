library(GGally)
library(viridis)
library(RColorBrewer)
library(data.table)
library(plyr)
library(dplyr)
library(propr)

source("/home/rebecca/code/misc/upper_first.R")
source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_data_fxn.R")
source("/home/rebecca/SCSN_meta_analysis/code/cell_class_full_name_fxn.R")

palette <- readRDS("/home/rebecca/SCSN_meta_analysis/palette_15.RDS")

model_cluster_misfits <- function(
  ct_expr, 
  datinfo, 
  data_type=c("author_data"), 
  expr_type=c("counts", "normalized_counts"), 
  which_genes=c("intersection", "union"), 
  sim_type, 
  prop_metric=c("rho", "phs"), 
  clust_method, 
  top_n=NULL, 
  n_clusters
){
  
  ct_data <- prep_CT_data(
    datinfo, 
    ct_df=ct_expr, 
    expr_type, 
    which_genes, 
    top_n, 
    ct_expr=NULL
  )
  
  datinfo <- ct_data[[1]]
  ct_expr <- ct_data[[2]]
  
  if(is.null(top_n)){
    
    n_genes <- nrow(ct_expr)
    
  } else {
    
    n_genes <- paste("top", nrow(ct_expr), sep="_")
    
  }
  
  n_cts <- length(unique(paste(datinfo$Study, datinfo$Cell_Type)))
  
  if(sim_type=="prop"){
    
    sim <- propr(ct_expr[,-c(1)], metric=prop_metric, symmetrize=T, p=1)@matrix
    sim_type <- paste(sim_type, prop_metric, sep="_")
    if(grepl("ph", prop_metric)){ # phs and phi are already distance measures
      sim <- 1-sim
    }
    
  } else {
    sim <- cor(ct_expr[,-c(1)], method=sim_type, use="pairwise.complete.obs")
  }
  
  clust <- stats::hclust(
    as.dist(1-sim), method=clust_method
  )
  clust_labs <- stats::cutree(clust, k=n_clusters)
  
  misfits <- c()
  
  for(i in 1:n_clusters){
    
    working_clust <- clust_labs[clust_labs==i]
    datinfo_clust <- datinfo[match(names(working_clust), datinfo$Label),]
    
    if(!identical(names(working_clust), datinfo_clust$Label)){
      stop("!identical(names(working_clust), datinfo_clust$Label)")
    }
    
    n_cts_per_class <- table(datinfo_clust$Cell_Class)
    clust_class <- names(n_cts_per_class[which.max(n_cts_per_class)])
    clust_misfits <- names(working_clust)[datinfo_clust$Cell_Class!=clust_class]
    
    if(length(clust_misfits)>0){
      misfits <- c(misfits, clust_misfits)
    }
    
  } ## for(i in 1:n_clusters){
  
  datinfo$Misfit_Status <- F
  datinfo$Misfit_Status[is.element(datinfo$Label, misfits)] <- T
  
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/mean_expr_cluster_misfit_status_vs_covariates_", data_type, "_", expr_type, "_", toupper(sim_type), "_", toupper(clust_method), "_", n_clusters, "_clusters_", n_cts, "_CTs_", n_distinct(datinfo$Dataset), "_datasets_", nrow(ct_expr), "_intersection_genes.pdf"), width=12, height=10)
  
  print(head(datinfo))
  
  color_var <- "Cell_Class"; color_lab <- "Cell Class"
  color_var_type <- "discrete"
  misfit_model_plot_template(df=datinfo, color_var, color_lab, color_var_type, n_cts)
  
  color_var <- "Study"; color_lab <- "Study"
  color_var_type <- "discrete"
  misfit_model_plot_template(df=datinfo, color_var, color_lab, color_var_type, n_cts)
  
  color_var <- "Median_Unique_Genes"; color_lab <- "Median # Unique Genes per Nucleus"
  color_var_type <- "continuous"
  misfit_model_plot_template(df=datinfo, color_var, color_lab, color_var_type, n_cts)
  
  color_var <- "Median_UMIs"; color_lab <- "Median # UMIs per Nucleus"
  color_var_type <- "continuous"
  misfit_model_plot_template(df=datinfo, color_var, color_lab, color_var_type, n_cts)
  
  color_var <- "No.Nuclei"; color_lab <- "# Nuclei"
  color_var_type <- "continuous"
  misfit_model_plot_template(df=datinfo, color_var, color_lab, color_var_type, n_cts)
  
  dev.off()
  
} ## model_cluster_misfits <- function(
  
misfit_model_plot_template <- function(
  df, color_var, color_lab, color_var_type=c("continous", "discrete"), n_cts
){
  
  r2 <- summary(lm(df$Misfit_Status~df[,color_var]))$adj.r.squared
  
  plot_title <- paste("Mean Expression Cluster Misfit Status vs.", color_lab)
  plot_sub <- paste0(
    "adj. R2 = ", signif(r2, 2), "\n", upper_first(gsub("_", " ", sim_type)), " Similarity, ", upper_first(clust_method), " Agglomeration", "\n", n_clusters, " Clusters, ", n_cts, " Cell Types"
  )
  
  if(color_var_type=="discrete"){
    
    p <-  
      ggplot(df, aes_string(x=color_var, y="Misfit_Status", color=color_var)) +
      geom_jitter(width=.3, size=.5) +
      theme_classic() + 
      theme(plot.title=element_text(hjust=.5, face="bold", size=15),
            plot.subtitle=element_text(hjust=.5, size=13, lineheight=1, margin=margin(t=4, b=8)),
            legend.title=element_text(size=12, face="bold"),
            legend.position="bottom",
            legend.box="vertical",
            axis.title.x=element_blank(), 
            axis.text.x=element_blank(),
            axis.ticks.x=element_line(size=0),
            axis.title.y=element_text(size=11, margin=margin(r=8)),
            axis.text.y=element_text(face="bold"),
            plot.margin = unit(c(2, 2, 2, 2), "cm")) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab("") + ylab("Cluster Misfit Status") + 
      guides(color=guide_legend(title=color_lab, override.aes=list(alpha=.8, size=3))) +
      scale_y_discrete(labels=c("F", "T"))
    
    if(n_distinct(df[,color_var])>15){
      p <- p + scale_color_manual(
        values=colorRampPalette(brewer.pal(9, "Set1"))(n_distinct(df[,color_var]))
      )
    } else {
      p <- p + scale_color_manual(values=palette[1:n_distinct(df[,color_var])])
    }
    
  } else {
    
    min_val <- min(df[,color_var], na.rm=T); max_val <- max(df[,color_var], na.rm=T)
    color_breaks <- c(
      round_any(quantile(df[,color_var], probs=seq(0, .5, by=.1)), 10),
      round_any(quantile(df[,color_var], probs=seq(.6, 1, by=.1)), 100)
    )
    color_breaks[color_breaks>1e6] <- max_val
    color_breaks[which.min(color_breaks)] <- min_val
    color_breaks <- unique(color_breaks)
    color_labs <- sapply(as.numeric(color_breaks), function(x) format(x, big.mark=","))
    
    p <- 
      ggplot(df, aes_string(x=color_var, y="Misfit_Status", color=color_var)) +
      geom_jitter(width=.4, size=.5) + 
      theme_classic() + 
      theme(
        plot.title=element_text(hjust=.5, face="bold", size=15),
        plot.subtitle=element_text(hjust=.5, size=13, lineheight=1, margin=margin(t=4, b=8)),
        legend.title=element_text(size=12),
        legend.position="bottom",
        legend.box="vertical",
        axis.title.y=element_text(size=11, margin=margin(r=8)),
        plot.margin = unit(c(2, 2, 2, 2), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + 
      xlab("") + ylab("Cluster Misfit Status") + 
      scale_color_viridis(
        guide=guide_colourbar(
          title=color_lab, 
          barwidth=25, 
          label.theme=element_text(
            size=8, angle=40, hjust=.8
          ),
          title.position="top",
          title.hjust=.54
        ),
        trans="log", 
        limits=c(min_val, max_val), 
        labels=color_labs, 
        breaks=color_breaks,
      ) +
      scale_y_discrete(labels=c("F", "T"))
    
  }
  
  print(p)
  
} ## misfit_model_plot_template <- function(

plot_cluster_misfits <- function(
  ct_expr, 
  datinfo, 
  data_type=c("author_data"), 
  expr_type=c("counts", "normalized_counts"), 
  which_genes=c("intersection", "union"), 
  sim_type, 
  prop_metric=c("rho", "phs"), 
  clust_method, 
  top_n=NULL, 
  n_clusters
){
  
  ct_data <- prep_CT_data(
    datinfo, 
    ct_df=ct_expr, 
    expr_type, 
    which_genes, 
    top_n, 
    ct_expr=NULL
  )
  
  datinfo <- ct_data[[1]]
  ct_expr <- ct_data[[2]]
  
  if(is.null(top_n)){
    
    n_genes <- nrow(ct_expr)
    
  } else {
    
    n_genes <- paste("top", nrow(ct_expr), sep="_")
    
  }
  
  n_cts <- length(unique(paste(datinfo$Study, datinfo$Cell_Type)))

  if(sim_type=="prop"){
    
    sim <- propr(ct_expr[,-c(1)], metric=prop_metric, symmetrize=T, p=1)@matrix
    sim_type <- paste(sim_type, prop_metric, sep="_")
    if(grepl("ph", prop_metric)){ # phs and phi are already distance measures
      sim <- 1-sim
    }
    
  } else {
    sim <- cor(ct_expr[,-c(1)], method=sim_type, use="pairwise.complete.obs")
  }
  
  clust <- stats::hclust(
    as.dist(1-sim), method=clust_method
  )
  clust_labs <- stats::cutree(clust, k=n_clusters)
  
  names(clust_labs) <- make.names(names(clust_labs))
  colnames(sim) <- rownames(sim) <- make.names(colnames(sim))
  colnames(ct_expr) <- make.names(colnames(ct_expr))
  
  ## Normalize mean expression values per cell type so values are comparable for plotting:
  
  ct_expr[,-c(1)] <- scale(ct_expr[,-c(1)], center=FALSE, scale=colSums(ct_expr[,-c(1)]))  
  
  for(i in 1:n_clusters){
    
    working_clust <- clust_labs[clust_labs==i]
    datinfo_clust <- datinfo[match(names(working_clust), datinfo$Label),]
    
    if(!identical(names(working_clust), datinfo_clust$Label)){
      stop("!identical(names(working_clust), datinfo_clust$Label)")
    }
  
    n_cts_per_class <- table(datinfo_clust$Cell_Class)
    clust_class <- names(n_cts_per_class[which.max(n_cts_per_class)])
    clust_class_name <- cell_class_full_name(clust_class)
    clust_misfits <- names(working_clust)[datinfo_clust$Cell_Class!=clust_class]

    if(length(clust_misfits)==0){
      print(paste(clust_class_name, "cluster", i, "has no misfits"))
      next
    }
    
    clust_idx <- match(names(working_clust), colnames(sim))
    clust_sim <- sim[clust_idx, clust_idx]
    n_neighbors <- 5
  
    pdf(file=paste0("figures/", data_type, "/", expr_type, "/", clust_class, "_mean_expr_cluster_no.", i, "_misfits_scatterplot_", data_type, "_", expr_type, "_", toupper(sim_type), "_", toupper(clust_method), "_",  n_clusters, "_clusters_", n_cts, "_CTs_", n_distinct(datinfo$Dataset), "_datasets_", nrow(ct_expr), "_intersection_genes.pdf"), width=12, height=10)
    
    for(j in 1:length(clust_misfits)){
      
     working_ct <- clust_misfits[j]
      datinfo_ct <- 
        datinfo_clust[match(working_ct, datinfo_clust$Label),]
      working_ct_class <- 
        datinfo$Cell_Class[is.element(datinfo$Label, working_ct)]
      working_ct_class_name <- 
        cell_class_full_name(working_ct_class)
      working_ct_name <- 
        paste0(
          gsub(".", "-", gsub("_", " ", datinfo_ct$Dataset), fixed=T), "\n",
          paste(working_ct_class, gsub(".", "-", gsub("_", " ", datinfo_ct$Cell_Type), fixed=T))
        )
      
      ## For each misfit cluster member: 
      
      ## 1) Compare its expression to its nearest cluster neighbors:
      
      ct_sim <- clust_sim[,match(working_ct, colnames(clust_sim))]
      clust_neighbor_sim <- 
        sort(ct_sim, decreasing=T)[2:(n_neighbors+1)]
      clust_neighbor_cts <- names(clust_neighbor_sim)
      clust_neighbor_fid <- reshape2::melt(
        data.frame(ct_expr[,match(clust_neighbor_cts, colnames(ct_expr))])
      )
      datinfo_clust_neighbors <- datinfo_clust[match(c(clust_neighbor_cts), datinfo_clust$Label),]
      df <- data.frame(
        rep(ct_expr[,match(working_ct, colnames(ct_expr))], n_neighbors), ## Working cell type mean expression
        clust_neighbor_fid
      )
      colnames(df) <- c(working_ct, "Cell_Type", "Mean_Expr")
      df$Dataset <- 
        datinfo_clust_neighbors$Dataset[match(df$Cell_Type, datinfo_clust_neighbors$Label)]
      
      plot_title <- paste0(
        gsub(working_ct_class, working_ct_class_name, working_ct_name), 
        "\nvs.\nNearest Mean Expression Cluster Neighbors"
      )
      plot_sub <- paste0(
        clust_class_name, " Cluster\n",
        upper_first(gsub("_", " ", sim_type)), " Similarity, ", upper_first(clust_method), " Agglomeration", "\n",
        format(nrow(ct_expr), big.mark=","), " genes"
      )
      x_lab <- 
        paste(
          gsub(working_ct_class, working_ct_class_name, working_ct_name), "Mean Expression"
        )
      y_lab <- paste("Nearest Cluster Neighbor Mean Expression")
      leg_labs <- gsub("_", " ", datinfo_clust_neighbors$Dataset)
      facet_labs <- paste0(datinfo_clust_neighbors$Cell_Class, "\n", gsub("_", " ", datinfo_clust_neighbors$Cell_Type))
      facet_labs <- gsub(".", "-", facet_labs, fixed=T)
      names(facet_labs) <- datinfo_clust_neighbors$Label
      sim_labs <- data.frame(
        Cell_Type=clust_neighbor_cts,
        Sim=paste("Sim:\n", signif(clust_neighbor_sim, 2)),
        Dataset=datinfo_clust_neighbors$Dataset[match(clust_neighbor_cts, datinfo_clust_neighbors$Label)]
      )
      
      scatterplot_template(
        df, working_ct, 
        n_neighbors, 
        plot_title, 
        plot_sub, 
        sim_labs,
        x_lab, y_lab, 
        leg_labs, 
        facet_labs
      )
      
      ## 2) Compare its nearest neighbors dataset-wide
      
      working_ct_sim <- sim[,which(is.element(colnames(sim), working_ct))]
      ct_neighbor_sim <- 
        sort(working_ct_sim, decreasing=T)[2:(n_neighbors+1)]
      ct_neighbor_cts <- names(ct_neighbor_sim)
      ct_neighbor_fid <- reshape2::melt(
        data.frame(ct_expr[,match(ct_neighbor_cts, colnames(ct_expr))])
      )
      datinfo_ct_neighbors <- datinfo[match(c(ct_neighbor_cts), datinfo$Label),]
      df <- data.frame(
        rep(ct_expr[,match(working_ct, colnames(ct_expr))], n_neighbors), ## Working cell type mean expression
        ct_neighbor_fid
      )
      colnames(df) <- c(working_ct, "Cell_Type", "Mean_Expr")
      df$Dataset <- 
        datinfo_ct_neighbors$Dataset[match(df$Cell_Type, datinfo_ct_neighbors$Label)]
      
      plot_title <- paste0(
        gsub(working_ct_class, working_ct_class_name, working_ct_name), 
        "\nvs.\nMNearest Mean Expression Neighbors"
      )
      plot_sub <- paste0(
        upper_first(gsub("_", " ", sim_type)), " Similarity\n", 
        upper_first(clust_method), " Agglomeration", "\n",
        format(nrow(ct_expr), big.mark=","), " genes"
      )
      y_lab <- paste("Nearst Neighbor Mean Expression")
      leg_labs <- gsub("_", " ", datinfo_ct_neighbors$Dataset)
      facet_labs <- paste0(datinfo_ct_neighbors$Cell_Class, "\n", gsub("_", " ", datinfo_ct_neighbors$Cell_Type))
      facet_labs <- gsub(".", "-", facet_labs, fixed=T)
      names(facet_labs) <- datinfo_ct_neighbors$Label
      sim_labs <- data.frame(
        Cell_Type=ct_neighbor_cts,
        Sim=paste("Sim:\n", signif(ct_neighbor_sim, 2)),
        Dataset=datinfo_ct_neighbors$Dataset[match(ct_neighbor_cts, datinfo_ct_neighbors$Label)]
      )
      
      scatterplot_template(
        df, working_ct, 
        n_neighbors, 
        plot_title, 
        plot_sub, 
        sim_labs,
        x_lab, y_lab, 
        leg_labs, 
        facet_labs
      )
      
      ## 2) Compare its expression to canonical cell types of the same class as "misfit" cell type

      class_idx <- 
        match(datinfo$Label[is.element(datinfo$Cell_Class, working_ct_class)], colnames(sim))
      class_sim <-  sim[class_idx, class_idx]
      class_sim <- (class_sim+1)/2
      connectivity <- sort(colSums(class_sim), decreasing=T)
      canonical_cts <- names(connectivity)[1:n_neighbors]
      canon_sim <- sim[match(canonical_cts, colnames(sim)), match(working_ct, colnames(sim))]
      canon_fid <- reshape2::melt(data.frame(
        ct_expr[,match(canonical_cts, colnames(ct_expr))]
      ))
      datinfo_canon <- 
        datinfo[match(canonical_cts, datinfo$Label),]
      df <- data.frame(
        rep(ct_expr[,match(working_ct, colnames(ct_expr))], n_neighbors),
        canon_fid
      )
      colnames(df) <- c(working_ct, "Cell_Type", "Mean_Expr")
      df$Dataset <- 
        datinfo_canon$Dataset[match(df$Cell_Type, datinfo_canon$Label)]
      
      plot_title <- paste0(
        gsub(working_ct_class, working_ct_class_name, working_ct_name), 
        "\nvs.\nMean Expression Canonical ", working_ct_class_name, " Cell Types"
      )
      plot_sub <- paste0(
        upper_first(gsub("_", " ", sim_type)), " Similarity\n", 
        upper_first(clust_method), " Agglomeration", "\n",
        format(nrow(ct_expr), big.mark=","), " genes"
      )
      x_lab <- paste(working_ct_name, "Mean Expression")
      y_lab <- paste("Canonical", working_ct_class, "Mean Expression")
      leg_labs <- gsub("_", " ", datinfo_canon$Dataset)
      facet_labs <- paste0(datinfo_canon$Cell_Class, "\n", gsub("_", " ", datinfo_canon$Cell_Type))
      facet_labs <- gsub(".", "-", facet_labs, fixed=T)
      names(facet_labs) <- datinfo_canon$Label
      sim_labs <- data.frame(
        Cell_Type=canonical_cts,
        Sim=paste("Sim:\n", signif(canon_sim, 2)),
        Dataset=datinfo_canon$Dataset[match(canonical_cts, datinfo_canon$Label)]
      )

      scatterplot_template(
        df, working_ct, 
        n_neighbors, 
        plot_title, 
        plot_sub, 
        sim_labs,
        x_lab, y_lab, 
        leg_labs, 
        facet_labs
      )

    } ## for(j in 1:length(clust_nonclass)){
    
    dev.off()
    
  } ## for(i in 1:n_clusters){
    
} ## plot_cluster_misfits <- function(

scatterplot_template <- function(
  df, working_ct, n_neighbors, 
  plot_title, plot_sub, sim_labs, 
  x_lab, y_lab, leg_labs, facet_labs
){
  
  p <- 
    ggplot(df, aes_string(x=working_ct, y="Mean_Expr", color="Dataset")) +
    geom_point(size=1) +
    geom_text(
      data=sim_labs, 
      aes(label=Sim, x=Inf, y=Inf, vjust=1.5, hjust=2), 
      size=3.5, fontface="bold"
    ) +
    theme_classic() + 
    theme(plot.title=element_text(hjust=.5, face="bold", size=15),
          plot.subtitle=element_text(hjust=.5, size=13, lineheight=1.1, margin=margin(t=6, b=12)),
          legend.title=element_text(size=12, face="bold"),
          legend.text=element_text(size=12),
          legend.position="bottom",
          legend.box="vertical",
          axis.ticks.x=element_line(size=0),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=13, face="bold", margin=margin(t=9)), 
          axis.text.y=element_blank(),
          axis.title.y=element_text(size=13, face="bold", margin=margin(r=9)),
          axis.ticks.y=element_line(size=0),
          plot.margin = unit(c(2, 2, 2, 2), "cm")) +
    labs(title=plot_title, subtitle=plot_sub) + 
    xlab(x_lab) +
    ylab(y_lab) +
    guides(color=guide_legend(
      title="Dataset", ncol=2, 
      title.position="top", 
      title.hjust=.5, 
      override.aes=list(size=4))
    ) +
    facet_grid(.~Cell_Type, labeller=labeller(Cell_Type=facet_labs)) +
    theme(
      strip.text.x=element_text(size=8, face="bold"),
      strip.background.x=element_rect(fill="white")
    ) +
    scale_color_manual(
      values=palette[1:n_neighbors],
      labels=sort(unique(leg_labs))
    )
    
  print(p)

} ## scatterplot_template <- function(

plot_highly_connected_class_members <- function(
  ct_expr, 
  datinfo, 
  data_type=c("author_data"), 
  expr_type=c("counts", "normalized_counts"), 
  cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"),
  which_genes=c("intersection", "union"), 
  sim_type, 
  prop_metric=c("rho", "phs"), 
  top_n=NULL
){
  
  ct_data <- prep_CT_data(
    datinfo, 
    ct_df=ct_expr, 
    expr_type, 
    which_genes, 
    top_n, 
    ct_expr=NULL
  )
  
  datinfo <- ct_data[[1]]
  ct_expr <- ct_data[[2]]
  
  if(is.null(top_n)){
    
    n_genes <- nrow(ct_expr)
    
  } else {
    
    n_genes <- paste("top", nrow(ct_expr), sep="_")
    
  }
  
  if(sim_type=="prop"){
    
    sim <- propr(ct_expr[,-c(1)], metric=prop_metric, symmetrize=T, p=1)@matrix
    sim_type <- paste(sim_type, prop_metric, sep="_")
    if(grepl("ph", prop_metric)){ # phs and phi are already distance measures
      sim <- 1-sim
    }
    
  } else {
    sim <- cor(ct_expr[,-c(1)], method=sim_type, use="pairwise.complete.obs")
  }
  
  colnames(ct_expr) <- make.names(colnames(ct_expr))
  
  ## Convert sim to adjacency matrix:
  sim <- (sim+1)/2
  
  ## Column normalize mean expression so values are comparable for plotting:
  
  ct_expr[,-c(1)] <- scale(ct_expr[,-c(1)], center=FALSE, scale=colSums(ct_expr[,-c(1)]))  
 
  pdf(file=paste0("figures/", data_type, "/", expr_type, "/mean_expr_highly_connected_CTs_per_cell_class_scatterplot_", data_type, "_", expr_type, "_", toupper(sim_type), "_", n_distinct(datinfo$Dataset), "_datasets_", nrow(ct_expr), "_intersection_genes.pdf"), width=10, height=11)
  
  for(i in 1:length(cell_classes)){
    
    cell_class <- cell_classes[i]
    class_name <- cell_class_full_name(cell_class)
    
    datinfo_class <- datinfo[is.element(datinfo$Cell_Class, cell_class),]

    class_idx <- match(datinfo$Label, colnames(sim))
    class_sim <- sim[class_idx, class_idx]
    
    ## Identifify most highly connected class cell types:
    
    n_neighbors <- 5
    connectivity <- sort(rowSums(class_sim), decreasing=T)
    canonical_cts <- names(connectivity)[1:n_neighbors]
    canon_fid <- ct_expr[,match(canonical_cts, colnames(ct_expr))]
    datinfo_canon <- datinfo[match(c(canonical_cts), datinfo$Label),]
    
    plot_title <- paste("Canonical", class_name, "Mean Expression")
    x_lab <- gsub(".", "-", gsub("_", " ", paste(datinfo_canon$First_Author, datinfo_canon$Year, datinfo_canon$Region_Code, datinfo_canon$Cell_Type)), fixed=T)
    colnames(canon_fid) <- make.unique(x_lab, sep="-")
    
    print(
      ggpairs(canon_fid, title=plot_title, axisLabels="none", labeller=label_wrap_gen(16)) +
        theme_minimal() +
        theme(
          plot.title=element_text(hjust=.5, face="bold")
        )
    )
    
  } # for(i in 1:n_clusters){
  
  dev.off()
  
} ## plot_highly_connected_cluster_members <- function(
