library(ggplot2)
library(plyr) ## round_any()
library(RColorBrewer) ## brewer.pal()
library(scales) ## comma()

source("/home/rebecca/SCSN_meta_analysis/code/brewer_fxn.R")

PC_scatter_template <- function(
  df, x, y, 
  plot_title, 
  plot_sub, 
  color_var, color_lab, 
  color_var_type=c("continuous", "discrete"), 
  color_leg_rows=NULL, ## Only for 'discrete' color_var_type
  color_leg_breaks=NULL, ## Only for 'continuous' color_var_type
  color_leg_labs=NULL, ## Only for 'continuous' color_var_type
  size_var, 
  size_lab, 
  size_leg_breaks, 
  size_leg_labs,
  size_range=c(.5, 7),
  label_outliers=F
){
  
  n_row <- 1
  if(n_distinct(df[,color_var])>5){
    n_row <- 2
  }
  
  y_min <- min(df[,y])-abs(min(df[,y])*.05)
  y_max <- max(df[,y])+abs(max(df[,y])*.05)

  p <- ggplot(df, aes(!!sym(x), !!sym(y), fill=!!sym(color_var), size=!!sym(size_var))) +
    geom_point(alpha=1, shape=21) + 
    theme_classic() + 
    theme(
      plot.title=element_text(hjust=.5, face="bold"), 
      plot.subtitle=element_text(hjust=.5, lineheight=1, margin=margin(t=2, b=10)), 
      legend.title=element_text(size=12),
      legend.text=element_text(size=10),
      legend.position="bottom", 
      legend.box="vertical",
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    ylim(y_min, y_max) +
    labs(title=plot_title, subtitle=plot_sub, color=color_lab, size=size_lab) +
    scale_size_continuous(
      name=gsub("per Cell Type", "", size_lab), 
      range=size_range, 
      labels=size_leg_labs, 
      breaks=size_leg_breaks
    ) +
    guides(size=guide_legend(
      title.position="top", title.hjust=.5, order=1, nrow=1, override.aes=list(fill="black")
    ))
  
  if(color_var_type=="discrete"){
    
    p <- p + guides(
      fill=guide_legend(order=2, nrow=color_leg_rows, override.aes=list(size=5, alpha=1, color="black"))
    )
    
    if(color_var=="Year"){
      
      p <- p + scale_fill_manual(values=colorRampPalette(brewer.pal(9, "YlOrRd"))(n_distinct(df[,color_var])))
      
    } else if(n_distinct(df[,color_var])==2){
      
      p <- p + scale_fill_manual(name=color_lab, values=c("#3492EA", "#D92F4B"))
      
    } else { ## n_distinct(df[,color_var])>=16
      
      p <- p + scale_fill_manual(name=color_lab, values=brewer_fxn(n_distinct(df[,color_var])))
    }
    
  } else { ## if(color_var_type=="discrete"){

    min_val <- min(df[,color_var], na.rm=T); max_val <- max(df[,color_var], na.rm=T)
    
    p <- p + guides(fill=guide_colourbar(barwidth=18, title.position="top", title.hjust=.5, frame.colour="black", frame.linewidth=1, label.theme=element_text(size=9, angle=45, hjust=.8, vjust=.9))) + 
      scale_fill_continuous(name=gsub("per Cell Type", "", color_lab), trans="log", type="viridis", limits=c(min_val, max_val), labels=color_leg_labs, breaks=color_leg_breaks) 
  
    
  } ## if(color_var_type=="discrete"){} else {
  
  if(grepl("Cell", color_var)){
    
    p <- p + geom_point(data=subset(df, is.element(df[,color_var], c("END", "PER"))), alpha=1, shape=21, show.legend=F) +
      geom_point(data=subset(df, is.element(df[,color_var], c("VSMC"))), alpha=1, shape=21, show.legend=F)
    
  } else if(grepl("Platform", color_var)){
    
    p <- p + geom_point(data=subset(df, is.element(df[,color_var], c("DroNc-seq"))), alpha=1, shape=21, show.legend=F)
    
  } else if(grepl("Region", color_var)){
    
    p <- p + geom_point(data=subset(df, is.element(df[,color_var], c("PC", "OC"))), alpha=1, shape=21, show.legend=F)
    
  }
  
  if(label_outliers){
    p <- p + geom_text_repel(aes(label=Outlier_Label), size=3.5, show.legend=F, max.overlaps=1000, seed=666, box.padding=.3)
  }
  
  print(p)
  
}


PC_scatter_template_v2 <- function(
  df, x, y, 
  plot_title, 
  plot_sub, 
  color_var, color_lab, 
  color_var_type=c("continuous", "discrete"), 
  color_leg_rows=NULL, ## Only for 'discrete' color_var_type
  color_leg_breaks=NULL, ## Only for 'continuous' color_var_type
  color_leg_labs=NULL, ## Only for 'continuous' color_var_type
  size_var, 
  size_lab, 
  size_leg_breaks, 
  size_leg_labs,
  size_range=c(.5, 7),
  label_outliers=F
){
  
  y_min <- min(df[,y])-abs(min(df[,y])*.05)
  y_max <- max(df[,y])+abs(max(df[,y])*.05)
  
  p <- ggplot(df, aes(!!sym(x), !!sym(y), fill=!!sym(color_var), size=!!sym(size_var))) +
    geom_point(alpha=.6) + 
    theme_classic() + 
    theme(
      plot.title=element_text(hjust=.5, face="bold"), 
      plot.subtitle=element_text(hjust=.5, lineheight=1, margin=margin(t=2, b=10)), 
      legend.title=element_text(size=12),
      legend.text=element_text(size=10),
      legend.position="bottom", legend.box="vertical",
      plot.margin=unit(c(2, 2, 2, 2), "cm")
    ) +
    ylim(y_min, y_max) +
    labs(title=plot_title, subtitle=plot_sub, color=color_lab, size=size_lab) +
    scale_size_continuous(name=gsub("per Cell Type", "", size_lab), range=size_range, labels=size_leg_labs, breaks=size_leg_breaks) +
    guides(size=guide_legend(title.position="top", title.hjust=.5, order=1, nrow=1))
  
  if(color_var_type=="discrete"){
    
    p <- p + guides(
      color=guide_legend(order=2, nrow=color_leg_rows, override.aes=list(size=5, alpha=1))
    )
    
    if(color_var=="Year"){
      
      p <- p + scale_color_manual(
        values=colorRampPalette(brewer.pal(9, "YlOrRd"))(n_distinct(df[,color_var]))
      )
      
    } else if(n_distinct(df[,color_var])==2){
      
      p <- p + scale_color_manual(name=color_lab, values=c("#3492EA", "#D92F4B"))
      
    } else { ## n_distinct(df[,color_var])>=16
      
      p <- p + scale_color_manual(name=color_lab, values=brewer_fxn(n_distinct(df[,color_var])))
    }
    
  } else { ## if(color_var_type=="discrete"){
    
    min_val <- min(df[,color_var], na.rm=T); max_val <- max(df[,color_var], na.rm=T)
    
    p <- p + guides(color=guide_colourbar(barwidth=18, title.position="top", title.hjust=.5, label.theme=element_text(size=9, angle=45, hjust=.8))) + 
      scale_color_continuous(name=gsub("per Cell Type", "", color_lab), trans="log", type="viridis", limits=c(min_val, max_val), labels=color_leg_labs, breaks=color_leg_breaks) 
    
    
  } ## if(color_var_type=="discrete"){} else {
  
  if(label_outliers){
    p <- p + geom_text_repel(aes(label=Outlier_Label), size=3.5, show.legend=F, max.overlaps=1000, seed=666, box.padding=.3)
  }
  
  print(p)
  
}

