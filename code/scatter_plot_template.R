library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(scales)

scatter_plot_template <- function(df, plot_title, plot_sub, x_var, x_lab, y_var, y_lab){ # x_coord, y_coord
  
  plot_title <- paste0(x_lab, " vs.\n", y_lab)
  if(nchar(x_lab)<20){
    plot_title <- paste(x_lab, "vs.", y_lab)
  }
  
  print(
    ggplot(df, aes(x=!!sym(x_var), y=!!sym(y_var), color=Study)) +
      geom_point(show.legend=F) +
      geom_smooth(method='lm', formula= y~x, color="black", alpha=.2, show.legend=F) +
      geom_text_repel(aes(label=Plot_Label), box.padding=.4, show.legend=F, max.overlaps=100) +
      theme_minimal() +
      theme(
        plot.title=element_text(hjust=0, face="bold", size=15, lineheight=1.1, margin=margin(b=7)),
        plot.subtitle=element_text(hjust=0, size=12, lineheight=1.1, margin=margin(b=6)),
        legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=12, color="black", margin=margin(t=18)),
        axis.text.x=element_text(size=12, color="black"),
        axis.title.y=element_text(size=12, color="black", margin=margin(r=15)),
        axis.text.y=element_text(size=9, color="black"),
        panel.border=element_rect(colour="black", fill=NA),
        plot.margin=unit(c(1, 2, 1, 2), "cm")
      ) +
      labs(title=plot_title, subtitle=plot_sub) + # subtitle=plot_sub
      xlab(x_lab) + ylab(y_lab) + 
      scale_x_continuous(labels=comma) +
      scale_color_manual(values=brewer_fxn(n_distinct(df$Study))) +
      coord_cartesian(expand=F, clip="off")
  )
  
}

# +
#   geom_label(
#     aes(label=r2, x=x_coord, y=y_coord, fontface="italic"),
#     label.padding=unit(.3, "lines"),
#     color="black", size=5, parse=T, label.size=.5
#   ) +