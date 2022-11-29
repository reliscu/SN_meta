library(dplyr)
library(ggplot2)

# CT_real_vs_rand_distro <- function(){
#   
#   datinfo <- datinfo %>% dplyr::select(Dataset, Plot_Label, Author_No.Cell_Types, Author_Median_Unique_Genes_PC)
#   
#   for(i in 1:length(real_sim)){
#     
#     i=8
#     real_sim1 <- real_sim[[i]]
#     rand_sim1 <- rand_sim[[i]]
#     
#     real_sim1$Distro <- "Cell Type"
#     rand_sim1$Distro <- "Random"
#     
#     df <- rbind(real_sim1, rand_sim1)
#     df$ABI1 <- F
#     df$ABI1[is.element(df$Study1, c("ABI 2019", "Bakken 2019", "Hodge 2018", "Luo 2019"))] <- T
#     df$ABI2 <- F
#     df$ABI2[is.element(df$Study2, c("ABI 2019", "Bakken 2019", "Hodge 2018", "Luo 2019"))] <- T
#     df <- df %>% dplyr::filter(ABI1!=ABI2) %>% dplyr::select(-c(ABI1, ABI2))
#     
#     df_temp <- df %>%
#       dplyr::group_by(Dataset1) %>%
#       dplyr::summarise(
#         Percent=sum(
#           Similarity[Distro=="Cell Type"]>max(Similarity[Distro=="Random"])
#         )/n()*100
#       ) 
#     
#     df_temp <- merge(df_temp, datinfo, by.x="Dataset1", by.y="Dataset")
#     df_temp <- df_temp %>% dplyr::arrange(Percent)
#     df_temp$Plot_Label <- factor(df_temp$Plot_Label, levels=df_temp$Plot_Label)
# 
#     ggplot(df_temp, aes(x=Plot_Label, y=Percent)) +
#       geom_bar(stat="identity", position="stack", linewidth=.4, color="black") +
#       theme_minimal() + 
#       theme(
#         plot.title=element_text(hjust=0, face="bold", size=17),
#         plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
#         legend.title=element_text(size=12),
#         legend.text=element_text(size=10),
#         legend.position="bottom",
#         legend.direction="vertical",
#         legend.box.just="left",
#         legend.key=element_rect(color="black", linewidth=1), 
#         axis.title.x=element_blank(), 
#         axis.text.x=element_text(size=11, angle=45, color="black", hjust=1),
#         axis.ticks.x=element_line(linewidth=.5),
#         axis.title.y=element_text(size=13, margin=margin(r=13)),
#         axis.text.y=element_text(color="black", size=10),
#         panel.border=element_rect(colour="black", fill=NA),
#         panel.grid.major.x=element_blank(),
#         plot.margin=unit(c(1, 1, 1, 1), "cm")
#       ) +
#       scale_y_continuous(limits=c(0, .2))
#     
#     ggplot(df_temp, aes(x=Percent, y=Author_No.Cell_Types)) +
#       geom_smooth(method='lm', formula= y~x, color="black", alpha=.2, show.legend=F) +
#       geom_point(size=.1) +
#       theme_minimal() + 
#       theme(
#         plot.title=element_text(hjust=0, face="bold", size=17),
#         plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
#         legend.title=element_text(size=12),
#         legend.text=element_text(size=10),
#         legend.position="bottom",
#         legend.direction="vertical",
#         legend.box.just="left",
#         legend.key=element_rect(color="black", linewidth=1), 
#         #axis.title.x=element_blank(), 
#         axis.text.x=element_text(size=11, angle=45, color="black", hjust=1),
#         axis.ticks.x=element_line(linewidth=.5),
#         axis.title.y=element_text(size=13, margin=margin(r=13)),
#         axis.text.y=element_text(color="black", size=10),
#         panel.border=element_rect(colour="black", fill=NA),
#         panel.grid.major.x=element_blank(),
#         plot.margin=unit(c(1, 1, 1, 1), "cm")
#       )
#     
#     ggplot(df_temp, aes(x=Percent, y=Author_Median_Unique_Genes_PC)) +
#       geom_smooth(method='lm', formula= y~x, color="black", alpha=.2, show.legend=F) +
#       geom_point(size=.1) +
#       theme_minimal() + 
#       theme(
#         plot.title=element_text(hjust=0, face="bold", size=17),
#         plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
#         legend.title=element_text(size=12),
#         legend.text=element_text(size=10),
#         legend.position="bottom",
#         legend.direction="vertical",
#         legend.box.just="left",
#         legend.key=element_rect(color="black", linewidth=1), 
#         #axis.title.x=element_blank(), 
#         axis.text.x=element_text(size=11, angle=45, color="black", hjust=1),
#         axis.ticks.x=element_line(linewidth=.5),
#         axis.title.y=element_text(size=13, margin=margin(r=13)),
#         axis.text.y=element_text(color="black", size=10),
#         panel.border=element_rect(colour="black", fill=NA),
#         panel.grid.major.x=element_blank(),
#         plot.margin=unit(c(1, 1, 1, 1), "cm")
#       )
#     
#     ggplot(df_temp, aes(x=Plot_Label, y=Percent)) +
#       geom_jitter(size=.1, width=.5) +
#       geom_boxplot(width=.1, outlier.shape=NA, fill="white") +
#       theme_minimal() +
#       theme(
#         plot.title=element_text(hjust=0, face="bold", size=17),
#         plot.subtitle=element_text(hjust=0, size=11, lineheight=1, margin=margin(t=4, b=11)),
#         legend.title=element_text(size=12),
#         legend.text=element_text(size=10),
#         legend.position="bottom",
#         legend.direction="vertical",
#         legend.box.just="left",
#         legend.key=element_rect(color="black", linewidth=1),
#         axis.title.x=element_blank(),
#         axis.text.x=element_text(size=11, angle=45, color="black", hjust=1),
#         axis.ticks.x=element_line(linewidth=.5),
#         axis.title.y=element_text(size=13, margin=margin(r=13)),
#         axis.text.y=element_text(color="black", size=10),
#         panel.border=element_rect(colour="black", fill=NA),
#         panel.grid.major.x=element_blank(),
#         plot.margin=unit(c(1, 1, 1, 1), "cm")
#       )
#   
#   }
#   
# 
# }
