## Ref: https://www.nature.com/articles/s41586-019-1506-7/figures/1
## Ref: https://celltypes.brain-map.org/rnaseq/human_m1_10X?selectedVisualization=Heatmap&colorByFeature=Cell+Type&colorByFeatureValue=GAD1
## Ref: https://celltypes.brain-map.org/rnaseq/human_ctx_smart-seq?selectedVisualization=Heatmap&colorByFeature=Cell+Type&colorByFeatureValue=GAD1

setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/14_neuron_subtypes_annotations")

library(dplyr)

source("/home/rebecca/SCSN_meta_analysis/code/prep_CT_stats_fxn.R")

datinfo <- read.csv("../13_CT_stats/datinfo_SCSN_meta_analysis_CT_stats.csv") %>% na_if("")

ref1 <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/cluster_hierarchy_tableS3.csv")
ref2 <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/cluster_labels_tableS2.csv")
ref3 <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/cluster_hierarchy.csv")

colnames(ref2)[1] <- "Cluster"
colnames(ref2)[4] <- "Level3"

ref2$Cluster <- gsub(" ", "_", ref2$Cluster)
ref2$Cluster <- gsub("-", "_", ref2$Cluster)
ref3$Cluster <- gsub(" ", "_", ref3$Cluster)
ref3$Cluster <- gsub("-", "_", ref3$Cluster)

ct_stats <- prep_CT_stats(datinfo)
ct_stats$Cell_Type <- gsub(" ", "_", gsub("-", "_", ct_stats$Cell_Type, fixed=T))

dim(ct_stats)
# [1] 1410    7

sum(is.element(ct_stats$Cell_Type, ref3$Cluster))
# [1] 241
sum(is.element(ct_stats$Cell_Type, ref2$Cluster))
# [1] 262
sum(is.element(ref3$Cluster, ref2$Cluster))
# [1] 0

ref2$Cluster[!is.element(ref2$Cluster, ct_stats$Cell_Type)]
# [1] 0
ref3$Cluster[!is.element(ref3$Cluster, ct_stats$Cell_Type)]
# [1] 0

ct_stats <- merge(ct_stats, ref3[,c(5, 2)], by.x="Cell_Type", by.y=1, all.x=T)
ct_stats <- merge(ct_stats, ref2[,c(1, 4)], by.x="Cell_Type", by.y=1, all.x=T)

colnames(ct_stats)[grep("Level2", colnames(ct_stats))] <- "Hodge_2018"
colnames(ct_stats)[grep("Level3", colnames(ct_stats))] <- "Bakken_2019"

ct_stats <- ct_stats[is.element(ct_stats$Cell_Class, c("INH", "EXC")),]
ct_stats <- ct_stats %>% dplyr::select(-c(No.Nuclei, Median_UMIs, Median_Unique_Genes))

idx_inh <- grep("INH", ct_stats$Cell_Class)
idx_exc <- grep("EXC", ct_stats$Cell_Class)

###################################################################################################
############################################# Level 1 #############################################
###################################################################################################

ct_stats$Class_Level1 <- NA

## INH: CGE vs. MGE

marker_list <- paste(c("ADARB2", "VIP", "PAX6", "LAMP5", "CGE"), collapse="|")

ct_stats$Cell_Type[intersect(idx_inh, grep(marker_list, ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level1[intersect(idx_inh, grep(marker_list, ct_stats$Cell_Type, ignore.case=T))] <- "CGE"

ct_stats$Class_Level1[grep("ADARB2", ct_stats$Hodge_2018_Annotation, ignore.case=T)] <- "CGE"

marker_list <- paste(c("LHX6", "SST", "PVALB", "PV", "MGE"), collapse="|")

ct_stats$Cell_Type[intersect(idx_inh, grep(marker_list, ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level1[intersect(idx_inh, grep(marker_list, ct_stats$Cell_Type, ignore.case=T))] <- "MGE"

ct_stats$Class_Level1[grep("LHX6", ct_stats$Hodge_2018_Annotation, ignore.case=T)] <- "MGE"

## EXC: FEZF2 vs. RORB

ct_stats$Cell_Type[intersect(idx_exc, grep("FEZF2", ct_stats$Cell_Type))] 
ct_stats$Class_Level1[intersect(idx_exc, grep("FEZF2", ct_stats$Cell_Type, ignore.case=T))] <- "FEZF2"

ct_stats$Cell_Type[intersect(idx_exc, grep("RORB", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level1[intersect(idx_exc, grep("RORB", ct_stats$Cell_Type, ignore.case=T))] <- "RORB"

###################################################################################################
############################################# Level 2 #############################################
###################################################################################################

ct_stats$Class_Level2 <- NA

## INH: LAMP5/PAX6, VIP, SST, PVALB

ct_stats$Cell_Type[intersect(idx_inh, grep("LAMP5|PAX6", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_inh, grep("LAMP5|PAX6", ct_stats$Cell_Type, ignore.case=T))] <- "LAMP5/PAX6"

ct_stats$Cell_Type[intersect(idx_inh, grep("VIP", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_inh, grep("VIP", ct_stats$Cell_Type, ignore.case=T))] <- "VIP"

ct_stats$Cell_Type[intersect(idx_inh, grep("SST", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_inh, grep("SST", ct_stats$Cell_Type, ignore.case=T))] <- "SST"

ct_stats$Cell_Type[intersect(idx_inh, grep("PVALB|PV", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_inh, grep("PVALB|PV", ct_stats$Cell_Type, ignore.case=T))] <- "PVALB"

ct_stats$Class_Level1[grep("LAMP5", ct_stats$Hodge_2018_Annotation, ignore.case=T)] <- "LAMP5/PAX6"
ct_stats$Class_Level1[grep("VIP", ct_stats$Hodge_2018_Annotation, ignore.case=T)] <- "VIP"
ct_stats$Class_Level1[grep("SST", ct_stats$Hodge_2018_Annotation, ignore.case=T)] <- "SST"
ct_stats$Class_Level1[grep("PVALB", ct_stats$Hodge_2018_Annotation, ignore.case=T)] <- "PVALB"

## EXC: L1-6

ct_stats$Cell_Type[intersect(idx_exc, grep("L1_[A-Z]", ct_stats$Cell_Type, ignore.case=T))]

ct_stats$Cell_Type[intersect(idx_exc, grep("L2_[A-Z]", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L2_[A-Z]", ct_stats$Cell_Type, ignore.case=T))] <- "L2"

ct_stats$Cell_Type[intersect(idx_exc, grep("L2_3", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L2_3", ct_stats$Cell_Type, ignore.case=T))] <- "L2-3"

ct_stats$Cell_Type[intersect(idx_exc, grep("L3_[A-Z]|_L3$", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L3_[A-Z]|_L3$", ct_stats$Cell_Type, ignore.case=T))] <- "L3"

ct_stats$Cell_Type[intersect(idx_exc, grep("L2_4", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L2_4", ct_stats$Cell_Type, ignore.case=T))] <- "L2-4"

ct_stats$Cell_Type[intersect(idx_exc, grep("L3_4", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L3_4", ct_stats$Cell_Type, ignore.case=T))] <- "L3-4"

ct_stats$Cell_Type[intersect(idx_exc, grep("L4_[A-Z]|L4$", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L4_[A-Z]|L4$", ct_stats$Cell_Type, ignore.case=T))] <- "L4"

ct_stats$Cell_Type[intersect(idx_exc, grep("L4_5", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L4_5", ct_stats$Cell_Type, ignore.case=T))] <- "L4-5"

ct_stats$Cell_Type[intersect(idx_exc, grep("L3_5", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L3_5", ct_stats$Cell_Type, ignore.case=T))] <- "L3-5"

ct_stats$Cell_Type[intersect(idx_exc, grep("L5_[A-Z]|_L5$", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L5_[A-Z]|_L5$", ct_stats$Cell_Type, ignore.case=T))] <- "L5"

ct_stats$Cell_Type[intersect(idx_exc, grep("L2_6", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L2_6", ct_stats$Cell_Type, ignore.case=T))] <- "L2-6"

ct_stats$Cell_Type[intersect(idx_exc, grep("L4_6", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L4_6", ct_stats$Cell_Type, ignore.case=T))] <- "L4-6"

ct_stats$Cell_Type[intersect(idx_exc, grep("L5_6|L56", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L5_6|L56", ct_stats$Cell_Type, ignore.case=T))] <- "L5-6"

ct_stats$Cell_Type[intersect(idx_exc, grep("L6_[A-Z]|L6$|L6b|L_6", ct_stats$Cell_Type, ignore.case=T))] 
ct_stats$Class_Level2[intersect(idx_exc, grep("L6_[A-Z]|L6$|L6b|L_6", ct_stats$Cell_Type, ignore.case=T))] <- "L6"

ct_stats$Class_Level1[grep("L2/3", ct_stats$Hodge_2018_Annotation, fixed=T)] <- "L2-3"
ct_stats$Class_Level1[grep("L5 ", ct_stats$Hodge_2018_Annotation, fixed=T)] <- "L5"
ct_stats$Class_Level1[grep("L6 ", ct_stats$Hodge_2018_Annotation, fixed=T)] <- "L6"
ct_stats$Class_Level1[grep("L5/6", ct_stats$Hodge_2018_Annotation, fixed=T)] <- "L2-3"

############################################# Save ############################################# 

idx <- apply(ct_stats[,5:ncol(ct_stats)], 1, function(x) sum(!is.na(x))>0)

ct_stats <- ct_stats[idx,]

fwrite(ct_stats, file=paste0("neuron_subtype_annotations_", nrow(ct_stats), "_CTs_", n_distinct(ct_stats$Dataset), "_datasets.csv"))
