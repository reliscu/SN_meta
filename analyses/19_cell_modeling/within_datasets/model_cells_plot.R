setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/19_cell_modeling/within_datasets")

library(data.table)

source("model_cells_plot_fxns.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/prelim_new/15_CT_summary_vecs/datinfo_SCSN_meta_analysis_summary_vecs.csv") %>% na_if("") 

data_type <- "author_data"
expr_type <- "normalized_counts"
summary_type <- "mean"
n_genes <- 8596

############################################# Predictor group 1 #############################################

cell_classes <- c("ASC", "MIC", "NEU", "OG", "OPC")
n_pred <- length(cell_classes)
df_paths <- Sys.glob(path=file.path("data", data_type, expr_type, paste0("*variance*", paste(cell_classes, collapse="_"), "_", summary_type, "*")))

df_list <- lapply(1:length(df_paths), function(i){
  temp <- fread(df_paths[i], data.table=F)
  if(sum(is.element(cell_classes, temp$MO_Cell_Class))==length(cell_classes)){
    ## Make sure columns for all data match:
    idx <- which(!is.element(colnames(temp), "Cell_Class"))
    temp <- temp[,idx]
    temp <- temp[,!grepl("Level", colnames(temp))]
    dataset <- sapply(strsplit(sapply(strsplit(df_paths[i], "/"), "[", 4), "_cell"), "[", 1)
    return(data.frame(Dataset=dataset, temp))
  }
})
df1 <- do.call(rbind, df_list)
df1 <- reshape2::melt(df1, value.name="R2", variable.name="Distro")
df1$Predictors <- paste(c("ASC", "MIC", "OG", "OPC", "NEU"), collapse=" | ")

df <- df1

# plot_r2_real_vs_rand(datinfo, df, data_type, expr_type, summary_type, n_pred, n_genes, order_by="real_r2")
# plot_r2_real_vs_rand(datinfo, df, data_type, expr_type, summary_type, n_pred, n_genes, order_by="delta")
# r2_scatterplot(datinfo, df, data_type, expr_type, summary_type, n_pred, n_genes)

############################################# Predictor group 2 #############################################

cell_classes <- c("ASC", "EXC", "INH", "MIC", "OG", "OPC") 
n_pred <- length(cell_classes)
df_paths <- Sys.glob(path=file.path("data", data_type, expr_type, paste0("*variance*", paste(cell_classes, collapse="_"), "_", summary_type, "*")))

df_list <- lapply(1:length(df_paths), function(i){
  temp <- fread(df_paths[i], data.table=F)
  if(sum(is.element(cell_classes, temp$MO_Cell_Class))==length(cell_classes)){
    ## Make sure columns for all data match:
    idx <- which(!is.element(colnames(temp), "Cell_Class"))
    temp <- temp[,idx]
    temp <- temp[,!grepl("Level", colnames(temp))]
    dataset <- sapply(strsplit(sapply(strsplit(df_paths[i], "/"), "[", 4), "_cell"), "[", 1)
    return(data.frame(Dataset=dataset, temp))
  }
})
df2 <- do.call(rbind, df_list)
df2 <- reshape2::melt(df2, value.name="R2", variable.name="Distro")
df2$Predictors <- paste(c("ASC", "MIC", "OG", "OPC", "EXC", "INH"), collapse=" | ")

############################################# Predictor group 3 #############################################

cell_classes <- c("ASC", "CGE", "EXC", "MGE", "MIC", "OG", "OPC") 
n_pred <- length(cell_classes)
df_paths <- Sys.glob(path=file.path("data", data_type, expr_type, paste0("*variance*", paste(cell_classes, collapse="_"), "_", summary_type, "*")))

df_list <- lapply(1:length(df_paths), function(i){
  temp <- fread(df_paths[i], data.table=F)
  if(sum(is.element(cell_classes, temp$MO_Cell_Class))==length(cell_classes)){
    ## Make sure columns for all data match:
    idx <- which(!is.element(colnames(temp), "Cell_Class"))
    temp <- temp[,idx]
    temp <- temp[,!grepl("Level", colnames(temp))]
    dataset <- sapply(strsplit(sapply(strsplit(df_paths[i], "/"), "[", 4), "_cell"), "[", 1)
    return(data.frame(Dataset=dataset, temp))
  }
})
df3 <- do.call(rbind, df_list)
df3 <- reshape2::melt(df3, value.name="R2", variable.name="Distro")
df3$Predictors <- paste(c("ASC", "MIC", "OG", "OPC", "EXC", "CGE", "MGE"), collapse=" | ")

############################################# Combine #############################################

# datasets <- intersect(unique(df1$Dataset), unique(df2$Dataset))
# 
# df1 <- df1[is.element(df1$Dataset, datasets),]
# df2 <- df2[is.element(df2$Dataset, datasets),]
# 
# df <- rbind(df1, df2)

datasets <- intersect(unique(df3$Dataset), intersect(unique(df1$Dataset), unique(df2$Dataset)))

df1 <- df1[is.element(df1$Dataset, datasets),]
df2 <- df2[is.element(df2$Dataset, datasets),]
df3 <- df3[is.element(df3$Dataset, datasets),]

df <- rbind(df1, df2, df3)

plot_r2_predictor_groups(datinfo, df, data_type, expr_type, summary_type, n_genes)
