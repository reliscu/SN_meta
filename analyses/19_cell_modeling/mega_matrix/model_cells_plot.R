setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/19_cell_modeling/mega_matrix")

library(data.table)

source("model_cells_plot_fxns.R")

data_type <- "author_data"
expr_type <- "normalized_counts"
summary_type <- "mean"
n_genes <- 8596

############################################# Plot coefficients #############################################

cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU")
n_pred <- length(cell_classes)
df <- fread(Sys.glob(path=file.path("data", data_type, expr_type, paste0("*coef*", paste(cell_classes, collapse="_"), "_", summary_type, "*"))), data.table=F) 

############################################# Predictor group 1 #############################################

cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU")
n_pred <- length(cell_classes)
df1 <- fread(Sys.glob(path=file.path("data", data_type, expr_type, paste0("*variance*", paste(cell_classes, collapse="_"), "_", summary_type, "*"))), data.table=F)
df1 <- reshape2::melt(df1, value.name="R2", variable.name="Distro")
df1$Predictors <- paste(cell_classes, collapse=" | ")

df <- df1

############################################# Predictor group 2 #############################################

cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "EXC", "INH")
n_pred <- length(cell_classes)
df2 <- fread(Sys.glob(path=file.path("data", data_type, expr_type, paste0("*variance*", paste(cell_classes, collapse="_"), "_", summary_type, "*"))), data.table=F)
df2 <- reshape2::melt(df2, value.name="R2", variable.name="Distro")
df2$Predictors <- paste(cell_classes, collapse=" | ")

############################################# Predictor group 3 #############################################

cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU", "EXC", "INH", "CGE", "MGE")
n_pred <- length(cell_classes)
df3 <- fread(Sys.glob(path=file.path("data", data_type, expr_type, paste0("*variance*", paste(cell_classes, collapse="_"), "_", summary_type, "*"))), data.table=F)
df3 <- reshape2::melt(df3, value.name="R2", variable.name="Distro")
df3$Predictors <- paste(cell_classes, collapse=" | ")

############################################# Combine #############################################

df <- rbind(df1, df2, df3)

############################################# Sequential predictors 1 #############################################

pred_list <- readRDS(Sys.glob(path=file.path("data", data_type, expr_type, paste0("*variance*", summary_type, "*RDS"))))
df <- reshape2::melt(pred_list, value.name="R2", variable.name="Predictors") %>% select(-c(L1))

############################################# Sequential predictors 2 #############################################

############################################# Left out data 1 #############################################

cell_classes <- c("ASC", "END", "MIC", "OG", "OPC", "PER", "VSMC", "NEU")
mdl_list <- readRDS(Sys.glob(path=file.path("data", data_type, expr_type, paste0("modeling_", paste(cell_classes, collapse="_"), "_", summary_type, "cell_*left_out*"))))
df <- do.call(rbind, mdl_list)
rownames(df) <- 1:nrow(df)

n_pred <- length(cell_classes)
