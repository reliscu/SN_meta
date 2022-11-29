library(dplyr)
library(data.table)

cell_class_comp <- function(
  datinfo, 
  cell_classes=c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC"), 
  na_cts=c("Unassigned", "drop.*", "Mix_", "u7", "Outlier", "*-MIX", "no class")
){
  
  cell_comp_df <- c()
  
  for(i in 1:nrow(datinfo)){
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F) %>% na_if("")
    
    if(is.element("MO_Cell_Class", colnames(cellinfo))){
      
      ## Don't include NA/ambiguous cell types in total CT counts:
      
      cellinfo <- cellinfo[!is.na(cellinfo$Cell_Type),]
      cellinfo <- cellinfo[!grepl(paste(na_cts, collapse="|"), cellinfo$Cell_Type),]
      cellinfo$MO_Cell_Class[!is.element(cellinfo$MO_Cell_Class, cell_classes)] <- "Other" 
      cell_comp <- data.frame(
        table(cellinfo$MO_Cell_Class)/nrow(cellinfo)
      )
      cell_comp$Dataset <- datinfo$Dataset[i]
      
      if(is.null(cell_comp_df)){
        cell_comp_df <- cell_comp
      } else {
        cell_comp_df <- rbind(cell_comp_df, cell_comp)
      }
      
    } ## if(is.element("MO_Cell_Class", colnames(cellinfo))){
    
  } ## for(i in 1:nrow(datinfo)){
  
  colnames(cell_comp_df) <- c("Cell_Class", "Proportion", "Dataset")
  
  fwrite(cell_comp_df, file=paste0("data/", data_type, "/cell_class_composition_", n_distinct(cell_comp_df$Dataset), "_datasets.csv"))
  
} ## cell_class_proportion <- function(