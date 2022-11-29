setwd("/home/rebecca/SCSN_meta_analysis/prelim_new/12_cell_class_annotations")

library(dplyr)
library(data.table)

## Add Gugene's cell class annotation info. to author barcode annotations:

datinfo <- read.csv("../11_dataset_mean_expr/datinfo_SCSN_meta_analysis_mean_expr.csv") %>% na_if("")

temp <- read.csv("/home/gugene/home_cluster/consensus/metadata/datinfo_SCSN_meta_analysis.csv") %>% dplyr::select(Dataset, sampleinfo)

temp$Dataset <- gsub("C2T", "CAT", temp$Dataset)

datinfo <- merge(datinfo, temp, by="Dataset") #%>% arrange(idx) %>% na_if("")

colnames(datinfo)[grep("sampleinfo", colnames(datinfo))] <- "MO_Barcode_Annotations"

for(i in 1:nrow(datinfo)){
  
  if(!is.na(datinfo$MO_Barcode_Annotations[i])){
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F)
    
    print(datinfo$Dataset[i])
    
    if(!is.element("MO_Cell_Class", colnames(cellinfo))){
      
      classinfo <- fread(datinfo$MO_Barcode_Annotations[i], data.table=F)

      if(sum(is.element(cellinfo$Cell_ID, classinfo$Cell_ID))!=nrow(cellinfo)){
        stop("sum(is.element(classinfo$Cell_ID, cellinfo$Cell_ID))!=nrow(classinfo)")
      }

      classinfo <- classinfo %>% dplyr::select(Cell_ID, celltype_standard)
      cellinfo <- merge(cellinfo, classinfo, by="Cell_ID", all.x=T)

      ## Confirm order or barcodes matches original order:

      barcodes <- fread(datinfo$Author_Barcodes[i], header=F, data.table=F)
      cellinfo <- cellinfo[match(barcodes[,1], cellinfo[,1]),]

      if(!identical(barcodes[,1], cellinfo$Cell_ID)){
        stop("!identical(barcodes[,1], cellinfo$Cell_ID)")
      }

      colnames(cellinfo)[ncol(cellinfo)] <- "MO_Cell_Class"

      cellinfo$MO_Cell_Class <- toupper(cellinfo$MO_Cell_Class)
      cellinfo$MO_Cell_Class <- gsub("AST", "ASC", cellinfo$MO_Cell_Class)
      cellinfo$MO_Cell_Class <- gsub("OLI", "OG", cellinfo$MO_Cell_Class)
      cellinfo$MO_Cell_Class <- gsub("PERI", "PER", cellinfo$MO_Cell_Class)
      cellinfo$MO_Cell_Class[grep("NEU", cellinfo$MO_Cell_Class)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("DOP", cellinfo$MO_Cell_Class)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("RELN", cellinfo$MO_Cell_Class)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("RELN", cellinfo$MO_Cell_Class)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("PURK", cellinfo$MO_Cell_Class)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("MSN", cellinfo$MO_Cell_Class)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("MSN", cellinfo$Cell_Type)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("SOX6", cellinfo$Cell_Type)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("Neu_FAT2", cellinfo$Cell_Type)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("neurons", cellinfo$Cell_Type)] <- "NEU"
      cellinfo$MO_Cell_Class[grep("GABA", cellinfo$Cell_Type)] <- "INH"
      cellinfo$MO_Cell_Class[grep("Glut", cellinfo$Cell_Type)] <- "EXC"
      cellinfo$MO_Cell_Class[grep("granule", cellinfo$Cell_Type)] <- "EXC"
      cellinfo$MO_Cell_Class[grep("NSC", cellinfo$MO_Cell_Class)] <- "EPEN"
      
      ## Fix a couple of inconsistencies:

      if(grepl("Jakel", datinfo$Dataset[i])){
        cellinfo$MO_Cell_Class[grep("COP", cellinfo$Cell_Type)] <- "COP"
      }

    }
    
    fwrite(cellinfo, file=paste0("data/author_data/", datinfo$Dataset[i], "_barcode_annotations_cell_class_annotated.csv"))
    
  }
  
}

## Add Tran data (Gugene made his annotations with older cell IDs I provided):

datinfo <- read.csv("../11_dataset_mean_expr/datinfo_SCSN_meta_analysis_mean_expr.csv") %>% na_if("")

idx <- grep("Tran", datinfo$First_Author)

for(i in idx){
  
  cat("\n")
  print(datinfo$Dataset[i])
  
  cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F)
  
  cellinfo$MO_Cell_Class <- NA
  cellinfo$MO_Cell_Class[grep("Ast", cellinfo$Cell_Type)] <- "ASC"
  cellinfo$MO_Cell_Class[grep("Endo", cellinfo$Cell_Type)] <- "END"
  cellinfo$MO_Cell_Class[grep("Excit", cellinfo$Cell_Type)] <- "EXC"
  cellinfo$MO_Cell_Class[grep("Inh", cellinfo$Cell_Type)] <- "INH"
  cellinfo$MO_Cell_Class[grep("Mic", cellinfo$Cell_Type)] <- "MIC"
  cellinfo$MO_Cell_Class[grep("Mur", cellinfo$Cell_Type)] <- "MUR"
  cellinfo$MO_Cell_Class[grep("Oligo", cellinfo$Cell_Type)] <- "OG"
  cellinfo$MO_Cell_Class[grep("OPC", cellinfo$Cell_Type)] <- "OPC"
  cellinfo$MO_Cell_Class[grep("Tcell", cellinfo$Cell_Type)] <- "IMM"
  cellinfo$MO_Cell_Class[grep("Macro", cellinfo$Cell_Type)] <- "IMM"
  cellinfo$MO_Cell_Class[grep("Neu_FAT2", cellinfo$Cell_Type)] <- "NEU"
  
  print(table(cellinfo$MO_Cell_Class, useNA=c("ifany")))
  
  print(table(cellinfo$Cell_Type[is.na(cellinfo$MO_Cell_Class)]))
  
  fwrite(cellinfo, file=paste0("data/author_data/", datinfo$Dataset[i], "_barcode_annotations_cell_class_annotated.csv"))
  
}

## Add info for Yang 2021:

datinfo <- read.csv("../11_dataset_mean_expr/datinfo_SCSN_meta_analysis_mean_expr.csv") %>% na_if("")

idx <- grep("Yang_2021_10x_Chromium_V3_Medial_frontal_gyrus", datinfo$Dataset)

cellinfo <- fread(datinfo$Author_Barcode_Annotations[idx], data.table=F)

cellinfo$MO_Cell_Class <- cellinfo$Cell_Class
cellinfo$MO_Cell_Class[grep("Ast", cellinfo$Cell_Class)] <- "ASC"
cellinfo$MO_Cell_Class[grep("Endo", cellinfo$Cell_Class)] <- "END"
cellinfo$MO_Cell_Class[grep("Exc", cellinfo$Cell_Class)] <- "EXC"
cellinfo$MO_Cell_Class[grep("Inh", cellinfo$Cell_Class)] <- "INH"
cellinfo$MO_Cell_Class[grep("Mic", cellinfo$Cell_Class)] <- "MIC"
cellinfo$MO_Cell_Class[grep("NRGN", cellinfo$Cell_Class)] <- "NEU"
cellinfo$MO_Cell_Class[grep("Oli", cellinfo$Cell_Class)] <- "OG"

# table(cellinfo$Cell_Type[is.na(cellinfo$MO_Cell_Class)])
# table(cellinfo$MO_Cell_Class, useNA="ifany")

fwrite(cellinfo, file=paste0("data/author_data/", datinfo$Dataset[idx], "_barcode_annotations_cell_class_annotated.csv"))

# ## Look for EPEND cell types 
# 
# datinfo <- read.csv("../11_dataset_mean_expr/datinfo_SCSN_meta_analysis_mean_expr.csv") %>% na_if("") %>% dplyr::filter(!is.na(Author_Barcode_Annotations))
# 
# cellinfo_list <- lapply(1:nrow(datinfo), function(i){
#   cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F)
#   if(is.element("MO_Cell_Class", colnames(cellinfo))){
#     print(datinfo$Dataset[i])
#     cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F)
#     cellinfo <- cellinfo %>% dplyr::select(
#       Cell_ID, Cell_Type, MO_Cell_Class
#     )
#     cellinfo$Dataset <- datinfo$Dataset[i]
#     return(cellinfo)
#   }
# })
# cellinfo <- do.call(rbind, cellinfo_list) %>% na_if("")
# 
# table(cellinfo$Cell_Type[is.na(cellinfo$MO_Cell_Class)])

## Add datinfo paths:

rm(list=ls())

datinfo <- read.csv("../11_dataset_mean_expr/datinfo_SCSN_meta_analysis_mean_expr.csv") %>% na_if("")

datinfo$Author_Barcode_Annotations_Class <- paste0(getwd(), "/data/author_data/", datinfo$Dataset, "_barcode_annotations_cell_class_annotated.csv")

datinfo$Author_Barcode_Annotations_Class[is.na(datinfo$Author_Barcode_Annotations)] <- NA

idx <- which(!is.na(datinfo$Author_Barcode_Annotations))

sum(sapply(datinfo$Author_Barcode_Annotations_Class[idx], file.exists)==F)

## Make sure cell annotations look ok:

for(i in 1:nrow(datinfo)){
  
  if(!is.na(datinfo$Author_Barcode_Annotations_Class[i])){
    
    cat("\n")
    print(datinfo$Dataset[i])
    
    cellinfo <- fread(datinfo$Author_Barcode_Annotations[i], data.table=F) %>% na_if("")
    
    # print(table(cellinfo$Cell_Type[is.na(cellinfo$MO_Cell_Class)]))
    
    print(table(cellinfo$Cell_Type[!is.element(cellinfo$MO_Cell_Class, c("ASC", "END", "EXC", "INH", "MIC", "NEU", "PER", "OG", "OPC", "VSMC"))]))
    
  }
  
}

fwrite(datinfo, file="datinfo_SCSN_meta_analysis_cell_class.csv")
