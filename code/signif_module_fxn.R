library(dplyr)
library(data.table)

signif_module_kme <- function(network_path, set_string){
  
  networks <- list.files(path=network_path, pattern="signum", full.names=T)
  idx <- unlist(lapply(networks, function(x) length(list.files(x))))>0
  networks <- networks[idx]
  
  mod_list <- lapply(1:length(networks), function(i){
    
    enrich_paths <- list.files(path=networks[i], pattern="GSHyperG_CONSENSUS", full.names=T)
    
    for(j in 1:length(enrich_paths)){
      
      enrich <- fread(enrich_paths[j], data.table=F)
      enrich <- enrich[grep(set_string, enrich$SetName, ignore.case=T),]
      
      if(nrow(enrich)>0){
        
        pvals <- enrich[,-c(1:3)]
        setnames <- enrich[,c(2)]
        mods <- colnames(enrich)[-c(1:3)]
        
        if(is.null(nrow(pvals))){
          
          sig_pvals <- pvals[pvals<.05]
          setnames <- setnames[pvals<.05]
          
        } else { ## if(is.null(nrow(pvals))){
          
          sig_pvals <- apply(pvals, 1, function(x) x[x<.05])
          mods <- unlist(apply(pvals, 1, function(x) mods[x<.05]))
          
        } ## if(is.null(nrow(pvals))){} else {
        
        if(length(sig_pvals)>0){
          
          setnames <- rep(setnames, lengths(sig_pvals))
          sig_pvals <- unlist(sig_pvals)
          module_def <- gsub(".csv", "", sapply(strsplit(enrich_paths[j], "_"), function(x) x[length(x)]))
          
          return(data.frame(SetName=setnames, Network=networks[i], Module=mods, Module_Def=module_def, Pval=sig_pvals))
          
        }
        
      } ## for(j in 1:length(enrich_paths)){
      
    } ## for(j in 1:length(enrich_paths)){

  }) ## for(i in 1:length(networks)){
  
  df <- do.call(rbind, mod_list) %>% slice_min(Pval)
  
  mod_eig <- fread(list.files(path=df$Network, pattern="eigengenes", full.names=T), data.table=F)
  
  kme <- fread(list.files(path=df$Network, pattern="kME", full.names=T), data.table=F)
  kme <- kme[,c(1:5, which(is.element(colnames(kme), paste0("kME", df$Module))))]
  kme$No.Samples <- nrow(mod_eig)
  kme$SetName <- df$SetName
  kme$Pval <- df$Pval
  
  return(kme)
  
} ## signif_module_kme <- function(
