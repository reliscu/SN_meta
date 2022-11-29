fidelity <- function(expr, ct_idx, prop_scaled=F){
  
  expr_ct <- as.matrix(expr[,ct_idx])
  sensitivity <- apply(expr_ct, 1, function(x) sum(x>0))/length(ct_idx)
  specificity <- rowSums(expr_ct)/rowSums(expr)
  fid <- sensitivity*specificity
  
  if(prop_scaled){
    
    ## Proportion of expression from gene vs. all genes in cell type cells:
    prop <- rowSums(expr_ct)/sum(colSums(expr_ct)) #prop <- prop/max(prop)
    fid <- fid*prop
    
  }
  
  return(fid)
  
}
