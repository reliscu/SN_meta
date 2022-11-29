accuracy <- function(cell_class, true_labs, clust_labs){
  
  clust_list <- split(true_labs, factor(clust_labs))
  
  ## Total # of class X members:
  
  total_class <- sum(sapply(clust_list, function(clust) sum(is.element(clust, cell_class))))
  
  ## Total # of non-class X members across all clusters:
  
  not_cell_class <- cell_classes[!is.element(cell_classes, cell_class)]
  total_not_class <- sum(sapply(clust_list, function(clust) sum(is.element(clust, not_cell_class))))
  
  ## Which cluster has maximum # of class X members?
  
  clust_idx <- which.max(sapply(clust_list, function(clust) sum(is.element(clust, cell_class))))
  
  # No. of class X clustered together (true positives)
  
  tp <- sum(is.element(clust_list[[clust_idx]], cell_class))
  
  # No. of NON-class X clustered with class X
  
  fp <- length(clust_list[[clust_idx]])-tp
  
  # No. of class X in a different cluster (false negatives)
  
  fn <- total_class-tp
  
  # No. NON-class X in a different cluster as class X (true negatives)
  
  tn <- total_not_class-fp
  acc <- (tp+tn)/(tp+tn+fp+fn)
  
  return(acc)
  
} ## accuracy <- function(

purity <- function(clust_labs, true_labs) {
  
  ## Calculates the purity of a clustering result given in x, w.r.t. true cluster/class labels given in y
  
  clust_list <- split(true_labs, factor(clust_labs))
  
  purity_idx <- sum(sapply(clust_list, function(z){max(sapply(split(z, factor(z)), length))}))/length(clust_labs)
  
  ## https://www.reddit.com/r/rstats/comments/28wa4y/calculate_the_purity_of_a_clustering_results_with/
  
}