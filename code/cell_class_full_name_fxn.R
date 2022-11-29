cell_class_full_name <- function(cell_class){
  
  cell_class <- gsub("ASC", "Astrocyte", cell_class)
  cell_class <- gsub("END", "Endothelial", cell_class)
  cell_class <- gsub("EXC", "Excitatory Neuron", cell_class)
  cell_class <- gsub("INH", "Inhibitory Neuron", cell_class)
  cell_class <- gsub("OG", "Oligodendrocyte", cell_class)
  cell_class <- gsub("MIC", "Microglia", cell_class)
  cell_class <- gsub("PER", "Pericyte", cell_class)
  cell_class <- gsub("NEU", "Neuron", cell_class)
  
  return(cell_class)
}