setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/4_marker_gene_comparison/marker_distribution")

source("marker_distribution_plot.R")

datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv")

cell_classes <- c("ASC", "END", "EXC", "INH", "MIC", "OG", "OPC", "PER", "VSMC")

remove_subclust <- T

## Restrict to cortical datasets:

datinfo <- datinfo %>% dplyr::filter(grepl("C$", Region_Code))

legend <- read.csv("/home/rebecca/SCSN_meta_analysis/gene_sets/SCSN_sets_legend_metadata.csv")

#plot_marker_distro(datinfo, legend)

remove_subclust <- T

gene_sets <- gene_sets_mapped

jaccard <- T
min_size <- 5

marker_overlap_heatmap(datinfo, cell_classes, legend, gene_sets, jaccard, pc_genes, min_size, remove_subclust=T)
