setwd("/home/rebecca/SCSN_meta_analysis/prelim/1_gene_sets")

library(dplyr)

source("/home/rebecca/code/misc/prep_gene_sets/prep_gene_sets_fxns.R")

set_dir <- "/home/rebecca/SCSN_meta_analysis/gene_sets/SCSN_gene_sets/"
mapping_tables_dir <- "/home/rebecca/omicon/mapping_tables/"
out_dir <- "/home/rebecca/SCSN_meta_analysis/gene_sets"
projectname <- "SCSN"
n_threads <- 5

MO_legend <- read.csv("/home/rebecca/SCSN_meta_analysis/gene_sets/SCSN_gene_sets/MyGeneSetsLEGEND.csv") %>% na_if("")

prep_MO_sets(projectname, set_dir, MO_legend, out_dir)

load("/home/rebecca/SCSN_meta_analysis/gene_sets/SCSN_sets.RData")

all.equal(MO_legend$SetID, names(MO_sets))

map_sets(projectname, gene_sets=MO_sets, legend=MO_legend, mapping_tables_dir, out_dir, n_threads)
