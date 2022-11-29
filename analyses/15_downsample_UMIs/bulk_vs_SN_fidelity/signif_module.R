setwd("/home/rebecca/SCSN_meta_analysis/analyses/cortical_regions/15_downsample_UMIs/bulk_vs_SN_fidelity")

source("/home/rebecca/SCSN_meta_analysis/code/signif_module_fxn.R")

network_path <- "/mnt/bdata/gugene/megadataset/cortex_megaset/cortex_Combat_by_platform_color_by_platform_Modules"
set_string <- "neuron"

kme <- signif_module_kme(network_path, set_string)

fwrite(kme, file=paste0("data/consensus_kME_", kme$SetName[1], "_enriched_", signif(kme$Pval[1], 2), "_pval_", kme$No.Samples[1], "_samples.csv"))
