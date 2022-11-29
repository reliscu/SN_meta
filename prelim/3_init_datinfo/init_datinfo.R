datinfo <- data.frame(
  First_Author=NA,
  Year=NA,
  PubMedID=NA,
  Disease=NA,
  Region=NA,
  Subregion=NA,
  Author_Region=NA,
  Region_Code=NA,
  FACS_Sorted=NA,
  Unbiased_Sampling=NA,
  Platform=NA,
  Directory=NA,
  Author_Counts=NA,
  Author_Barcodes=NA,
  Author_Genes=NA,
  Author_Barcode_Annotations=NA
)


########################################################################################################
############################################# Agarwal 2020 #############################################
########################################################################################################

datinfo$First_Author[1] <- "Agarwal"
datinfo$Year[1] <- 2020
datinfo$PubMedID[1] <- 32826893
datinfo$Disease[1] <- "Normal"
datinfo$Region[1] <- "Cortex"
datinfo$Subregion[1] <- "Frontal cortex"
datinfo$Author_Region[1] <- "Middle frontal gyrus"
datinfo$Platform[1] <- "10x Chromium V2"
datinfo$Directory[1] <- "/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/"
datinfo$Region_Code[1]="FC"
datinfo$FACS_Sorted[1]="Y"
datinfo$Unbiased_Sampling[1]="Y"
datinfo$Author_Counts[1] <- "/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/cortex/matrix_cortex_annotated.mtx"
datinfo$Author_Barcodes[1] <- "/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/cortex/barcodes_cortex_annotated.tsv"
datinfo$Author_Genes[1] <- "/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/cortex/genes_cortex.tsv"
datinfo$Author_Barcode_Annotations[1] <- "/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/cortex/author_barcode_annotations_cortex.csv"

datinfo <- rbind(datinfo, data.frame(First_Author="Agarwal", Year=2020, PubMedID=32826893, Disease="Normal", Region="Basal ganglia", Subregion="Subsantia nigra", Author_Region="Subsantia nigra", Region_Code="MID", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="10x Chromium V2", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/SN/matrix_SN_annotated.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/SN/barcodes_SN_annotated.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/SN/genes_SN.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/agarwal_2020/SN/author_barcode_annotations_SN.csv"))

###################################################################################################
############################################# ABI LGN 2018 ########################################
###################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="ABI", Year=2018, PubMedID=NA, Disease="Normal", Region="Forebrain", Subregion="Thalamus", Author_Region="Lateral geniculate nucleus", Region_Code="DI", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="SMART-Seq v4", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2018_LGN/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2018_LGN/expression.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2018_LGN/barcodes.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2018_LGN/genes.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2018_LGN/author_barcode_annotations.csv"))

################################################################################################################
############################################# ABI cortical regions 2019 ########################################
################################################################################################################

## ACC

datinfo <- rbind(datinfo, data.frame(First_Author="ABI", Year=2019, PubMedID=NA, Disease="Normal", Region="Cortex", Subregion="Anterior cingulate cortex", Author_Region="Anterior cingulate cortex", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="SMART-Seq v4", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/ACC/expression_ACC.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/ACC/barcodes_ACC.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/ACC/genes_ACC.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/ACC/author_barcode_annotations_ACC.csv"))

## A1

datinfo <- rbind(datinfo, data.frame(First_Author="ABI", Year=2019, PubMedID=NA, Disease="Normal", Region="Cortex", Subregion="Temporal lobe", Author_Region="A1", Region_Code="TC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="SMART-Seq v4", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/A1/expression_A1.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/A1/barcodes_A1.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/A1/genes_A1.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/A1/author_barcode_annotations_A1.csv"))

## S1

datinfo <- rbind(datinfo, data.frame(First_Author="ABI", Year=2019, PubMedID=NA, Disease="Normal", Region="Cortex", Subregion="Parietal lobe", Author_Region="S1", Region_Code="SC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="SMART-Seq v4", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/S1/expression_S1.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/S1/barcodes_S1.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/S1/genes_S1.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/S1/author_barcode_annotations_S1.csv"))

## V1

datinfo <- rbind(datinfo, data.frame(First_Author="ABI", Year=2019, PubMedID=NA, Disease="Normal", Region="Cortex", Subregion="Occipital lobe", Author_Region="V1", Region_Code="OC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="SMART-Seq v4", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/V1/expression_V1.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/V1/barcodes_V1.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/V1/genes_V1.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/ABI_2019_cortical_regions/V1/author_barcode_annotations_V1.csv"))

######################################################################################################
############################################# Bakken 2019 #############################################
######################################################################################################

## SSv4

datinfo <- rbind(datinfo, data.frame(First_Author="Bakken", Year=2019, PubMedID=34616062, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="M1", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="SMART-Seq v4", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/SSv4/NeMO/matrix_SSv4_annotated.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/SSv4/NeMO/barcodes_SSv4_annotated.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/SSv4/NeMO/genes_SSv4.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/SSv4/NeMO/author_barcode_annotations_SSv4.csv"))

## 10x

datinfo <- rbind(datinfo, data.frame(First_Author="Bakken", Year=2019, PubMedID=34616062, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="M1", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/10x/NeMO/matrix_10x.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/10x/NeMO/barcodes_10x.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/10x/NeMO/genes_10x.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/10x/NeMO/author_barcode_annotations_10x.csv"))

######################################################################################################
############################################# Habib 2017 #############################################
######################################################################################################

## PFC

datinfo <- rbind(datinfo, data.frame(First_Author="Habib", Year=2017, PubMedID=28846088, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="Prefrontal cortex", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="DroNc-seq", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/PFC/expression_PFC.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/PFC/barcodes_PFC.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/PFC/genes_PFC.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/PFC/author_barcode_annotations_PFC.csv"))

## Hippocampus

datinfo <- rbind(datinfo, data.frame(First_Author="Habib", Year=2017, PubMedID=28846088, Disease="Normal", Region="Hippocampus", Subregion="Hippocampus", Author_Region="CA1 | CA3 | Dentate gyrus", Region_Code="HIP", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="DroNc-seq", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/hippocampus/expression_hippocampus.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/hippocampus/barcodes_hippocampus.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/hippocampus/genes_hippocampus.tsv",  Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/habib_2017/hippocampus/author_barcode_annotations_hippocampus.csv"))

######################################################################################################
############################################# Hodge 2018 #############################################
######################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Hodge", Year=2018, PubMedID=31435019, Disease="Normal", Region="Cortex", Subregion="Temporal lobe", Author_Region="Middle temporal gyrus", Region_Code="TC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="SMART-Seq v4", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/ABI/expression.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/ABI/barcodes.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/ABI/genes.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/hodge_2018/ABI/author_barcode_annotations.csv"))

##########################################################################################################
############################################# Khrameeva 2020 #############################################
##########################################################################################################

## ACC

datinfo <- rbind(datinfo, data.frame(First_Author="Khrameeva", Year=2020, PubMedID=32424074, Disease="Normal", Region="Cortex", Subregion="Anterior cingulate cortex", Author_Region="BA24", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="Y",Platform="10x Chromium V2", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/ACC/matrix_ACC.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/ACC/barcodes_ACC.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/ACC/genes_ACC.tsv", Author_Barcode_Annotations=NA))

## Cerebellum

datinfo <- rbind(datinfo, data.frame(First_Author="Khrameeva", Year=2020, PubMedID=32424074, Disease="Normal", Region="Cerebellum", Subregion="Cerebellum", Author_Region="Cerebellum", Region_Code="CB", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="10x Chromium V2", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/cerebellum/matrix_cerebellum.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/cerebellum/barcodes_cerebellum.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/cerebellum/genes_cerebellum.tsv", Author_Barcode_Annotations=NA))

## Caudate nucleus

datinfo <- rbind(datinfo, data.frame(First_Author="Khrameeva", Year=2020, PubMedID=32424074, Disease="Normal", Region="Basal ganglia", Subregion="Striatum", Author_Region="Caudate nucleus", Region_Code="STR", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="10x Chromium V2", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/CN/matrix_CN.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/CN/barcodes_CN.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/khrameeva_2020/CN/genes_CN.tsv", Author_Barcode_Annotations=NA))

#####################################################################################################
############################################# Lake 2018 #############################################
#####################################################################################################

## BA6

datinfo <- rbind(datinfo, data.frame(First_Author="Lake", Year=2018, PubMedID=29227469, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="BA6", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="snDrop-seq", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/FCX/BA6/expression_BA6.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/FCX/BA6/barcodes_BA6.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/FCX/BA6/genes_BA6.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/FCX/BA6/author_barcode_annotations_BA6.csv"))

## BA10

datinfo <- rbind(datinfo, data.frame(First_Author="Lake", Year=2018, PubMedID=29227469, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="BA10", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="snDrop-seq", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/FCX/BA10/expression_BA10.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/FCX/BA10/barcodes_BA10.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/FCX/BA10/genes_BA10.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/FCX/BA10/author_barcode_annotations_BA10.csv"))

## V1

datinfo <- rbind(datinfo, data.frame(First_Author="Lake", Year=2018, PubMedID=29227469, Disease="Normal", Region="Cortex", Subregion="V1", Author_Region="BA17", Region_Code="OC", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="snDrop-seq", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/V1/expression_V1.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/V1/barcodes_V1.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/V1/genes_V1.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/V1/author_barcode_annotations_V1.csv"))

## Cerebellum

datinfo <- rbind(datinfo, data.frame(First_Author="Lake", Year=2018, PubMedID=29227469, Disease="Normal", Region="Cerebellum", Subregion="Cerebellum", Author_Region="Lateral cerebellar hemisphere", Region_Code="CB", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="snDrop-seq", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/cerebellum/expression_cerebellum.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/cerebellum/barcodes_cerebellum.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/cerebellum/genes_cerebellum.tsv",  Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/lake_2018/cerebellum/author_barcode_annotations_cerebellum.csv"))

####################################################################################################
############################################# Lee 2020 #############################################
####################################################################################################

## Caudate nucleus

datinfo <- rbind(datinfo, data.frame(First_Author="Lee", Year=2020, PubMedID=32681824, Disease="Normal", Region="Basal ganglia", Subregion="Striatum", Author_Region="Caudate nucleus", Region_Code="STR", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/CN/matrix_CN_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/CN/barcodes_CN_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/CN/genes_CN_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/CN/author_barcode_annotations_CN.csv"))

## Putamen

datinfo <- rbind(datinfo, data.frame(First_Author="Lee", Year=2020, PubMedID=32681824, Disease="Normal", Region="Basal ganglia", Subregion="Striatum", Author_Region="Putamen", Region_Code="STR", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/putamen/matrix_putamen_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/putamen/barcodes_putamen_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/putamen/genes_putamen_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/lee_2020/putamen/author_barcode_annotations_putamen.csv"))

####################################################################################################
############################################# Li 2018 #############################################
####################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Li", Year=2018, PubMedID=30545854, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="BA9", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V2", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/li_2018/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/li_2018/matrix.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/li_2018/barcodes.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/li_2018/genes.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/li_2018/author_barcode_annotations.csv"))

#######################################################################################################
############################################# Mathys 2019 #############################################
#######################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Mathys", Year=2019, PubMedID=31042697, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="BA10", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/mathys_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/mathys_2019/skinnider_2021/controls/matrix_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/mathys_2019/skinnider_2021/controls/barcodes_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/mathys_2019/skinnider_2021/controls/genes_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/mathys_2019/skinnider_2021/controls/author_barcode_annotations_controls.csv"))

#########################################################################################################
############################################# Schirmer 2019 #############################################
#########################################################################################################

## Motor cortex

datinfo <- rbind(datinfo, data.frame(First_Author="Schirmer", Year=2019, PubMedID=31316211, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="Motor/pre-motor cortex", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V2", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/motor/matrix_motor_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/motor/barcodes_motor_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/motor/genes_motor_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/motor/author_barcode_annotations_motor_controls.csv"))

## Parietal cortex

datinfo <- rbind(datinfo, data.frame(First_Author="Schirmer", Year=2019, PubMedID=31316211, Disease="Normal", Region="Cortex", Subregion="Parietal cortex", Author_Region="Parietal cortex", Region_Code="PC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V2", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/parietal/matrix_parietal_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/parietal/barcodes_parietal_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/parietal/genes_parietal_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/parietal/author_barcode_annotations_parietal_controls.csv"))

## PFC

datinfo <- rbind(datinfo, data.frame(First_Author="Schirmer", Year=2019, PubMedID=31316211, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="Prefrontal cortex", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V2", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/PFC/matrix_PFC_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/PFC/barcodes_PFC_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/PFC/genes_PFC_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/UCSC/PFC/author_barcode_annotations_PFC_controls.csv"))

#####################################################################################################
############################################# Tran 2020 #############################################
#####################################################################################################

## Hippocampus

datinfo <- rbind(datinfo, data.frame(First_Author="Tran", Year=2020, PubMedID=34582785, Disease="Normal", Region="Hippocampus", Subregion="Hippocampus", Author_Region="Hippocampus", Region_Code="HIP", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/hippocampus/matrix_hippocampus.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/hippocampus/barcodes_hippocampus.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/hippocampus/genes_hippocampus.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/hippocampus/author_barcode_annotations_hippocampus.csv"))

## DLPFC

datinfo <- rbind(datinfo, data.frame(First_Author="Tran", Year=2020, PubMedID=34582785, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="Dorsolateral prefrontal cortex", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/DLPFC/matrix_DLPFC.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/DLPFC/barcodes_DLPFC.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/DLPFC/genes_DLPFC.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/DLPFC/author_barcode_annotations_DLPFC.csv"))

## ACC

datinfo <- rbind(datinfo, data.frame(First_Author="Tran", Year=2020, PubMedID=34582785, Disease="Normal", Region="Cortex", Subregion="Anterior cingulate cortex", Author_Region="Subgenual anterior cingulate cortex", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/ACC/matrix_ACC.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/ACC/barcodes_ACC.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/ACC/genes_ACC.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/ACC/author_barcode_annotations_ACC.csv"))

## Amygdala

datinfo <- rbind(datinfo, data.frame(First_Author="Tran", Year=2020, PubMedID=34582785, Disease="Normal", Region="Amygdala", Subregion="Amygdala", Author_Region="Amygdala", Region_Code="AMY", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/amygdala/matrix_amygdala.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/amygdala/barcodes_amygdala.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/amygdala/genes_amygdala.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/amygdala/author_barcode_annotations_amygdala.csv"))

## NAc

datinfo <- rbind(datinfo, data.frame(First_Author="Tran", Year=2020, PubMedID=34582785, Disease="Normal", Region="Basal ganglia", Subregion="Striatum", Author_Region="Nucleus accumbens", Region_Code="STR", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/NAc/matrix_NAc.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/NAc/barcodes_NAc.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/NAc/genes_NAc.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/tran_2020/NAc/author_barcode_annotations_NAc.csv"))

##########################################################################################################
############################################# Velmeshev 2019 #############################################
##########################################################################################################

## PFC

datinfo <- rbind(datinfo, data.frame(First_Author="Velmeshev", Year=2019, PubMedID=31097668, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="Prefrontal cortex", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/PFC/matrix_PFC_adult_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/PFC/barcodes_PFC_adult_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/PFC/genes_PFC_adult_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/PFC/author_barcode_annotations_PFC_adult_controls.csv"))

## ACC

datinfo <- rbind(datinfo, data.frame(First_Author="Velmeshev", Year=2019, PubMedID=31097668, Disease="Normal", Region="Cortex", Subregion="Anterior cingulate cortex", Author_Region="BA24", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/ACC/matrix_ACC_adult_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/ACC/barcodes_ACC_adult_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/ACC/genes_ACC_adult_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/ACC/author_barcode_annotations_ACC_adult_controls.csv"))

######################################################################################################
############################################# Yang COVID #############################################
######################################################################################################

## Choroid plexus

datinfo <- rbind(datinfo, data.frame(First_Author="Yang", Year=2021, PubMedID=34153974, Disease="Normal", Region="Choroid plexus", Subregion="Choroid plexus", Author_Region="Choroid plexus", Region_Code="CP", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/Stanford/choroid/matrix_choroid_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/Stanford/choroid/barcodes_choroid_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/Stanford/choroid/genes_choroid_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/Stanford/choroid/author_barcode_annotations_choroid_controls.csv"))

## FCX

datinfo <- rbind(datinfo, data.frame(First_Author="Yang", Year=2021, PubMedID=34153974, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="Medial frontal gyrus", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/Stanford/PFC/matrix_PFC_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/Stanford/PFC/barcodes_PFC_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/Stanford/PFC/genes_PFC_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_COVID/Stanford/PFC/author_barcode_annotations_PFC_controls.csv"))

###################################################################################################
############################################# Yang AD #############################################
###################################################################################################

## FCX

datinfo <- rbind(datinfo, data.frame(First_Author="Yang", Year=2021, PubMedID=35165441, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="Superior frontal cortex", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/FCX/matrix_FCX_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/FCX/barcodes_FCX_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/FCX/genes_FCX_controls.tsv", Author_Barcode_Annotations=NA))

## Hippocampus

datinfo <- rbind(datinfo, data.frame(First_Author="Yang", Year=2021, PubMedID=35165441, Disease="Normal", Region="Hippocampus", Subregion="Hippocampus", Author_Region="Hippocampus", Region_Code="HIP", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/hippocampus/matrix_hippocampus_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/hippocampus/barcodes_hippocampus_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2021_AD/hippocampus/genes_hippocampus_controls.tsv", Author_Barcode_Annotations=NA))

#####################################################################################################
############################################# Nagy 2020 #############################################
#####################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Nagy", Year=2020, PubMedID=32341540, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="BA9", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V2", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/nagy_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/nagy_2020/controls/matrix_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/nagy_2020/controls/barcodes_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/nagy_2020/controls/genes_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/nagy_2020/controls/author_barcode_annotations_controls.csv"))

########################################################################################################
############################################# Grubman 2019 #############################################
########################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Grubman", Year=2019, PubMedID=31768052, Disease="Normal", Region="Cortex", Subregion="Temporal lobe", Author_Region="Entorhinal cortex", Region_Code="TC", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/grubman_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/grubman_2019/controls/expression_controls.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/grubman_2019/controls/barcodes_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/grubman_2019/controls/genes_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/grubman_2019/controls/author_barcode_annotations_controls.csv"))

#####################################################################################################
############################################# Franjic 2022 #############################################
#####################################################################################################

## CA1

datinfo <- rbind(datinfo, data.frame(First_Author="Franjic", Year=2022, PubMedID=34798047, Disease="Normal", Region="Hippocampus", Subregion="CA1", Author_Region="CA1", Region_Code="HIP", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/CA1/matrix_CA1.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/CA1/barcodes_CA1.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/CA1/genes_CA1.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/CA1/author_barcode_annotations_CA1.csv"))

## CA2-4

datinfo <- rbind(datinfo, data.frame(First_Author="Franjic", Year=2022, PubMedID=34798047, Disease="Normal", Region="Hippocampus", Subregion="CA2-4", Author_Region="CA2-4", Region_Code="HIP", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/CA2-4/matrix_CA24.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/CA2-4/barcodes_CA24.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/CA2-4/genes_CA24.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/CA2-4/author_barcode_annotations_CA24.csv"))

## DG

datinfo <- rbind(datinfo, data.frame(First_Author="Franjic", Year=2022, PubMedID=34798047, Disease="Normal", Region="Hippocampus", Subregion="Dentate gyrus", Author_Region="Dentate gyrus", Region_Code="HIP", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/DG/matrix_DG.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/DG/barcodes_DG.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/DG/genes_DG.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/DG/author_barcode_annotations_DG.csv"))

## EC

datinfo <- rbind(datinfo, data.frame(First_Author="Franjic", Year=2022, PubMedID=34798047, Disease="Normal", Region="Cortex", Subregion="Temporal lobe", Author_Region="Entorhinal cortex", Region_Code="TC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/EC/matrix_EC.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/EC/barcodes_EC.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/EC/genes_EC.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/EC/author_barcode_annotations_EC.csv"))

## Subiculum

datinfo <- rbind(datinfo, data.frame(First_Author="Franjic", Year=2022, PubMedID=34798047, Disease="Normal", Region="Subiculum", Subregion="Subiculum", Author_Region="Subiculum", Region_Code="HIP", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/subiculum/matrix_subiculum.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/subiculum/barcodes_subiculum.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/subiculum/genes_subiculum.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/franjic_2022/subiculum/author_barcode_annotations_subiculum.csv"))

#########################################################################################################
############################################# Badanjak 2020 #############################################
#########################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Badanjak", Year=2020, PubMedID=34805149, Disease="Normal", Region="Midbrain", Subregion="Midbrain", Author_Region="Ventral midbrain", Region_Code="MID", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/badanjak_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/badanjak_2020/controls/expression_controls.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/badanjak_2020/controls/barcodes_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/badanjak_2020/controls/genes_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/badanjak_2020/controls/author_barcode_annotations_controls.csv"))

########################################################################################################
############################################# Fullard 2021 #############################################
########################################################################################################

## DLPFC

datinfo <- rbind(datinfo, data.frame(First_Author="Fullard", Year=2021, PubMedID=34281603, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="Dorsolateral prefrontal cortex", Region_Code="FC", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/DLPFC/matrix_PFC_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/DLPFC/barcodes_DLPFC_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/DLPFC/genes_DLPFC_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/DLPFC/author_barcode_annotations_DLPFC.csv"))

## Choroid

datinfo <- rbind(datinfo, data.frame(First_Author="Fullard", Year=2021, PubMedID=34281603, Disease="Normal", Region="Choroid plexus", Subregion="Choroid plexus", Author_Region="Choroid plexus", Region_Code="CP", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/choroid/matrix_choroid_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/choroid/barcodes_choroid_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/choroid/genes_choroid_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/choroid/author_barcode_annotations_choroid.csv"))

## Medulla

datinfo <- rbind(datinfo, data.frame(First_Author="Fullard", Year=2021, PubMedID=34281603, Disease="Normal", Region="Brainstem", Subregion="Medulla oblongata", Author_Region="Medulla oblongata", Region_Code="MED", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/medulla/matrix_medulla_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/medulla/barcodes_medulla_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/medulla/genes_medulla_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/fullard_2021/medulla/author_barcode_annotations_medulla.csv"))

########################################################################################################
############################################# Absinta 2021 #############################################
########################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Absinta", Year=2021, PubMedID=34497421, Disease="Normal", Region="White matter", Subregion="White matter", Author_Region="White matter", Region_Code="WM", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/controls/expression_controls.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/controls/barcodes_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/controls/genes_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/absinta_2021/controls/author_barcode_annotations_controls.csv"))

#####################################################################################################
############################################# Lau 2020 #############################################
#####################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Lau", Year=2020, PubMedID=32989152, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="BA6 | BA8 | BA9", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/lau_2020/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/lau_2020/controls/matrix_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/lau_2020/controls/barcodes_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/lau_2020/controls/genes_controls.tsv", Author_Barcode_Annotations=NA))

#########################################################################################################
############################################# Morabito 2021 #############################################
#########################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Morabito", Year=2021, PubMedID=34239132, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="Prefrontal cortex", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/morabito_2021/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/morabito_2021/controls/matrix_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/morabito_2021/controls/barcodes_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/morabito_2021/controls/genes_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/morabito_2021/controls/author_barcode_annotations_controls.csv"))

####################################################################################################
############################################# Luo 2019 #############################################
####################################################################################################

## snmC2T-seq BA10

datinfo <- rbind(datinfo, data.frame(First_Author="Luo", Year=2019, PubMedID=35419551, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="BA10", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="snmCAT-seq", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/snmC2T/matrix.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/snmC2T/barcodes.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/snmC2T/genes.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/snmC2T/author_barcode_annotations.csv"))

## 10x BA44-45

datinfo <- rbind(datinfo, data.frame(First_Author="Luo", Year=2019, PubMedID=35419551, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="BA44-45", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/10x/BA44-45/matrix_BA44.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/10x/BA44-45/barcodes_BA44.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/10x/BA44-45/genes_BA44.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/10x/BA44-45/author_barcode_annotations_BA44.csv"))

## 10x BA46

datinfo <- rbind(datinfo, data.frame(First_Author="Luo", Year=2019, PubMedID=35419551, Disease="Normal", Region="Cortex", Subregion="Frontal cortex", Author_Region="BA46", Region_Code="FC", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/10x/BA46/matrix_BA46.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/10x/BA46/barcodes_BA46.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/10x/BA46/genes_BA46.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/luo_2019/10x/BA46/author_barcode_annotations_BA46.csv"))

#####################################################################################################
############################################# Yang 2022 #############################################
#####################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Yang", Year=2022, PubMedID=35349784, Disease="Normal", Region="Pons", Subregion="Pons", Author_Region="Trigeminal ganglion", Region_Code="TRI", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2022/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2022/matrix.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2022/barcodes.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2022/genes.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/yang_2022/author_barcode_annotations.csv"))

######################################################################################################
############################################# Jakel 2019 #############################################
######################################################################################################

datinfo <- rbind(datinfo, data.frame(First_Author="Jakel", Year=2022, PubMedID=30747918, Disease="Normal", Region="White matter", Subregion="White matter", Author_Region="White matter", Region_Code="WM", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/controls/expression_controls.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/controls/barcodes_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/controls/genes_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/jakel_2019/controls/author_barcode_annotations_controls.csv"))

#######################################################################################################
############################################# Kamath 2021 #############################################
#######################################################################################################

## Caudate nucleus

datinfo <- rbind(datinfo, data.frame(First_Author="Kamath", Year=2021, PubMedID=35513515, Disease="Normal", Region="Basal ganglia", Subregion="Striatum", Author_Region="Caudate nucleus", Region_Code="STR", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/CN/matrix_CN_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/CN/barcodes_CN_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/CN/genes_CN_controls.tsv", Author_Barcode_Annotations=NA))

## SN

datinfo <- rbind(datinfo, data.frame(First_Author="Kamath", Year=2021, PubMedID=35513515, Disease="Normal", Region="Basal ganglia", Subregion="Subsantia nigra", Author_Region="Substantia nigra pars compacta", Region_Code="MID", FACS_Sorted="Y", Unbiased_Sampling="N", Platform="10x Chromium V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/SN/matrix_SN_controls.mtx", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/SN/barcodes_SN_controls.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/SN/genes_SN_controls.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/kamath_2021/SN/author_barcode_annotations_SN_controls.csv"))

######################################################################################################
############################################# Ayhan 2021 #############################################
######################################################################################################

## Anterior hippocampus

datinfo <- rbind(datinfo, data.frame(First_Author="Ayhan", Year=2021, PubMedID=34051145, Disease="Epilepsy", Region="Hippocampus", Subregion="Hippocampus", Author_Region="Anterior hippocampus", Region_Code="HIP", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V2/V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/anterior_hippocampus/expression_anterior_hippocampus.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/anterior_hippocampus/barcodes_anterior_hippocampus.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/anterior_hippocampus/genes_anterior_hippocampus.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/anterior_hippocampus/author_barcode_annotations_anterior_hippocampus.csv"))

## Posterior hippocampus

datinfo <- rbind(datinfo, data.frame(First_Author="Ayhan", Year=2021, PubMedID=34051145, Disease="Epilepsy", Region="Hippocampus", Subregion="Hippocampus", Author_Region="Posterior hippocampus", Region_Code="HIP", FACS_Sorted="N", Unbiased_Sampling="Y", Platform="10x Chromium V2/V3", Directory="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/", Author_Counts="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/posterior_hippocampus/expression_posterior_hippocampus.csv", Author_Barcodes="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/posterior_hippocampus/barcodes_posterior_hippocampus.tsv", Author_Genes="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/posterior_hippocampus/genes_posterior_hippocampus.tsv", Author_Barcode_Annotations="/home/shared/scsn.expr_data/human_expr/postnatal/ayhan_2021/posterior_hippocampus/author_barcode_annotations_posterior_hippocampus.csv"))

############################################# Save #############################################

datinfo$Plot_Label <- paste(datinfo$First_Author, datinfo$Region_Code)
datinfo$Plot_Label[grepl("^BA", datinfo$Author_Region)&nchar(datinfo$Author_Region)<5] <- paste(datinfo$Plot_Label[grepl("^BA", datinfo$Author_Region)&nchar(datinfo$Author_Region)<5], datinfo$Author_Region[grepl("^BA", datinfo$Author_Region)&nchar(datinfo$Author_Region)<5])
datinfo$Plot_Label[grep("Ento", datinfo$Author_Region)] <- paste(datinfo$Plot_Label[grep("Ento", datinfo$Author_Region)], "EC")
datinfo$Plot_Label[grep("HIP", datinfo$Region_Code)] <- paste(datinfo$Plot_Label[grep("HIP", datinfo$Region_Code)], datinfo$Author_Region[grep("HIP", datinfo$Region_Code)])
datinfo$Plot_Label <- gsub(" | ", " ", datinfo$Plot_Label, fixed=T)
datinfo$Plot_Label <- gsub("Subiculum", "SUB", datinfo$Plot_Label)
datinfo$Plot_Label <- gsub("Dentate gyrus", "DG", datinfo$Plot_Label)
datinfo$Plot_Label <- gsub(" Hippocampus", "", datinfo$Plot_Label)
datinfo$Plot_Label <- gsub("Anterior hippocampus", "ANT", datinfo$Plot_Label)
datinfo$Plot_Label <- gsub("Posterior hippocampus", "POST", datinfo$Plot_Label)
datinfo$Plot_Label[grep("STR", datinfo$Region_Code)] <- paste(datinfo$Plot_Label[grep("STR", datinfo$Region_Code)], datinfo$Author_Region[grep("STR", datinfo$Region_Code)])
datinfo$Plot_Label <- gsub("Nucleus accumbens", "NAc", datinfo$Plot_Label)
datinfo$Plot_Label <- gsub("Putamen", "PUT", datinfo$Plot_Label)
datinfo$Plot_Label <- gsub("Caudate nucleus", "CN", datinfo$Plot_Label)
datinfo$Plot_Label[intersect(grep("Schirmer", datinfo$First_Author), grep("Motor", datinfo$Author_Region))] <- 
  gsub("FC", "FC M1", datinfo$Plot_Label[intersect(grep("Schirmer", datinfo$First_Author), grep("Motor", datinfo$Author_Region))])
datinfo$Plot_Label[intersect(grep("Schirmer", datinfo$First_Author), grep("Prefrontal", datinfo$Author_Region))] <- 
  gsub("FC", "PFC", datinfo$Plot_Label[intersect(grep("Schirmer", datinfo$First_Author), grep("Prefrontal", datinfo$Author_Region))])
datinfo$Plot_Label[intersect(grep("Bakken", datinfo$First_Author), grep("10x", datinfo$Platform))] <- 
  gsub("FC", "FC 10x", datinfo$Plot_Label[intersect(grep("Bakken", datinfo$First_Author), grep("10x", datinfo$Platform))])
datinfo$Plot_Label[intersect(grep("Bakken", datinfo$First_Author), grep("SMART", datinfo$Platform))] <- 
  gsub("FC", "FC SSv4", datinfo$Plot_Label[intersect(grep("Bakken", datinfo$First_Author), grep("SMART", datinfo$Platform))])
datinfo$Plot_Label[intersect(grep("Tran", datinfo$First_Author), grep("Dorso", datinfo$Author_Region))] <- 
  gsub("FC", "FC DLPFC", datinfo$Plot_Label[intersect(grep("Tran", datinfo$First_Author), grep("Dorso", datinfo$Author_Region))])
datinfo$Plot_Label[intersect(grep("Tran", datinfo$First_Author), grep("cingulate", datinfo$Author_Region))] <- 
  gsub("FC", "FC ACC", datinfo$Plot_Label[intersect(grep("Tran", datinfo$First_Author), grep("cingulate", datinfo$Author_Region))])
datinfo$Plot_Label[intersect(grep("Yang", datinfo$First_Author), grep("Medial", datinfo$Author_Region))] <- 
  gsub("FC", "FC MED", datinfo$Plot_Label[intersect(grep("Yang", datinfo$First_Author), grep("Medial", datinfo$Author_Region))])
datinfo$Plot_Label[intersect(grep("Yang", datinfo$First_Author), grep("Superior", datinfo$Author_Region))] <- 
  gsub("FC", "FC SUP", datinfo$Plot_Label[intersect(grep("Yang", datinfo$First_Author), grep("Superior", datinfo$Author_Region))])
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Velmeshev FC")] <- "Velmeshev PFC"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Morabito FC")] <- "Morabito PFC"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Schirmer M1 FC")] <- "Schirmer FC M1"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Luo FC")] <- "Luo FC BA44-45"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Hodge TC")] <- "Hodge TC MTG"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Habib FC")] <- "Habib PFC"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Fullard FC")] <- "Fullard DLPFC"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Bakken FC 10x")] <- "Bakken FC M1 10x"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Bakken FC SSv4")] <- "Bakken FC M1 SSv4"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "Agarwal FC")] <- "Agarwal FC MFG"
datinfo$Plot_Label[is.element(datinfo$Plot_Label, "ABI FC")] <- "ABI FC ACC"

datinfo$Dataset <- gsub(" ", "_", paste(datinfo$First_Author, datinfo$Year, datinfo$Platform, datinfo$Author_Region, sep="_"))
datinfo$Dataset <- gsub("/", "_", datinfo$Dataset)
datinfo$Dataset <- gsub("_|_", "_", datinfo$Dataset, fixed=T)

datinfo$Study <- paste(datinfo$First_Author, datinfo$Year)
datinfo$Study[is.element(datinfo$PubMedID, "35165441")] <- paste(datinfo$Study[is.element(datinfo$PubMedID, "35165441")], " ")

## Add metadata about the sampling strategies used for each dataset/study:

datinfo$Filter_Size <- NA

datinfo$Filter_Size[grep("Absinta", datinfo$Dataset)] <- 40
datinfo$Filter_Size[grep("Agarwal", datinfo$Dataset)] <- 35
datinfo$Filter_Size[grep("Franjic", datinfo$Dataset)] <- 35
datinfo$Filter_Size[grep("Grubman", datinfo$Dataset)] <- 40
datinfo$Filter_Size[grep("Habib", datinfo$Dataset)] <- 35
datinfo$Filter_Size[grep("Hodge", datinfo$Dataset)] <- 30
datinfo$Filter_Size[grep("Kamath", datinfo$Dataset)] <- 70 
datinfo$Filter_Size[grep("Khrameeva", datinfo$Dataset)] <- 30
datinfo$Filter_Size[grep("Lake", datinfo$Dataset)] <- 50
datinfo$Filter_Size[grep("Lee", datinfo$Dataset)] <- 40
datinfo$Filter_Size[grep("Li", datinfo$Dataset)] <- 40
datinfo$Filter_Size[grep("Mathys", datinfo$Dataset)] <- 40
datinfo$Filter_Size[grep("Morabito", datinfo$Dataset)] <- 70
datinfo$Filter_Size[grep("Nagy", datinfo$Dataset)] <- 30
datinfo$Filter_Size[grep("Velmeshev", datinfo$Dataset)] <- 30
datinfo$Filter_Size[grep("Yang", datinfo$Dataset)] <- 40

datinfo$Mean_PMI <- NA

datinfo$Mean_PMI[grep("Absinta", datinfo$Dataset)] <- 7
datinfo$Mean_PMI[intersect(grep("Agarwal", datinfo$Dataset), grep("FC", datinfo$Region_Code))] <- 34.7
datinfo$Mean_PMI[intersect(grep("Agarwal", datinfo$Dataset), grep("MID", datinfo$Region_Code))] <- 33.2
datinfo$Mean_PMI[grep("ABI", datinfo$Dataset)] <- 22.7
datinfo$Mean_PMI[grep("Badanjak", datinfo$Dataset)] <- 15.3
datinfo$Mean_PMI[grep("Badanjak", datinfo$Dataset)] <- 15.3
datinfo$Mean_PMI[intersect(grep("Bakken", datinfo$Dataset), grep("10x", datinfo$Platform))] <- 14
datinfo$Mean_PMI[intersect(grep("Bakken", datinfo$Dataset), grep("SMART", datinfo$Platform))] <- 24
datinfo$Mean_PMI[grep("Franjic", datinfo$Dataset)] <- 12.7
datinfo$Mean_PMI[intersect(grep("Franjic", datinfo$Dataset), grep("Dentate", datinfo$Author_Region))] <- 15.7
datinfo$Mean_PMI[grep("Fullard", datinfo$Dataset)] <- 15.1
datinfo$Mean_PMI[grep("Grubman", datinfo$Dataset)] <- 36.7
datinfo$Mean_PMI[intersect(grep("Habib", datinfo$Dataset), grep("FC", datinfo$Region_Code))] <- 10.4
datinfo$Mean_PMI[intersect(grep("Habib", datinfo$Dataset), grep("HIP", datinfo$Region_Code))] <- 11.9
datinfo$Mean_PMI[grep("Hodge", datinfo$Dataset)] <- 22.5
datinfo$Mean_PMI[intersect(grep("Kamath", datinfo$Dataset), grep("MID", datinfo$Region_Code))] <- 16.5
datinfo$Mean_PMI[intersect(grep("Kamath", datinfo$Dataset), grep("STR", datinfo$Region_Code))] <- 17.9
datinfo$Mean_PMI[intersect(grep("Lake", datinfo$Dataset), grep("CB", datinfo$Region_Code))] <- 17.3
datinfo$Mean_PMI[intersect(grep("Lake", datinfo$Dataset), grep("BA6", datinfo$Author_Region))] <- 14
datinfo$Mean_PMI[intersect(grep("Lake", datinfo$Dataset), grep("BA10", datinfo$Author_Region))] <- 17.5
datinfo$Mean_PMI[intersect(grep("Lake", datinfo$Dataset), grep("OC", datinfo$Region_Code))] <- 16
datinfo$Mean_PMI[grep("Lau", datinfo$Dataset)] <- 28.7
datinfo$Mean_PMI[intersect(grep("Lee", datinfo$Dataset), grep("Caudate", datinfo$Author_Region))] <- 11.5
datinfo$Mean_PMI[intersect(grep("Lee", datinfo$Dataset), grep("Putamen", datinfo$Author_Region))] <- 11.6
datinfo$Mean_PMI[grep("Li", datinfo$Dataset)] <- 13
datinfo$Mean_PMI[intersect(grep("Luo", datinfo$Dataset), grep("BA10", datinfo$Author_Region))] <- 13.5
datinfo$Mean_PMI[intersect(grep("Luo", datinfo$Dataset), grep("BA4", datinfo$Author_Region))] <- 12
datinfo$Mean_PMI[grep("Morabito", datinfo$Dataset)] <- 4.7
datinfo$Mean_PMI[grep("Nagy", datinfo$Dataset)] <- 34
datinfo$Mean_PMI[intersect(grep("Schirmer", datinfo$Dataset), grep("Motor", datinfo$Author_Region))] <- 19.5
datinfo$Mean_PMI[intersect(grep("Schirmer", datinfo$Dataset), grep("Pre", datinfo$Author_Region))] <- 22.5
datinfo$Mean_PMI[intersect(grep("Schirmer", datinfo$Dataset), grep("Par", datinfo$Author_Region))] <- 22
datinfo$Mean_PMI[intersect(grep("Tran", datinfo$Dataset), grep("AMY", datinfo$Region_Code))] <- 26
datinfo$Mean_PMI[intersect(grep("Tran", datinfo$Dataset), grep("cingulate", datinfo$Author_Region))] <- 26
datinfo$Mean_PMI[intersect(grep("Tran", datinfo$Dataset), grep("Dorso", datinfo$Author_Region))] <- 29.7
datinfo$Mean_PMI[intersect(grep("Tran", datinfo$Dataset), grep("HIP", datinfo$Region_Code))] <- 29.3
datinfo$Mean_PMI[intersect(grep("Tran", datinfo$Dataset), grep("STR", datinfo$Region_Code))] <- 26.9
datinfo$Mean_PMI[intersect(grep("Velmeshev", datinfo$Dataset), grep("BA24", datinfo$Author_Region))] <- 17.5
datinfo$Mean_PMI[intersect(grep("Velmeshev", datinfo$Dataset), grep("Pre", datinfo$Author_Region))] <- 18.9
datinfo$Mean_PMI[intersect(grep("Yang", datinfo$Dataset), grep("Superior", datinfo$Author_Region))] <- 5.8
datinfo$Mean_PMI[intersect(grep("Yang", datinfo$Dataset), grep("HIP", datinfo$Region_Code))] <- 5.7

#fwrite(datinfo, file="/home/rebecca/SCSN_meta_analysis/prelim_new/3_init_datinfo/datinfo_SCSN_meta_analysis_init.csv")


