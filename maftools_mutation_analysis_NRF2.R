
library(maftools)
library(stringr)
library(here)
#my_maf = read.maf(maf = "TCGA.HNSC.varscan.5296cf00-4d8c-4db3-80d7-930a4b44f90d.DR-10.0.somatic.maf", clinicalData = NULL)

#GISTIC HNSC
#all.lesions = here("gistic", "all_lesions.conf_99.txt")
#amp.genes = here("gistic", "amp_genes.conf_99.txt")
#del.genes = here("gistic", "del_genes.conf_99.txt")
#scores.gis = here("gistic", "scores.gistic")

#Read GISTIC results along with MAF
#my_maf = read.maf(
#  maf = "TCGA.HNSC.varscan.5296cf00-4d8c-4db3-80d7-930a4b44f90d.DR-10.0.somatic.maf",
#  gisticAllLesionsFile = all.lesions,
#  gisticAmpGenesFile = amp.genes,
#  gisticDelGenesFile = del.genes,
#  gisticScoresFile = scores.gis,
#  isTCGA = TRUE,
#  verbose = FALSE, 
#  clinicalData = NULL)

load("merged_clinical_OC.RData")
load("cn_pway_data_NRF2_GEMM.RData")
holder <- names(cn_pway_data)[1:16]
values <- c()
sample_ids <- c()
cn_genes <- c()
for(i in holder){
  values <- c(values, cn_pway_data[[i]])
  sample_ids <- c(sample_ids, cn_pway_data$sample_id)
  for(ii in 1:length(cn_pway_data$sample_id)){
    cn_genes <- c(cn_genes, str_sub(i, start = 1L, end = -7L))
  }
}

custom_cn_table_raw <- data.frame("Gene" = cn_genes, "Sample_Name" = sample_ids, "Values" = values )
custom_cn_table_raw <- subset(custom_cn_table_raw, custom_cn_table_raw$Values %in% c(2, -2))

amp_del <- c()
for(i in custom_cn_table_raw$Values){
  if(i == 2){
    amp_del <- c(amp_del, "Amp")
  }else if(i == -2){
    amp_del <- c(amp_del, "Del")
  }
}


custom.cn.data <- data.frame("Gene" = custom_cn_table_raw$Gene, "Sample_Name" = custom_cn_table_raw$Sample_Name, "CN" = amp_del)

my_maf = read.maf(maf = "TCGA.HNSC.varscan.5296cf00-4d8c-4db3-80d7-930a4b44f90d.DR-10.0.somatic.maf", clinicalData = NULL, cnTable = custom.cn.data)
my_maf@data$Tumor_Sample_Barcode <- str_sub(my_maf@data$Tumor_Sample_Barcode, start = 1L, end = 12L)

my_maf_mut_only = read.maf(maf = "TCGA.HNSC.varscan.5296cf00-4d8c-4db3-80d7-930a4b44f90d.DR-10.0.somatic.maf", clinicalData = NULL)
my_maf_mut_only@data$Tumor_Sample_Barcode <- str_sub(my_maf_mut_only@data$Tumor_Sample_Barcode, start = 1L, end = 12L)

#oad("clin_fact_merged_HPVpos_inclusive_TCGA.RData")

clinical_factors_merged <- merged_clinical_OC

clin_lt_zero <- subset(clinical_factors_merged, clinical_factors_merged$pway_active_cent)
clin_gt_zero <- subset(clinical_factors_merged, !(clinical_factors_merged$pway_active_cent))

oc_maf <- subsetMaf(maf = my_maf, tsb = as.character(clinical_factors_merged$sample_id) , mafObj = TRUE)
oc_maf@data <- droplevels(oc_maf@data)


oc_maf_mut_only <- subsetMaf(maf = my_maf_mut_only, tsb = as.character(clinical_factors_merged$sample_id) , mafObj = TRUE)
oc_maf_mut_only@data <- droplevels(oc_maf_mut_only@data)

#load("TCGA_NFkB_PCA_data.RData")

#TCGA_NFkB_PCA_data$sample_id <- gsub(TCGA_NFkB_PCA_data$sample, pattern = "[.]", replacement = "-")

#NFkB_active <- subset(TCGA_NFkB_PCA_data, TCGA_NFkB_PCA_data$NFkB_mod_PCA > median(TCGA_NFkB_PCA_data$NFkB_mod_PCA))
#NFkB_inactive <- subset(TCGA_NFkB_PCA_data, TCGA_NFkB_PCA_data$NFkB_mod_PCA <= median(TCGA_NFkB_PCA_data$NFkB_mod_PCA) )

nrf2_active_maf <- subsetMaf(maf = oc_maf, tsb =  as.character(clin_lt_zero$sample_id) , mafObj = TRUE)
#nfkb_active_maf <- read.maf(data.frame(nfkb_active_maf@data))
nrf2_inactive_maf <- subsetMaf(maf = oc_maf, tsb = as.character(clin_gt_zero$sample_id) , mafObj = TRUE)
#nfkb_inactive_maf <- read.maf(data.frame(nfkb_inactive_maf@data))

nrf2_active_maf_mut_only <- subsetMaf(maf = oc_maf_mut_only, tsb =  as.character(clin_lt_zero$sample_id) , mafObj = TRUE)
#nfkb_active_maf <- read.maf(data.frame(nfkb_active_maf@data))
nrf2_inactive_maf_mut_only <- subsetMaf(maf = oc_maf_mut_only, tsb = as.character(clin_gt_zero$sample_id) , mafObj = TRUE)
#nfkb_inactive_maf <- read.maf(data.frame(nfkb_inactive_maf@data))

gene_summary_df <- getGeneSummary(nrf2_active_maf)
plot_me_act <- gene_summary_df$Hugo_Symbol[1:10]
gene_summary_df <- getGeneSummary(nrf2_inactive_maf)
plot_me_inact <- c(plot_me, gene_summary_df$Hugo_Symbol[1:10])

plot_me <- c("NFE2L2", "KEAP1", "CUL3", "CDKN2A", "TP53")
#plot_me <- c(plot_me_act, plot_me_inact)
library(RColorBrewer)
col = RColorBrewer::brewer.pal(n = 9, name = 'Set1')
col = c(col, "black")

names(col) = c('Amp','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del', 'Frame_Shift_Del', 'Del')
bg <- rgb(0,0,0, 0.0)


#quartz()
pdf(file = "coOncoplot_NFKB_MEDIAN_thresh_PCA_CNV_MEDIAN.pdf", height = 6, width = 12)
coOncoplot(nrf2_inactive_maf, nrf2_active_maf, genes = plot_me, m1Name = "NFkB Inactive (N=35)", m2Name = "NFkB Active (N=26)", colors = col, bgCol = bg, borderCol = bg, annotationFontSize = 8)
dev.off()


pdf(file = "Oncoplot_OC_NRF2_PWAY_GEMM.pdf", height = 6, width = 12)
oncoplot(oc_maf, genes = plot_me, colors = col, bgCol = bg, borderCol = bg, annotationFontSize = 8, cohortSize = 268 )
dev.off()



print(mafCompare(nrf2_inactive_maf, nrf2_active_maf, useCNV = FALSE, minMut = 2))

comparison_df <- mafCompare(nfkb_inactive_maf, nfkb_active_maf, m1Name = "NFkB_InacmakePie(active_signatures)tive", m2Name = "NFkB_Active", useCNV = FALSE, minMut = 2)

write.table(comparison_df$results, file = "gene_comparison_table_nfkb.txt", sep = "\t", quote = FALSE)

my.tnm = trinucleotideMatrix(maf = hpv_opscc_maf, ref_genome = 'BSgenome.Hsapiens.UCSC.hg38')

apobec_data <- my.tnm$APOBEC_scores

apobec_data$nfkb_active <- apobec_data$Tumor_Sample_Barcode %in% nfkb_active_maf@data$Tumor_Sample_Barcode

#math_data <- math.score(hpv_opscc_maf)
#apobec_data <- merge(math_data, apobec_data, by = "Tumor_Sample_Barcode")

apobec_data_nfkb_active <- subset(apobec_data, apobec_data$nfkb_active == TRUE)
apobec_data_nfkb_inactive <- subset(apobec_data, apobec_data$nfkb_active == FALSE)

#print(wilcox.test(apobec_data_nfkb_active$MATH, apobec_data_nfkb_inactive$MATH))
print(wilcox.test(apobec_data_nfkb_active$n_mutations, apobec_data_nfkb_inactive$n_mutations))
print(ks.test(apobec_data_nfkb_active$n_mutations, apobec_data_nfkb_inactive$n_mutations))
print(wilcox.test(apobec_data_nfkb_active$fraction_APOBEC_mutations, apobec_data_nfkb_inactive$fraction_APOBEC_mutations))
print(ks.test(apobec_data_nfkb_active$fraction_APOBEC_mutations, apobec_data_nfkb_inactive$fraction_APOBEC_mutations))

nmf_matrix <- my.tnm$nmf_matrix

nmf_matrix_nfkb_active <- nmf_matrix[row.names(nmf_matrix) %in% clin_lt_zero$sample_id, ]
nmf_matrix_nfkb_inactive <- nmf_matrix[row.names(nmf_matrix) %in% clin_gt_zero$sample_id, ]

mytest <- rbind(colSums(nmf_matrix_nfkb_active), colSums(nmf_matrix_nfkb_inactive))
mytest <- data.frame(mytest)
names(mytest) <- colnames(nmf_matrix)
row.names(mytest) <- c("active", "inactive")

library(deconstructSigs)
active_signatures = whichSignatures(tumor.ref = mytest, 
                                        signatures.ref = signatures.cosmic, 
                                        sample.id = "active" ,
                                        contexts.needed = TRUE,
                                        tri.counts.method = 'default')

plotSignatures(active_signatures)
makePie(active_signatures)



inactive_signatures = whichSignatures(tumor.ref = mytest, 
                                    signatures.ref = signatures.cosmic, 
                                    sample.id = "inactive" ,
                                    contexts.needed = TRUE,
                                    tri.counts.method = 'default')

plotSignatures(inactive_signatures)
makePie(inactive_signatures)

library(ggplot2)
library(gridExtra)


apobec_data$log_mut_burden <- log(apobec_data$n_mutations, base = 10)

box_burden <- ggplot(apobec_data, aes(y = log_mut_burden, x = nfkb_active, color = nfkb_active )) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.15))+
  scale_color_grey(start = 1-0.5, end = 1-0.95)+
  theme_classic()

box_apobec <- ggplot(apobec_data, aes(y = fraction_APOBEC_mutations, x = nfkb_active, color = nfkb_active )) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.15))+
  scale_color_grey(start = 1-0.5, end = 1-0.95)+
  theme_classic()

pdf(file = "box_plots_MEDIAN_thesh.pdf", width = 5.5, height = 3)
grid.arrange(box_burden, box_apobec, ncol = 2, nrow = 1)
dev.off()


  

hn_merged_alterations <- read.table(file = "merged_pway_alt_complete_prediction_hpv_neg_hnsc_edittps.csv", sep = ",", header = TRUE)

oc_merged_alterations <- subset(hn_merged_alterations, hn_merged_alterations$sample_id %in% clinical_factors_merged$sample_id)

write.table(file = "merged_oc_NRF2_pway_alterations.csv", oc_merged_alterations, sep = ",", row.names = FALSE)

