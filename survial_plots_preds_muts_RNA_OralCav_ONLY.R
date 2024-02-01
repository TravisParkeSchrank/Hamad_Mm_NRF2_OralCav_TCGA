library("edgeR")
library("limma")
library("Glimma")
library("Homo.sapiens")
library("cancerclass")
library("caret")
library("AdaSampling")
library("rlist")
library("survival")


library("stringr")
library("ggplot2")
library("ggExtra")
library("stringr")

seed = 8
criter = 0.3
train_partition = 0.33


set.seed(seed)

truncate <- function(x){
	return(substr(x, 1, 12))
}





load("norm_lcpm.RData")
load("pway_alt_data_complete.RData")
load("avg_all_ada_pred_df.RData")

pred_cent <- read.table("merged_pway_alt_complete_prediction_hnsc_PAN.csv", header = TRUE, sep = ",")

avg_all_ada_pred_df$sample_id <- gsub("[.]", "-", avg_all_ada_pred_df$Group.1)

merged <- merge(pred_cent, avg_all_ada_pred_df, by = "sample_id")

pway_alt_data_complete <- merged

no_pway_changes <- subset(pway_alt_data_complete, is.na(pway_alt_data_complete$cul3_position))
no_pway_changes <- subset(no_pway_changes, is.na(no_pway_changes$nfe2l2_position))
no_pway_changes <- subset(no_pway_changes, is.na(no_pway_changes$keap1_position))
no_pway_changes<- subset(no_pway_changes, no_pway_changes$NFE2L2_value < 2)
no_pway_changes<- subset(no_pway_changes, no_pway_changes$CUL3_value > -2)
no_pway_changes<- subset(no_pway_changes, no_pway_changes$KEAP1_value > -2)
no_pway_changes<- subset(no_pway_changes, no_pway_changes$splice_event == 0)
no_pway_changes_ids <- no_pway_changes$sample_id
#no_pway_changes_ids <- gsub("[.]", ".", c( as.character(no_pway_changes_ids)))
no_pway_changes_ids <- paste(no_pway_changes_ids, ".01", sep = "")
no_pway_changes_ids <- truncate(no_pway_changes_ids)

merged$no_pway_change <- merged$sample_id %in% no_pway_changes_ids

merged_no_pway_change <- subset(merged, merged$no_pway_change == TRUE)

clinical_xena <- read.table("Survival_SupplementalTable_S1_20171025_xena_sp.txt", sep = "\t", header = TRUE)
normals <- substring(as.character(clinical_xena$sample), 14, 15) %in% c("11", "06")
clinical_xena$noramals <- normals
clinical_xena$tumor <- !normals
clinical_xena <- subset(clinical_xena, clinical_xena$tumor == TRUE)
clinical_xena <- subset(clinical_xena, clinical_xena$cancer.type.abbreviation  == "HNSC")
clinical_xena$sample_id <- clinical_xena$X_PATIENT


merged_clinical <- merge(merged, clinical_xena, by = "sample_id")


#merged_clinical <- subset(merged_clinical, merged_clinical$clinical_stage %in% c("Stage I", "Stage II", "Stage III", "Stage IVA", "Stage IVB"))

merged_clinical$early_stage <- merged_clinical$clinical_stage %in% c("Stage I", "Stage II")
merged_clinical$advanced_stage <- !(merged_clinical$clinical_stage %in% c("Stage I", "Stage II"))


merged_clinical$zeta <- merged_clinical$zeta * -1

merged_clinical$pway_active_cent <- merged_clinical$zeta > -0.5

merged_clinical$pway_active_ada <- merged_clinical$ada_p_pos >= 0.3



firehose_clinical <- read.table(file = "hnsc_firehose_clinical.txt", header = TRUE, sep = "\t")

firehose_clinical$sample_id <- toupper(firehose_clinical$bcr_patient_barcode)

merged_clinical <- merge(merged_clinical, firehose_clinical, by = "sample_id")

#oral_cavity <- c("alveolar ridge", "buccal mucosa", "floor of mouth",  "hard palate", "oral cavity", "oral tongue")
oral_cavity <- c("alveolar ridge", "buccal mucosa", "floor of mouth",  "oral cavity", "oral tongue")
#oral_cavity <- c( "oral tongue")
merged_clinical <- subset(merged_clinical, merged_clinical$anatomic_neoplasm_subdivision %in% oral_cavity)

merged_clinical$current_smoker <- merged_clinical$tobacco_smoking_history %in% c(2)
merged_clinical$race_binary <- merged_clinical$race.x %in% c("BLACK OR AFRICAN AMERICAN")


print("centroid analysis")

CoxFit <- coxph(Surv(OS.time, OS) ~ pway_active_cent, data = merged_clinical);
print(summary(CoxFit));

CoxFit <- coxph(Surv(OS.time, OS) ~ pway_active_cent+advanced_stage+current_smoker+race_binary, data = merged_clinical);
print(summary(CoxFit));

CoxFit <- coxph(Surv(OS.time, OS) ~ current_smoker, data = merged_clinical);
print(summary(CoxFit));

CoxFit <- coxph(Surv(DSS.time, DSS) ~ pway_active_cent, data = merged_clinical);
print(summary(CoxFit));

CoxFit <- coxph(Surv(DSS.time, DSS) ~ pway_active_cent+advanced_stage+current_smoker+race_binary, data = merged_clinical);
print(summary(CoxFit));

CoxFit <- coxph(Surv(PFI.time, PFI) ~ pway_active_cent, data = merged_clinical);
print(summary(CoxFit));

CoxFit <- coxph(Surv(PFI.time, PFI) ~ pway_active_cent+advanced_stage+current_smoker+race_binary, data = merged_clinical);
print(summary(CoxFit));


SurvA <- merged_clinical[merged_clinical$pway_active_cent == TRUE,]$DSS.time;
SurvB <- merged_clinical[merged_clinical$pway_active_cent == TRUE,]$DSS;
Fit_7 <- survfit(Surv(SurvA, SurvB) ~ 1);

SurvA <- merged_clinical[merged_clinical$pway_active_cent == FALSE,]$DSS.time;
SurvB <- merged_clinical[merged_clinical$pway_active_cent == FALSE,]$DSS;
Fit_1 <- survfit(Surv(SurvA, SurvB) ~ 1);
    
print("Plot Univariate DSS")

print(summary(Fit_1, c(5*12, 10*12, 15*12)));
pdf(file = "DSS_HNSC_OCSCC_cent.pdf", width = 4, height = 4);
plot(Fit_1, conf.int=FALSE, ylim = c(0.0,1), xlim = c(0.0,5000), xlab = "Time(days)", ylab = "Prop. Surviving", cex.lab = 1.0, axes=FALSE, main = "");
	box(lwd = 1.5, cex = 1.5);
   axis(side = 1, lwd = 1.1, cex.axis = 1.6);
	axis(side = 2, lwd = 1.1, cex.axis = 1.6);

   #lines(Fit_2, conf.int=FALSE,lwd = 1.5, col="darkgrey", mark.time = T);
   lines(Fit_1, conf.int=FALSE,lwd = 1.5, col="black", mark.time = T);
   lines(Fit_7, conf.int=FALSE,lwd = 1.5, col="red", mark.time = T);
   #lines(Fit_8, conf.int=FALSE,lwd = 1.5,col="darkred", mark.time = T);	
	legend(40, 0.2, legend = c("NRF2 Inactive","NRF2 Active"), 
	bty = "n",
	cex = c(0.75),
	lty = c(1),
	col = c("black", "red")
	
	);
dev.off()


print("ada sampling analysis")

CoxFit <- coxph(Surv(OS.time, OS) ~ pway_active_ada, data = merged_clinical);
print(summary(CoxFit));

CoxFit <- coxph(Surv(DSS.time, DSS) ~ pway_active_ada, data = merged_clinical);
print(summary(CoxFit));

CoxFit <- coxph(Surv(PFI.time, PFI) ~ pway_active_ada, data = merged_clinical);
print(summary(CoxFit));

SurvA <- merged_clinical[merged_clinical$pway_active_ada == TRUE,]$DSS.time;
SurvB <- merged_clinical[merged_clinical$pway_active_ada == TRUE,]$DSS;
Fit_7 <- survfit(Surv(SurvA, SurvB) ~ 1);

SurvA <- merged_clinical[merged_clinical$pway_active_ada == FALSE,]$DSS.time;
SurvB <- merged_clinical[merged_clinical$pway_active_ada == FALSE,]$DSS;
Fit_1 <- survfit(Surv(SurvA, SurvB) ~ 1);
    
print("Plot Univariate DSS")

print(summary(Fit_1, c(5*12, 10*12, 15*12)));
pdf(file = "DSS_OCSCC_ES_ada.pdf", width = 4, height = 4);
plot(Fit_1, conf.int=FALSE, ylim = c(0.0,1), xlim = c(0.0,5000), xlab = "Time(days)", ylab = "Prop. Surviving", cex.lab = 1.0, axes=FALSE, main = "");
	box(lwd = 1.5, cex = 1.5);
   axis(side = 1, lwd = 1.1, cex.axis = 1.6);
	axis(side = 2, lwd = 1.1, cex.axis = 1.6);

   #lines(Fit_2, conf.int=FALSE,lwd = 1.5, col="darkgrey", mark.time = T);
   lines(Fit_1, conf.int=FALSE,lwd = 1.5, col="black", mark.time = T);
   lines(Fit_7, conf.int=FALSE,lwd = 1.5, col="red", mark.time = T);
   #lines(Fit_8, conf.int=FALSE,lwd = 1.5,col="darkred", mark.time = T);	
	legend(40, 0.2, legend = c("NRF2 Inactive","NRF2 Active"), 
	bty = "n",
	cex = c(0.75),
	lty = c(1),
	col = c("black", "red")
	
	);
dev.off()
    
    
    
   
p=ggplot(merged_clinical, aes(x=zeta, y=ada_p_pos, color = no_pway_change)) +
       geom_point() +
       scale_color_manual(values = c("blue", "red")) +
       theme(legend.position="left")
       
correlate_with_hist <- ggMarginal(p, type="histogram", col = "blue", fill = "blue")


correlate_with_density <- ggMarginal(p, type="density")







CoxFit <- coxph(Surv(DSS.time, DSS) ~ no_pway_change, data = merged_clinical);
print(summary(CoxFit));

library(survminer)
library(wesanderson)


fit <- survfit(Surv(OS.time, OS) ~ no_pway_change, data = merged_clinical)
print(summary(fit))
my_surv_plot_OS_MUTS <- ggsurvplot(
  fit, 
  data = merged_clinical, 
  size = 1,                 # change line size
  palette = c(wes_palette("Zissou1", 21, type = "continuous")[19],wes_palette("Zissou1", 21, type = "continuous")[3]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("NRF2act", "WT"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)


fit <- survfit(Surv(DSS.time, DSS) ~ no_pway_change, data = merged_clinical)
print(summary(fit))
my_surv_plot_DSS_MUTS <- ggsurvplot(
  fit, 
  data = merged_clinical, 
  size = 1,                 # change line size
  palette = c(wes_palette("Zissou1", 21, type = "continuous")[19],wes_palette("Zissou1", 21, type = "continuous")[3]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("NRF2act", "WT"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)


fit <- survfit(Surv(PFI.time, PFI) ~ no_pway_change, data = merged_clinical)
print(summary(fit))
my_surv_plot_PFI_MUTS <- ggsurvplot(
  fit, 
  data = merged_clinical, 
  size = 1,                 # change line size
  palette = c(wes_palette("Zissou1", 21, type = "continuous")[19],wes_palette("Zissou1", 21, type = "continuous")[3]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("NRF2act", "WT"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)


fit <- survfit(Surv(OS.time, OS) ~ pway_active_cent, data = merged_clinical)
print(summary(fit))
my_surv_plot_OS_CENT <- ggsurvplot(
  fit, 
  data = merged_clinical, 
  size = 1,                 # change line size
  palette = c(wes_palette("Zissou1", 21, type = "continuous")[3],wes_palette("Zissou1", 21, type = "continuous")[19]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("WT", "NRF2act"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)



fit <- survfit(Surv(PFI.time, PFI) ~ pway_active_cent, data = merged_clinical)
print(summary(fit))
my_surv_plot_PFI_CENT <- ggsurvplot(
  fit, 
  data = merged_clinical, 
  size = 1,                 # change line size
  palette = rev(c(wes_palette("Zissou1", 21, type = "continuous")[19],wes_palette("Zissou1", 21, type = "continuous")[3])),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("WT", "NRF2act"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)


fit <- survfit(Surv(DSS.time, DSS) ~ pway_active_cent, data = merged_clinical)
print(summary(fit))
my_surv_plot_DSS_CENT <- ggsurvplot(
  fit, 
  data = merged_clinical, 
  size = 1,                 # change line size
  palette = c(wes_palette("Zissou1", 21, type = "continuous")[3],wes_palette("Zissou1", 21, type = "continuous")[19]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("WT", "NRF2act"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

####quartz("test", 2.75, 4)

merged_clinical_OC <- merged_clinical
save(merged_clinical_OC, file = "merged_clinical_OC.RData")