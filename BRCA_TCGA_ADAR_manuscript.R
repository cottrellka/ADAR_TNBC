library(RTCGA)
library(Biobase)
library(affy)
library(genefilter)
library(XML)
library(reshape2)
library(survival)
library(survminer)
require("survival")
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(TCGAretriever)
library(extrafont)
library(DESeq2)
library(tidyverse)
library(edgeR)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_science <- function (base_size = 12, base_family = "Arial Rounded MT Bold") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}

setwd("/Users/cottr/Box Sync/TCGA/")



# read in RNAseq data for breast cancer samples from TCGA - 
RNAseq_BRCA <- fread("/Users/cottr/Box Sync/TCGA/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt")

#first row has relevant column names, make col.names match the first row
colnames(RNAseq_BRCA) <- paste(colnames(RNAseq_BRCA), RNAseq_BRCA[1,], sep= "-")

#hybridization reference has gene name and id as a string, split it and make a data.frame of the split
out <- strsplit(as.character(RNAseq_BRCA$`Hybridization REF`), "\\|")
test <- as.data.frame(do.call(rbind, out))

#get just the columns with raw counts
RNAseq_BRCA <- dplyr::select(RNAseq_BRCA, contains("raw_count"))

#Make a new column for gene name using the split from above
RNAseq_BRCA$Gene  <- test$V1

#remove the first row
RNAseq_BRCA <- RNAseq_BRCA[-1,]

#get rid of genes with ? or duplicate names
RNAseq_BRCA <- subset(RNAseq_BRCA, Gene != "?")
RNAseq_BRCA <- subset(RNAseq_BRCA, Gene != "SLC35E2")

#make a new data.frame with RNAseq_BRCA data with row.names = gene name
RNAseq_BRCA <- data.frame(RNAseq_BRCA, row.names = RNAseq_BRCA$Gene)

#get rid of the the gene name column
RNAseq_BRCA$Gene <- NULL

#make a new data.frame for RNAseq_BRCA where all of the raw counts are numeric
RNAseq_BRCA_df <- as.data.frame(row.names = row.names(RNAseq_BRCA), lapply(RNAseq_BRCA, function(x) as.numeric(as.character(x))))

#the next few lines are based off of this example: https://www.biostars.org/p/153013/
#quick check to see how many normal and primary samples there are. 
#In TCGA barcode Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29.
#the 14th character will be 0 for tumor and 1 for normal
table(substr(colnames(RNAseq_BRCA_df),14,14))

#index tumor (t) or normal (n) based on barcode value
n_index <- which(substr(colnames(RNAseq_BRCA_df),14,14) == '1')
t_index <- which(substr(colnames(RNAseq_BRCA_df),14,14) == '0')

#get CPM values for RNAseq data
logCPM <- cpm(RNAseq_BRCA_df, log=TRUE)

#define scal function, calculates a modified z-score, see link above for more details
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}

#apply scal function to the logCPM data 
z_rna <- scal(logCPM,logCPM[,n_index])

#transpose z_rna and make new RNAseq_BRCA data.frame
RNAseq_BRCA <- as.data.frame(t(z_rna))


#make a column for sample barcode based on row.names
RNAseq_BRCA$Sample <- row.names(RNAseq_BRCA)

#split the sample barcode
out <- strsplit(as.character(RNAseq_BRCA$Sample), "\\.")
test <- as.data.frame(do.call(rbind, out))

#make the patient barcode using the first three values of the split
test$bcr_patient_barcode  <- do.call(paste, c(test[c("V1", "V2", "V3")], sep = "-")) 

#add the patient_barcode to the data.frame
RNAseq_BRCA$patient_barcode <- test$bcr_patient_barcode

#make the sample barcode using the first three values of the split
test$Sample_Barcode <- do.call(paste, c(test[c("V1", "V2", "V3","V4")], sep = "-")) 

#add the sample_barcode to the data.frame
RNAseq_BRCA$Sample_Barcode <- test$Sample_Barcode

#add the sample type
RNAseq_BRCA$Sample_Type <- test$V4

#rename the sample types, add DROP for the B samples for each normal, primary or metastatic
RNAseq_BRCA$Sample_Type <- gsub("11A", "Normal", RNAseq_BRCA$Sample_Type)
RNAseq_BRCA$Sample_Type <- gsub("11B", "DROP", RNAseq_BRCA$Sample_Type)

RNAseq_BRCA$Sample_Type <- gsub("01A", "Primary", RNAseq_BRCA$Sample_Type)
RNAseq_BRCA$Sample_Type <- gsub("01B", "DROP", RNAseq_BRCA$Sample_Type)

RNAseq_BRCA$Sample_Type <- gsub("06A", "Metastatic", RNAseq_BRCA$Sample_Type)
RNAseq_BRCA$Sample_Type <- gsub("06B", "DROP", RNAseq_BRCA$Sample_Type)

#drop the b samples
RNAseq_BRCA <- subset(RNAseq_BRCA, !Sample_Type == "DROP")

#subset only the primary tumor samples
RNAseq_BRCA_primary <- subset(RNAseq_BRCA, Sample_Type == "Primary")

#load annotation file from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4911051/ 
Annotation <- fread("journal.pone.0157368.s008.txt")

#keep only columns 2-4 of annotation file
Annotation <- Annotation[,2:4]

#merge with RNAseq data, patient barcode works here because there are no normal or metastic samples
RNAseq_BRCA_primary <- merge(RNAseq_BRCA_primary, Annotation, by.y = "BARCODE", by.x = "patient_barcode", all.x = TRUE)

#subset only the normal samples
RNAseq_BRCA_normal <- subset(RNAseq_BRCA, Sample_Type == "Normal")

#defne TNBC and PAM50 as Normal or Normal Breast
RNAseq_BRCA_normal$TNBC <- "Normal"
RNAseq_BRCA_normal$PAM50 <- "Normal Breast"

#bind primary and normal data
RNAseq_BRCA_primary_normal <- rbind(RNAseq_BRCA_normal, RNAseq_BRCA_primary)

#make plots of interest and perform hypothesis testing
ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(TNBC)), aes(as.factor(TNBC), ADAR)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "NO", "YES"), labels=c("NO" = "Non TNBC", "YES" = "TNBC")) + ylab("ADAR Expression\n Z-Score") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16) 
ggsave("TNBC_v_ADAR.tiff", width = 5, height = 3, units = "in")
ggplot(subset(RNAseq_BRCA_primary_normal, TNBC %in% c("YES", "NO")) , aes(as.factor(TNBC), EIF2AK2)) + xlab("") + 
  scale_x_discrete(limits= c("NO", "YES"), labels=c("NO" = "Non TNBC", "YES" = "TNBC")) + ylab("PKR Expression\n Z-Score") + 
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16)
ggsave("TNBC_v_PKR.tiff")

ggplot(RNAseq_BRCA_primary_normal , aes(as.factor(TNBC), EIF2AK2)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "NO", "YES"), labels=c("Normal" = "Normal", "NO" = "Non TNBC", "YES" = "TNBC")) + ylab("PKR Expression\n Z-Score") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16)
ggsave("BRCA_v_PKR.tiff", width = 5, height = 3, units = "in")

ks.test(RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$TNBC == "YES"], RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$TNBC == "NO"])
t.test(RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$TNBC == "YES"], RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$TNBC == "NO"])


ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(PAM50)), aes(as.factor(PAM50), ADAR)) + xlab("") + ylab("ADAR Expression\n Z-Score") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16) +
  scale_x_discrete(limits= c("Normal Breast", "Normal", "Basal", "Her2", "LumA", "LumB"), labels=c("Normal\nBreast", "PAM50\nNormal", "Basal", "Her2", "LumA", "LumB"))
ggsave("PAM50_v_ADAR.tiff")

ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(PAM50)), aes(as.factor(PAM50), EIF2AK2)) + xlab("") + ylab("PKR Expression\n Z-Score") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16) + 
  scale_x_discrete(limits= c("Normal Breast", "Normal", "Basal", "Her2", "LumA", "LumB"), labels=c("Normal\nBreast", "PAM50\nNormal", "Basal", "Her2", "LumA", "LumB"))
ggsave("PAM50_v_PKR.tiff")

ks.test(RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "Basal"], RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "Her2"])
ks.test(RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "Basal"], RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "LumA"])
ks.test(RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "Basal"], RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "LumB"])
ks.test(RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "LumA"], RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "Her2"])
ks.test(RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "LumB"], RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "Her2"])
ks.test(RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "LumA"], RNAseq_BRCA_primary_normal$EIF2AK2[RNAseq_BRCA_primary_normal$PAM50 == "LumB"])

#read core ISGs list
ISG <- read.delim("Core_ISGs.txt")

#select RNAseq data for only Core ISGs
RNAseq_BRCA_primary_normal_ISG <- dplyr::select(RNAseq_BRCA_primary_normal, one_of(as.character(ISG$ISGs)))

#make a data.frame for Core ISG Score (ISG_median - median z-score of all ISGs per sample)
ISG_scores <- data.frame(ISG_median = rowMedians(as.matrix(RNAseq_BRCA_primary_normal_ISG)), Sample_Barcode = RNAseq_BRCA_primary_normal$Sample_Barcode)

#merge ISG_scores with RNAseq data
RNAseq_BRCA_primary_normal <- merge(RNAseq_BRCA_primary_normal, ISG_scores, by = "Sample_Barcode")

#make plots of interest and perform hypothesis tests
ggplot(subset(RNAseq_BRCA_primary_normal, TNBC %in% c("YES", "NO")) , aes(as.factor(TNBC), ISG_median)) + xlab("") + 
  scale_x_discrete(limits= c("NO", "YES"), labels=c("NO" = "Non TNBC", "YES" = "TNBC")) + ylab("Core ISG Score") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16)
ggsave("TNBC_v_ISG.tiff")

ggplot(RNAseq_BRCA_primary_normal , aes(as.factor(TNBC), ISG_median)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "NO", "YES"), labels=c("Normal" = "Normal", "NO" = "Non TNBC", "YES" = "TNBC")) + ylab("Core ISG Score") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16)
ggsave("BRCA_v_ISG.tiff", width = 5, height = 3, units = "in")

ks.test(RNAseq_BRCA_primary_normal$ISG_median[RNAseq_BRCA_primary_normal$TNBC == "YES"], RNAseq_BRCA_primary_normal$ISG_median[RNAseq_BRCA_primary_normal$TNBC == "NO"])
t.test(RNAseq_BRCA_primary_normal$ISG_median[RNAseq_BRCA_primary_normal$TNBC == "YES"], RNAseq_BRCA_primary_normal$ISG_median[RNAseq_BRCA_primary_normal$TNBC == "NO"])

ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(PAM50)), aes(as.factor(PAM50), ISG_median)) + xlab("") + 
   ylab("ISG Core Score") + geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science()
ggsave("PAM50_v_ISG.tiff")


#select RNAseq data for only Core ISGs
RNAseq_BRCA_ISG <- dplyr::select(RNAseq_BRCA, one_of(as.character(ISG$ISGs)))

#make a data.frame for Core ISG Score (ISG_median - median z-score of all ISGs per sample)
ISG_scores <- data.frame(ISG_median = rowMedians(as.matrix(RNAseq_BRCA_ISG)), Sample_Barcode = RNAseq_BRCA$Sample_Barcode)

#merge ISG_scores with RNAseq data
RNAseq_BRCA<- merge(RNAseq_BRCA, ISG_scores, by = "Sample_Barcode")

#make plots of interest and perform hypothesis tests
ggplot(subset(RNAseq_BRCA, !is.na(Sample_Type)), aes(as.factor(Sample_Type), ADAR)) + xlab("") + ylab("ADAR Expression Z-Score") + geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science()
ggsave("Sample_v_ADAR.tiff")

ggplot(subset(RNAseq_BRCA, !is.na(Sample_Type)), aes(as.factor(Sample_Type), EIF2AK2)) + xlab("") + ylab("EIF2AK2 Expression Z-Score") + geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science()
ggsave("Sample_v_PKR.tiff")

ggplot(subset(RNAseq_BRCA, !is.na(Sample_Type)), aes(as.factor(Sample_Type), ISG_median)) + xlab("") + ylab("ISG Core Score") + geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science()
ggsave("Sample_v_ISG.tiff")


#make a data.frame with just Core_ISG_Score, barcodes and TNBC status
TCGA_isg <- data.frame(Core_ISG_Score = RNAseq_BRCA_primary_normal$ISG_median, patient_barcode = RNAseq_BRCA_primary_normal$patient_barcode, TNBC = RNAseq_BRCA_primary_normal$TNBC, Sample_Barcode = RNAseq_BRCA_primary_normal$Sample_Barcode)


#read breast cancer clinical data
Clinical_BRCA <- readTCGA("/Users/cottr/Box Sync/TCGA/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.merged_only_clinical_clin_format.txt", dataType = 'clinical')

#use the survivalTCGA function to extract survival data
BRCA_survival <- survivalTCGA(Clinical_BRCA)

#merge survival data with primary tumor sample RNAseq
BRCA <- merge(BRCA_survival, RNAseq_BRCA_primary, by.x = "bcr_patient_barcode", by.y="patient_barcode")

#remove data with times <0
BRCA <- subset(BRCA, !times < 0)

#determine ADAR expression cutoff using surv_cutpoint function

BRCA.surv_rnaseq.cut <- surv_cutpoint(
  BRCA,
  time = "times",
  event = "patient.vital_status",
  variables = "ADAR"
)
summary(BRCA.surv_rnaseq.cut)
plot(BRCA.surv_rnaseq.cut, "ADAR")

#score patients by ADAR expression above or below cutoff
BRCA$ADAR_score <- ifelse(BRCA$ADAR < BRCA.surv_rnaseq.cut$cutpoint$cutpoint, "low", "high")

#fit the survival curve and plot
fit_ADAR <- survfit(Surv(times, patient.vital_status)
                    ~ ADAR_score , data = BRCA)


ggsurvplot(fit_ADAR, data = BRCA, size = 1,                 # change line size
           palette =
             c("#E7B800", "#2E9FDF"),# custom color palettes
           conf.int = TRUE,          # Add confidence interval
           pval = TRUE,              # Add p-value
           legend.title = "ADAR\nExpression",
           legend = "right",
           legend.labs = c("High", "Low"),
           xlab = "Days",
           risk.table = FALSE,        # Add risk table
           ggtheme = theme_science(base_size = 16)      # Change ggplot2 theme
)
ggsave("ADAR_cutoff.tiff")




#read isoform RNAseq data from FireBrowse
RNAseq_BRCA_isoforms <- fread("BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt")

#select only isoforms of interest (ARF, INFK4A, p150, p110)
RNAseq_BRCA_isoforms <- subset(RNAseq_BRCA_isoforms, grepl("uc001ffh|uc001ffj|uc001ffi|uc001ffk|uc003zpl|uc003zpk", RNAseq_BRCA_isoforms$`Hybridization REF`))

#rename to protein name
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc001ffh.2", "p150-1", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc001ffi.2", "p150-2", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc001ffj.2", "p150-3", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc001ffk.2", "p110", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc003zpl.2", "ARF", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc003zpk.2", "INK4A", RNAseq_BRCA_isoforms$`Hybridization REF`)

#make new data.frame with isoforms names as row.names
RNAseq_BRCA_isoforms <- data.frame(RNAseq_BRCA_isoforms, row.names = RNAseq_BRCA_isoforms$`Hybridization REF`)

#remove hybridization reference column
RNAseq_BRCA_isoforms$Hybridization.REF <- NULL

#transpose data.frame
RNAseq_BRCA_isoforms <- as.data.frame(t(RNAseq_BRCA_isoforms))

#make a new column for sample based on row.names
RNAseq_BRCA_isoforms$Sample <- row.names(RNAseq_BRCA_isoforms)

#split the sample identifier
out <- strsplit(as.character(RNAseq_BRCA_isoforms$Sample), "\\.")
test <- as.data.frame(do.call(rbind, out))

#build a patient_barcode
test$bcr_patient_barcode  <- do.call(paste, c(test[c("V1", "V2", "V3")], sep = "-")) 

#make a new column with the patient_barcode
RNAseq_BRCA_isoforms$patient_barcode <- test$bcr_patient_barcode

#build a sample barcode
test$Sample_Barcode <- do.call(paste, c(test[c("V1", "V2", "V3","V4")], sep = "-")) 

#make a new column wit the sample barcode
RNAseq_BRCA_isoforms$Sample_Barcode <- test$Sample_Barcode

#make a new column with sample type
RNAseq_BRCA_isoforms$Sample_Type <- test$V4

#rename the sample types as above
RNAseq_BRCA_isoforms$Sample_Type <- gsub("11A", "Normal", RNAseq_BRCA_isoforms$Sample_Type)
RNAseq_BRCA_isoforms$Sample_Type <- gsub("11B", "DROP", RNAseq_BRCA_isoforms$Sample_Type)

RNAseq_BRCA_isoforms$Sample_Type <- gsub("01A", "Primary", RNAseq_BRCA_isoforms$Sample_Type)
RNAseq_BRCA_isoforms$Sample_Type <- gsub("01B", "DROP", RNAseq_BRCA_isoforms$Sample_Type)

RNAseq_BRCA_isoforms$Sample_Type <- gsub("06A", "Metastatic", RNAseq_BRCA_isoforms$Sample_Type)
RNAseq_BRCA_isoforms$Sample_Type <- gsub("06B", "DROP", RNAseq_BRCA_isoforms$Sample_Type)

#remove B samples
RNAseq_BRCA_isoforms <- subset(RNAseq_BRCA_isoforms, !Sample_Type == "DROP")

#calculate total p150 expresssion, sum(p150-1,p150-2,p150-3)
RNAseq_BRCA_isoforms$p150 <- as.numeric(as.character(RNAseq_BRCA_isoforms$`p150-1`)) + 
  as.numeric(as.character(RNAseq_BRCA_isoforms$`p150-2`)) +
  as.numeric(as.character(RNAseq_BRCA_isoforms$`p150-3`))

#make p110 numeric
RNAseq_BRCA_isoforms$p110 <- as.numeric(as.character(RNAseq_BRCA_isoforms$p110))

#subset only primary tumor samples
RNAseq_BRCA_isoforms_primary <- subset(RNAseq_BRCA_isoforms , Sample_Type == "Primary")

#merge with annotation data
RNAseq_BRCA_isoforms_primary <- merge(RNAseq_BRCA_isoforms_primary, Annotation, by.y = "BARCODE", by.x = "patient_barcode", all.x = TRUE)

#subset only normal samples
RNAseq_BRCA_isoforms_normal <- subset(RNAseq_BRCA_isoforms, Sample_Type == "Normal")

#fill in TNBC and PAM50 status for normal samples
RNAseq_BRCA_isoforms_normal$TNBC <- "Normal"
RNAseq_BRCA_isoforms_normal$PAM50 <- "Normal Breast"

#bind primary and normal
RNAseq_BRCA_isoforms <- rbind(RNAseq_BRCA_isoforms_normal, RNAseq_BRCA_isoforms_primary)

#make a summary data.frame
TCGA_summary <- merge(RNAseq_BRCA_isoforms, TCGA_isg, by.y = "Sample_Barcode", by.x = "Sample_Barcode", all = T)

#remove the extra TNBC column
TCGA_summary$TNBC.x <- NULL
names(TCGA_summary)[names(TCGA_summary) == 'TNBC.y'] <- 'TNBC'

#remove the extra patient_barcode column
TCGA_summary$patient_barcode.x<- NULL
names(TCGA_summary)[names(TCGA_summary) == 'patient_barcode.y'] <- 'patient_barcode'

#write the summary file
write.csv(TCGA_summary, "TCGA_summary.csv", row.names = FALSE)


#make plots of interest and perform hypothesis testing
ggplot(subset(RNAseq_BRCA_isoforms, !is.na(TNBC)), aes(as.factor(TNBC), p150)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "NO", "YES"), labels=c("NO" = "Non TNBC", "YES" = "TNBC")) + ylab("p150 RSEM") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16)
ggsave("TNBC_v_p150.tiff")

ggplot(subset(RNAseq_BRCA_isoforms, !is.na(TNBC)), aes(as.factor(TNBC), p110)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "NO", "YES"), labels=c("NO" = "Non TNBC", "YES" = "TNBC")) + ylab("p110 RSEM") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16)
ggsave("TNBC_v_p110.tiff")

ks.test(RNAseq_BRCA_isoforms_primary$EIF2AK2[RNAseq_BRCA_isoforms_primary$TNBC == "YES"], RNAseq_BRCA_isoforms_primary$EIF2AK2[RNAseq_BRCA_isoforms_primary$TNBC == "NO"])
t.test(RNAseq_BRCA_isoforms_primary$EIF2AK2[RNAseq_BRCA_isoforms_primary$TNBC == "YES"], RNAseq_BRCA_isoforms_primary$EIF2AK2[RNAseq_BRCA_isoforms_primary$TNBC == "NO"])


ggplot(subset(RNAseq_BRCA_isoforms_primary, !is.na(PAM50)), aes(as.factor(PAM50), p150)) + xlab("") + ylab("p150 RSEM") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16)
ggsave("PAM50_v_p150.tiff")

ggplot(subset(RNAseq_BRCA_isoforms_primary, !is.na(PAM50)), aes(as.factor(PAM50), p110)) + xlab("") + ylab("p110 RSEM") +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16)
ggsave("PAM50_v_p110.tiff")

ggplot(subset(RNAseq_BRCA_isoforms_primary, !is.na(TNBC)), aes(p150, p110, colour = as.factor(TNBC))) + geom_point()


ggplot(subset(RNAseq_BRCA_isoforms, !is.na(Sample_Type)), aes(as.factor(Sample_Type), p150)) + xlab("") + ylab("p150") + geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science()
ggsave("Sample_v_p150.tiff")

ggplot(subset(RNAseq_BRCA_isoforms, !is.na(Sample_Type)), aes(as.factor(Sample_Type), p110)) + xlab("") + ylab("p110") + geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science()
ggsave("Sample_v_p110.tiff")

