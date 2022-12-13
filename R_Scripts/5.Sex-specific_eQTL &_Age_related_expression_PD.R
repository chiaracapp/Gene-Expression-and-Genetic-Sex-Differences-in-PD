### ### ### ### ### ### ### ### 
###  Sex-specific eQTLs PD  ###
### ### ### ### ### ### ### ### 


### Packages required:
library("DESeq2")
library("sva")
library("dplyr")
library("rstatix")
library("corrplot")
library("circlize")
library("limma")
library("ggplot2")
library("ggpubr")
library("lsmeans")


load("/dds_84_qSvs.Rda")
coldata <- colData(dds)
coldata <- as.data.frame(coldata)

### Genotype file as YG IDs, add them to coldata
YG_IDs <- read.csv2("/Sample_ID_New.csv", sep = ",")
YG_IDs <- YG_IDs[YG_IDs$X %in% rownames(coldata),]
coldata <- merge(coldata, YG_IDs, by.x = "row.names", by.y = "X")
rownames(coldata) <- coldata$Row.names
coldata$Row.names <- NULL

###  ENSG and genes ID
Genes_ID <- read.csv(file="/Genes_ID.csv")

### Load file with 90 SNPs from Nalls et.al 2019
Risk_genes <- read.csv2("Nalls_2019_RiskSNP.csv", sep = ",")

### Delete NA
Risk_genes <- Risk_genes[!is.na(Risk_genes$Ensembl.ID),]

### Keep only genes present in dds: genes actually detected with RNA-Seq
Risk_genes_T <- Risk_genes[Risk_genes$Ensembl.ID %in% Genes_ID$gene,]

### Remove duplicates to count the number of expressed genes
Risk_genes_T_notDup <- Risk_genes_T[!duplicated(Risk_genes_T$Ensembl.ID),] 

### Load genotypes
Genotypes <- read.csv2("/Genotypes_PDsnps.txt", sep = "")

### Keep only 83 samples
Genotypes <- Genotypes[Genotypes$FID %in% coldata$Sample_Name,] 
rownames(Genotypes) <- Genotypes$FID

### Keep only SNPs of interest in genotype
genotypes_col <- colnames(Genotypes)
genotypes_col_1 <- stringr::str_extract(genotypes_col, "[^.]*.[^.]*")
colnames(Genotypes) <- genotypes_col_1

X <- rep("X", length(Risk_genes_T$SNP))
Risk_genes_T$uniqID <- paste0(X,Risk_genes_T$CHR, ".", Risk_genes_T$BP )

### Make formula for each combination ENSG00000117020 ~ X1.243717682
Risk_genes_T$formula <- paste0(Risk_genes_T$Ensembl.ID, " ~ ", Risk_genes_T$uniqID)

### Remove row with rs76763715 X1.155305634 all samples are 0
Risk_genes_T <- Risk_genes_T[!(Risk_genes_T$SNP == "rs76763715"),]

### Remove row with rs34637584, no genotype
Risk_genes_T <- Risk_genes_T[!(Risk_genes_T$SNP == "rs34637584"),]

Genotypes <- Genotypes[,colnames(Genotypes) %in% Risk_genes_T$uniqID] 
intersect(Risk_genes_T$uniqID ,colnames(Genotypes))
setdiff(Risk_genes_T$uniqID,colnames(Genotypes))
setdiff(colnames(Genotypes), Risk_genes_T$uniqID)

# Transform counts
vsd <- vst(dds, blind = FALSE) 
vsd_age <- vsd # for later
vsd_age_1 <- vsd # for later

# Adjust counts for covariates
covariates <- coldata[,c(14:18)] # qSVs
assay(vsd) <- limma::removeBatchEffect(assay(vsd),batch = vsd$Braak_aSyn_stage, batch2 = scale(vsd$Age_death), covariates =  as.matrix(covariates),  model.matrix(~ vsd$Sex))
assay(vsd)

Adjusted_counts <- assay(vsd)

### Keep only genes of interest
Adjusted_counts <- Adjusted_counts[rownames(Adjusted_counts) %in% Risk_genes_T$Ensembl.ID,]

### Divide samples in 2 groups
Adjusted_counts <- t(Adjusted_counts)
Adjusted_counts <- as.data.frame(Adjusted_counts)
Adjusted_counts$Sex <- coldata$Sex
Adjusted_counts$YGID <- coldata$Sample_Name
Adjusted_counts <- Adjusted_counts[!(Adjusted_counts$YGID == "YG498"),] # remove sample without genotype YG498

Adjusted_counts <- merge(Adjusted_counts, Genotypes, by.x = "YGID", by.y = "row.names")

Adjusted_counts_F <- Adjusted_counts[Adjusted_counts$Sex == "F",]
Adjusted_counts_M <- Adjusted_counts[Adjusted_counts$Sex == "M",]


#### Females 

### Remove Snps with only one genotype
Risk_genes_T_F <- Risk_genes_T[!(Risk_genes_T$SNP == "rs114138760"),]
Risk_genes_T_F <- Risk_genes_T_F[!(Risk_genes_T_F$SNP == "rs75859381"),]
Risk_genes_T_F <- Risk_genes_T_F[!(is.na(Risk_genes_T_F$SNP)),]

results_F <- lapply(Risk_genes_T_F$formula, function(x) lm(x, data = Adjusted_counts_F))
results_F <- setNames(results_F,Risk_genes_T_F$formula )

### Extract p-value and add to table
p_val_F <- numeric()

for (x in results_F) {
  out_pval_F <- ((pf(summary(x)$fstatistic["value"], summary(x)$fstatistic["numdf"], summary(x)$fstatistic["dendf"],lower.tail = FALSE)))
  p_val_F <- c(p_val_F, out_pval_F)
}


Risk_genes_T_F$p_value_F <- p_val_F

Risk_genes_T_F$p_adj_F <- p.adjust(Risk_genes_T_F$p_value_F, "BH", n= length(Risk_genes_T_F$p_value_F))
Risk_genes_T_F$sig_F <- ifelse(Risk_genes_T_F$p_value_F < 0.05, "Significant", "Not Significant")
sum(Risk_genes_T_F$sig_F == "Significant") 
sum(Risk_genes_T_F$p_adj_F < 0.05) 


#### Males

### Remove Snps with only one genotype
Risk_genes_T_M <- Risk_genes_T[!(Risk_genes_T$SNP %in% c("rs117896735")),]
results_M <- lapply(Risk_genes_T_M$formula, function(x) lm(x, data = Adjusted_counts_M))
results_M <- setNames(results_M, Risk_genes_T_M$formula )

### Extract p value and add to table
p_val_M <- numeric()
for (x in results_M) {
  out_pval_M <- ((pf(summary(x)$fstatistic["value"], summary(x)$fstatistic["numdf"], summary(x)$fstatistic["dendf"],lower.tail = FALSE)))
  p_val_M <- c(p_val_M, out_pval_M)
}

Risk_genes_T_M$p_value_M <- p_val_M
Risk_genes_T_M$p_adj_M <- p.adjust(Risk_genes_T_M$p_value_M, "BH", n= length(Risk_genes_T_M$p_value_M))
Risk_genes_T_M$sig_M <- ifelse(Risk_genes_T_M$p_value_M < 0.05, "Significant", "Not Significant")
sum(Risk_genes_T_M$sig_M == "Significant") #7
sum(Risk_genes_T_M$p_adj_M < 0.05) 


write.csv(as.data.frame(Risk_genes_T_M), file="/M_PD_eQTLs.csv",  row.names = FALSE)
write.csv(as.data.frame(Risk_genes_T_F), file="/F_PD_eQTLs.csv",  row.names = FALSE)


###### Compare F and M slopes for the significant eQTLS

### Make formula for each combination ENSG00000117020 ~ X1.243717682*Sex
Risk_genes_T$formula_slopes <- paste0(Risk_genes_T$Ensembl.ID, " ~ ", Risk_genes_T$uniqID,"*","Sex")
Risk_genes_T$formula_slopes

results_slopes <- lapply(Risk_genes_T$formula_slopes, function(x) lm(x, data = Adjusted_counts))

### Significant in males
slope_VAMP4  <- lstrends(results_slopes[[4]], "Sex",var = "X1.171719769")
slope_SATB1  <- lstrends(results_slopes[[14]], "Sex",var = "X3.18361759")
slope_CAB39L  <- lstrends(results_slopes[[53]], "Sex",var = "X13.49927732")
slope_VPS13C  <- lstrends(results_slopes[[59]], "Sex",var = "X15.61997385")
slope_FAM171A2  <- lstrends(results_slopes[[67]], "Sex",var = "X17.42434630")
slope_CRLS1  <- lstrends(results_slopes[[77]], "Sex",var = "X20.6006041")
slope_MED13  <- lstrends(results_slopes[[90]], "Sex",var = "X17.59917366")

### Significant in females
slope_KCNIP3  <- lstrends(results_slopes[[10]], "Sex",var = "X2.96000943")
slope_KPNA1  <- lstrends(results_slopes[[16]], "Sex",var = "X3.122196892")
slope_GAK  <- lstrends(results_slopes[[20]], "Sex",var = "X4.925376")
# CAB39L already checked
slope_GCH1  <- lstrends(results_slopes[[56]], "Sex",var = "X14.55348869")
slope_GALC  <- lstrends(results_slopes[[58]], "Sex",var = "X14.88464264")


Slopes <- list(slope_VAMP4,slope_SATB1,slope_CAB39L,slope_VPS13C,slope_FAM171A2,slope_CRLS1, slope_MED13, slope_KCNIP3,slope_KPNA1, slope_GAK, slope_GCH1,slope_GALC)

### Compare slopes
results_slopes_1 <- lapply(Slopes, pairs)
results_pvalues <- lapply(results_slopes_1, summary)

### Extract p-values
p_val_slopes <- lapply(results_pvalues, function(x)(x$p.value))
names(p_val_slopes) <- c("VAMP4","SATB1","CAB39L","VPS13C","FAM171A2","CRLS1", "MED13", "KCNIP3","KPNA1", "GAK", "GCH1","GALC")

trends <-  as.data.frame(p_val_slopes)
trends <- t(trends)
trends <-  as.data.frame(trends)
trends$Genes <- rownames(trends)

write.csv(as.data.frame(trends), file="/All_PD_eQTLs_trends.csv",  row.names = FALSE)



### ### ### ### ### ### ### ### ### ### ### ### 
###   Sex-specific age-related expression   ###
### ### ### ### ### ### ### ### ### ### ### ### 


Adjusted_counts_1 <- assay(vsd_age_1)

### Keep only risk genes
Adjusted_counts_1 <- Adjusted_counts_1[rownames(Adjusted_counts_1) %in% Risk_genes_T$Ensembl.ID,]  


### Divide samples in 2 groups
Adjusted_counts_1 <- t(Adjusted_counts_1)
Adjusted_counts_1 <- as.data.frame(Adjusted_counts_1)
Adjusted_counts_1$Sex <- coldata$Sex
Adjusted_counts_1$Age <- coldata$Age_death
Adjusted_counts_1 <- Adjusted_counts_1[!(Adjusted_counts_1$Age == "102"),]


### Make formula for each combination ENSG00000117020 ~ Age*Sex
Risk_genes_T$formula_Age <- paste0(Risk_genes_T$Ensembl.ID, " ~ ", "Age","*","Sex")
Risk_genes_T$formula_Age


### Linear model
results_Age <- lapply(Risk_genes_T$formula_Age, function(x) lm(x, data = Adjusted_counts_1))

### Obtain slopes
Slopes_Age <- lapply(results_Age, function(x) lstrends(x, "Sex", var = "Age"))

### Compare slopes
results_Age_F <- lapply(Slopes_Age, pairs)
results_Age_pvalues <- lapply(results_Age_F, summary)


### Extract p-values
p_val_F <- lapply(results_Age_pvalues, function(x)(x$p.value))

Risk_genes_T$pvalue_AGE <- p_val_F
Risk_genes_T$p_adj_AGE <- p.adjust(Risk_genes_T$pvalue_AGE, "BH", n= length(Risk_genes_T$pvalue_AGE)) 
Risk_genes_T$sig_AGE <- ifelse(Risk_genes_T$pvalue_AGE < 0.05, "Significant", "Not Significant")

Risk_genes_T <- as.data.frame(Risk_genes_T)
Risk_genes_T <- apply(Risk_genes_T,2,as.character)

write.csv(Risk_genes_T, file="/PD_AGE.csv")


##### Does the genotype influence the age-related expression?

### Adjust counts for covariates do not correct for age
assay(vsd_age) <- limma::removeBatchEffect(assay(vsd_age),batch = vsd_age$Braak_aSyn_stage, batch2 = NULL, covariates =  as.matrix(covariates),  model.matrix(~ vsd_age$Sex))

Adjusted_counts_age <- assay(vsd_age)

### Keep only interesting genes
Adjusted_counts_age <- Adjusted_counts_age[rownames(Adjusted_counts_age) %in% Risk_genes_T$Ensembl.ID,]  ### 86 genes

### Keep only genes significant in Age analysis
Sig_age <- read.csv2("/PD_AGE.csv", sep = ",")
Sig_age <- Sig_age[!(Sig_age$sig_AGE=="Not Significant"),]

Adjusted_counts_age <- Adjusted_counts_age[rownames(Adjusted_counts_age) %in% Sig_age$Ensembl.ID,] 

### Divide samples in 2 groups
Adjusted_counts_age <- t(Adjusted_counts_age)
Adjusted_counts_age <- as.data.frame(Adjusted_counts_age)
Adjusted_counts_age$Sex <- coldata$Sex
Adjusted_counts_age$YGID <- coldata$Sample_Name
Adjusted_counts_age$Age <- coldata$Age_death
Adjusted_counts_age <- Adjusted_counts_age[!(Adjusted_counts_age$YGID == "YG498"),] # remove sample without genotype YG498


Adjusted_counts_age <- merge(Adjusted_counts_age, Genotypes, by.x = "YGID", by.y = "row.names")
str(Adjusted_counts_age)

Adjusted_counts_age_F <- Adjusted_counts_age[Adjusted_counts_age$Sex == "F",]
Adjusted_counts_age_M <- Adjusted_counts_age[Adjusted_counts_age$Sex == "M",]

Adjusted_counts_age_F[,13:93] <- lapply(Adjusted_counts_age_F[,13:93], factor)
Adjusted_counts_age_M[,13:93] <- lapply(Adjusted_counts_age_M[,13:93], factor)

### Make formula for each combination 
Risk_genes_T_Age <- Risk_genes_T[Risk_genes_T$Ensembl.ID %in% Sig_age$Ensembl.ID,] 

Risk_genes_T_Age$formula_age <- paste0(Risk_genes_T_Age$Ensembl.ID, " ~ ","Age","*", Risk_genes_T_Age$uniqID)
Risk_genes_T_Age$formula_age

### Linear model
results_age_F <- lapply(Risk_genes_T_Age$formula_age, function(x) lm(x, data = Adjusted_counts_age_F))
results_age_M <- lapply(Risk_genes_T_Age$formula_age, function(x) lm(x, data = Adjusted_counts_age_M))

### In Females

### Obtain slopes
slope_1 <- lstrends(results_age_F[[ 1 ]], "X1.205723572",var = "Age")
slope_2 <- lstrends(results_age_F[[ 2 ]], "X1.232664611",var = "Age")
slope_3 <- lstrends(results_age_F[[ 3 ]], "X2.169110394",var = "Age")
slope_4 <- lstrends(results_age_F[[ 4 ]], "X4.77110365",var = "Age")
slope_5 <- lstrends(results_age_F[[ 5 ]], "X4.170583157",var = "Age")
slope_6 <- lstrends(results_age_F[[ 6 ]], "X11.10558777",var = "Age")
slope_7 <- lstrends(results_age_F[[ 7 ]], "X14.88464264",var = "Age")
slope_8 <- lstrends(results_age_F[[ 8 ]], "X17.76425480",var = "Age")
slope_9 <- lstrends(results_age_F[[ 9 ]], "X4.17968811",var = "Age")


Slopes_F_Age <- list(slope_1,slope_2,slope_3,slope_4,slope_5,slope_6,slope_7,slope_8,slope_9)


### Compare slopes
results_age_F_F <- lapply(Slopes_F_Age, pairs)
results_pvalues_Age <- lapply(results_age_F_F, summary)

### extract p vlaues
p_val_F_Age <- lapply(results_pvalues_Age, function(x)(x$p.value))
names(p_val_F_Age) <- Risk_genes_T_Age$formula_age


p_values_F_Age <- as.data.frame(do.call(rbind, p_val_F_Age))

p_values_F_Age$sig <- ifelse((p_values_F_Age$V1 < 0.05|p_values_F_Age$V2 < 0.05|p_values_F_Age$V3 < 0.05), "Significant", "Not Significant")

### In males

### Obtain slopes
slope_1_M <- lstrends(results_age_M[[ 1 ]], "X1.205723572",var = "Age")
slope_2_M <- lstrends(results_age_M[[ 2 ]], "X1.232664611",var = "Age")
slope_3_M <- lstrends(results_age_M[[ 3 ]], "X2.169110394",var = "Age")
slope_4_M <- lstrends(results_age_M[[ 4 ]], "X4.77110365",var = "Age")
slope_5_M <- lstrends(results_age_M[[ 5 ]], "X4.170583157",var = "Age")
slope_6_M <- lstrends(results_age_M[[ 6 ]], "X11.10558777",var = "Age")
slope_7_M <- lstrends(results_age_M[[ 7 ]], "X14.88464264",var = "Age")
slope_8_M <- lstrends(results_age_M[[ 8 ]], "X17.76425480",var = "Age")
slope_9_M <- lstrends(results_age_M[[ 9 ]], "X4.17968811",var = "Age")


Slopes_M_Age <- list(slope_1_M,slope_2_M,slope_3_M,slope_4_M,slope_5_M,slope_6_M,slope_7_M,slope_8_M,slope_9_M)

### Compare slopes
results_age_M_F <- lapply(Slopes_M_Age, pairs)
results_pvalues_M_Age <- lapply(results_age_M_F, summary)

### extract p vlaues
p_val_M_Age <- lapply(results_pvalues_M_Age, function(x)(x$p.value))
names(p_val_M_Age) <- Risk_genes_T_Age$formula_age


p_values_M_Age <- as.data.frame(do.call(rbind, p_val_M_Age))
p_values_M_Age$sig <- ifelse((p_values_M_Age$V1 < 0.05|p_values_M_Age$V2 < 0.05|p_values_M_Age$V3 < 0.05), "Significant", "Not Significant")

Age_results <- cbind(p_values_F_Age,p_values_M_Age)
colnames(Age_results_Age) <- c("0vs1_F", "0vs2_F","1vs2_F", "Sig_F", "0vs1_M", "0vs2_M","1vs2_M", "Sig_M")

write.csv(as.data.frame(Age_results), file="/PD_Geno_Age.csv",  row.names = TRUE)

