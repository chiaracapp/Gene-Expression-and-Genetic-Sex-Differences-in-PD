### ### ### ### ### ### ### ### ### ### ### 
###    Selected eQTLs Pleiotropy loci   ###
### ### ### ### ### ### ### ### ### ### ### 


### Packages required:
library("DESeq2")
library("sva")
library("dplyr")
library("rstatix")
library("corrplot")
library("circlize")
library("limma")
library("ggplot2")
library("lsmeans")
library("ggpubr")

load("/dds_84_qSvs.Rda")
coldata <- colData(dds)
coldata <- as.data.frame(coldata)

### Genotype file as YG IDs, add them to coldata
YG_IDs <- read.csv2("/Sample_ID_New.csv", sep = ",")
YG_IDs <- YG_IDs[YG_IDs$X %in% rownames(coldata),]

coldata <- merge(coldata, YG_IDs, by.x = "row.names", by.y = "X")
rownames(coldata) <- coldata$Row.names
coldata$Row.names <- NULL

### Load risk genes
sex_genes <- read.csv2("/Expressed_genes_Conj.csv", sep = ",")

### Remove HLA genes
sex_genes <- sex_genes[!grepl("HLA", sex_genes$symbol),]
sum(sex_genes$Trait == "Menopause") #49
sum(sex_genes$Trait == "Menarche") #68

### Load genotypes
Genotypes <- read.csv2("/Conjsnps.txt", sep = "")

### Keep only 83 samples
Genotypes <- Genotypes[Genotypes$FID %in% coldata$Sample_Name,]
rownames(Genotypes) <- Genotypes$FID

### Make formula for each combination ENSG00000117020 ~ X1.243717682
sex_genes$formula <- paste0(sex_genes$ensg, " ~ ", sex_genes$Geno)

### Transform counts
vsd <- vst(dds, blind = FALSE) 
vsd_Age_1 <- vsd
vsd_age <- vsd

# Adjust counts for covariates
covariates <- coldata[,c(14:18)] # qSVs
assay(vsd_1) <- limma::removeBatchEffect(assay(vsd_1),batch = vsd_1$Braak_aSyn_stage, batch2 = scale(vsd_1$Age_death), covariates =  as.matrix(covariates),  model.matrix(~ vsd_1$Sex))
assay(vsd_1)

Adjusted_counts <- assay(vsd_1)

### Keep only sex genes
Adjusted_counts <- Adjusted_counts[rownames(Adjusted_counts) %in% sex_genes$ensg,]  

### Divide samples in 2 groups
Adjusted_counts <- t(Adjusted_counts)
Adjusted_counts <- as.data.frame(Adjusted_counts)
Adjusted_counts$Sex <- coldata$Sex
Adjusted_counts$YGID <- coldata$Sample_Name
Adjusted_counts <- Adjusted_counts[!(Adjusted_counts$YGID == "YG498"),] # remove sample without genotype YG498


Adjusted_counts <- merge(Adjusted_counts, Genotypes, by.x = "YGID", by.y = "row.names")
str(Adjusted_counts)

Adjusted_counts_F <- Adjusted_counts[Adjusted_counts$Sex == "F",]
Adjusted_counts_M <- Adjusted_counts[Adjusted_counts$Sex == "M",]


#### Females 
results_F <- lapply(sex_genes$formula, function(x) lm(x, data = Adjusted_counts_F))
results_F <- setNames(results_F,sex_genes$formula )

### Extract p-value and add to table
p_val_F <- numeric()

for (x in results_F) {
  out_pval <- ((pf(summary(x)$fstatistic["value"], summary(x)$fstatistic["numdf"], summary(x)$fstatistic["dendf"],lower.tail = FALSE)))
  p_val_F <- c(p_val_F, out_pval)
}

sex_genes$p_value_F <- p_val_F

sex_genes$p_adj_F <- p.adjust(sex_genes$p_value_F, "BH", n= length(sex_genes$p_value_F))
sex_genes$sig_F <- ifelse(sex_genes$p_value_F < 0.05, "Significant", "Not Significant")
sum(sex_genes$sig_F == "Significant") #1
sum(sex_genes$p_adj_F < 0.05) #0


#### Males
results_M <- lapply(sex_genes$formula, function(x) lm(x, data = Adjusted_counts_M))
results_M <- setNames(results_M, sex_genes$formula )

# Extract p-value and add to table
p_val_M <- numeric()

for (x in results_M) {
  out_pval_M <- ((pf(summary(x)$fstatistic["value"], summary(x)$fstatistic["numdf"], summary(x)$fstatistic["dendf"],lower.tail = FALSE)))
  p_val_M <- c(p_val_M, out_pval_M)
}

sex_genes$p_value_M <- p_val_M
sex_genes$p_adj_M <- p.adjust(sex_genes$p_value_M, "BH", n= length(sex_genes$p_value_M))
sex_genes$sig_M <- ifelse(sex_genes$p_value_M < 0.05, "Significant", "Not Significant")
sum(sex_genes$sig_M == "Significant") 
sum(sex_genes$p_adj_M < 0.05) 

write.csv(as.data.frame(sex_genes), file="/Conj_eQTLs.csv",  row.names = FALSE)


###### Compare F and M slopes for the significant eQTLS

### Make formula for each combination ENSG00000117020 ~ X1.243717682*Sex
sex_genes$formula_slopes <- paste0(sex_genes$ensg, " ~ ", sex_genes$Geno,"*","Sex")
sex_genes$formula_slopes

results_slopes <- lapply(sex_genes$formula_slopes, function(x) lm(x, data = Adjusted_counts))

slope_ELK4  <- lstrends(results_slopes[[13]], "Sex",var = "X1.205739266.C.T_C")
slope_MCMBP  <- lstrends(results_slopes[[31]], "Sex",var = "X10.121708929.T.C_C")
slope_ERGIC2  <- lstrends(results_slopes[[42]], "Sex",var = "X12.30895251.T.C_C")
slope_SEC23IP  <- lstrends(results_slopes[[46]], "Sex",var = "X12.30895251.T.C_C")
slope_ZSWIM7  <- lstrends(results_slopes[[55]], "Sex",var = "X17.16081958.C.T_T")
slope_NIT2  <- lstrends(results_slopes[[92]], "Sex",var = "X3.100472848.T.C_C")

Slopes <- list(slope_ELK4,slope_MCMBP,slope_ERGIC2,slope_SEC23IP,slope_ZSWIM7,slope_NIT2)

### Compare slopes
results_slopes_1 <- lapply(Slopes, pairs)
results_pvalues <- lapply(results_slopes_1, summary)

### Extract p-values
p_val_slopes <- lapply(results_pvalues, function(x)(x$p.value))
names(p_val_slopes) <- c("ELK4","MCMBP","ERGIC2","SEC23IP","ZSWIM7","NIT2")

trends <-  as.data.frame(p_val_slopes)
trends <- t(trends)

write.csv(as.data.frame(trends), file="/Conj_eQTLs_trends.csv",  row.names = TRUE)


### ### ### ### ### ### ### ### ### ### ### ### 
###   Sex-specific age-related expression   ###
### ### ### ### ### ### ### ### ### ### ### ### 

Adjusted_counts_Age_1 <- assay(vsd_Age_1)

### Keep only 50 genes
Adjusted_counts_Age_1 <- Adjusted_counts_Age_1[rownames(Adjusted_counts_Age_1) %in% sex_genes$ensg,]  

### Divide samples in 2 groups
Adjusted_counts_Age_1 <- t(Adjusted_counts_Age_1)
Adjusted_counts_Age_1 <- as.data.frame(Adjusted_counts_Age_1)
Adjusted_counts_Age_1$Sex <- coldata$Sex
Adjusted_counts_Age_1$Age <- coldata$Age_death
Adjusted_counts_Age_1 <- Adjusted_counts_Age_1[!(Adjusted_counts_Age_1$Age == "102"),]

### Make formula for each combination ENSG00000117020 ~ Age*Sex
sex_genes$formula_Age_1 <- paste0(sex_genes$ensg, " ~ ", "Age","*","Sex")

### linear model
results_Age_1 <- lapply(sex_genes$formula_Age_1, function(x) lm(x, data = Adjusted_counts_Age_1))

### Obtain slopes
Slopes_Age_1 <- lapply(results_Age_1, function(x) lstrends(x, "Sex", var = "Age"))

### Compare slopes
results_F_Age_1 <- lapply(Slopes_Age_1, pairs)
results_pvalues_Age_1 <- lapply(results_F_Age_1, summary)


### Extract p-values
p_val_F_Age_1 <- lapply(results_pvalues_Age_1, function(x)(x$p.value))

sex_genes$pvalue_AGE <- p_val_F_Age_1
sex_genes$p_adj_AGE <- p.adjust(sex_genes$pvalue_AGE, "BH", n= length(sex_genes$pvalue_AGE)) ## nothing significant after multiple testing correction
sex_genes$sig_AGE <- ifelse(sex_genes$pvalue_AGE < 0.05, "Significant", "Not Significant")

sex_genes <- as.data.frame(sex_genes)
sex_genes <- apply(sex_genes,2,as.character)
sex_genes_1 <- as.data.frame(sex_genes) 
sum(sex_genes_1$sig_AGE == "Significant") 

write.csv(sex_genes, file="/Conj_AGE.csv")

##### Does the genotype influence the age-related expression?

### Adjust counts for covariates do not correct for age
assay(vsd_age) <- limma::removeBatchEffect(assay(vsd_age),batch = vsd_age$Braak_aSyn_stage, batch2 = NULL, covariates =  as.matrix(covariates),  model.matrix(~ vsd_age$Sex))

Adjusted_counts_age <- assay(vsd_age)

### Keep only interesting genes
Adjusted_counts_age <- Adjusted_counts_age[rownames(Adjusted_counts_age) %in% sex_genes$ensg,] 

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

Adjusted_counts_age_F[,101:117] <- lapply(Adjusted_counts_age_F[,101:117], factor)
Adjusted_counts_age_M[,101:117] <- lapply(Adjusted_counts_age_M[,101:117], factor)

### Make formula for each combination 
sex_genes$formula_age <- paste0(sex_genes$ensg , " ~ ","Age","*", sex_genes$Geno)

### Linear model
results_age_F <- lapply(sex_genes$formula_age, function(x) lm(x, data = Adjusted_counts_age_F))
results_age_M <- lapply(sex_genes$formula_age, function(x) lm(x, data = Adjusted_counts_age_M))


### Obtain slopes

slope_1 <- lstrends(results_age_F[[ 1 ]], "X1.205652472.G.A_A",var = "Age")
slope_2 <- lstrends(results_age_F[[ 2 ]], "X1.205652472.G.A_A",var = "Age")
slope_3 <- lstrends(results_age_F[[ 3 ]], "X1.205652472.G.A_A",var = "Age")
slope_4 <- lstrends(results_age_F[[ 4 ]], "X1.205652472.G.A_A",var = "Age")
slope_5 <- lstrends(results_age_F[[ 5 ]], "X1.205652472.G.A_A",var = "Age")
slope_6 <- lstrends(results_age_F[[ 6 ]], "X1.205652472.G.A_A",var = "Age")
slope_7 <- lstrends(results_age_F[[ 7 ]], "X1.205735891.G.A_A",var = "Age")
slope_8 <- lstrends(results_age_F[[ 8 ]], "X1.205735891.G.A_A",var = "Age")
slope_9 <- lstrends(results_age_F[[ 9 ]], "X1.205735891.G.A_A",var = "Age")
slope_10 <- lstrends(results_age_F[[ 10 ]], "X1.205735891.G.A_A",var = "Age")
slope_11 <- lstrends(results_age_F[[ 11 ]], "X1.205735891.G.A_A",var = "Age")
slope_12 <- lstrends(results_age_F[[ 12 ]], "X1.205735891.G.A_A",var = "Age")
slope_13 <- lstrends(results_age_F[[ 13 ]], "X1.205739266.C.T_C",var = "Age")
slope_14 <- lstrends(results_age_F[[ 14 ]], "X1.205739266.C.T_C",var = "Age")
slope_15 <- lstrends(results_age_F[[ 15 ]], "X1.205739266.C.T_C",var = "Age")
slope_16 <- lstrends(results_age_F[[ 16 ]], "X1.205739266.C.T_C",var = "Age")
slope_17 <- lstrends(results_age_F[[ 17 ]], "X1.205739266.C.T_C",var = "Age")
slope_18 <- lstrends(results_age_F[[ 18 ]], "X1.205739266.C.T_C",var = "Age")
slope_19 <- lstrends(results_age_F[[ 19 ]], "X1.243674772.A.G_G",var = "Age")
slope_20 <- lstrends(results_age_F[[ 20 ]], "X1.243674772.A.G_G",var = "Age")
slope_21 <- lstrends(results_age_F[[ 21 ]], "X1.243674772.A.G_G",var = "Age")
slope_22 <- lstrends(results_age_F[[ 22 ]], "X1.243717682.A.G_A",var = "Age")
slope_23 <- lstrends(results_age_F[[ 23 ]], "X1.243717682.A.G_A",var = "Age")
slope_24 <- lstrends(results_age_F[[ 24 ]], "X1.243717682.A.G_A",var = "Age")
slope_25 <- lstrends(results_age_F[[ 25 ]], "X1.243717682.A.G_A",var = "Age")
slope_26 <- lstrends(results_age_F[[ 26 ]], "X1.243717682.A.G_A",var = "Age")
slope_27 <- lstrends(results_age_F[[ 27 ]], "X1.243717682.A.G_A",var = "Age")
slope_28 <- lstrends(results_age_F[[ 28 ]], "X10.121462047.G.A_A",var = "Age")
slope_29 <- lstrends(results_age_F[[ 29 ]], "X10.121462047.G.A_A",var = "Age")
slope_30 <- lstrends(results_age_F[[ 30 ]], "X10.121708929.T.C_C",var = "Age")
slope_31 <- lstrends(results_age_F[[ 31 ]], "X10.121708929.T.C_C",var = "Age")
slope_32 <- lstrends(results_age_F[[ 32 ]], "X10.121708929.T.C_C",var = "Age")
slope_33 <- lstrends(results_age_F[[ 33 ]], "X10.121708929.T.C_C",var = "Age")
slope_34 <- lstrends(results_age_F[[ 34 ]], "X12.123404966.C.A_A",var = "Age")
slope_35 <- lstrends(results_age_F[[ 35 ]], "X12.123404966.C.A_A",var = "Age")
slope_36 <- lstrends(results_age_F[[ 36 ]], "X12.123404966.C.A_A",var = "Age")
slope_37 <- lstrends(results_age_F[[ 37 ]], "X12.123404966.C.A_A",var = "Age")
slope_38 <- lstrends(results_age_F[[ 38 ]], "X12.123404966.C.A_A",var = "Age")
slope_39 <- lstrends(results_age_F[[ 39 ]], "X12.123404966.C.A_A",var = "Age")
slope_40 <- lstrends(results_age_F[[ 40 ]], "X12.123404966.C.A_A",var = "Age")
slope_41 <- lstrends(results_age_F[[ 41 ]], "X12.30895251.T.C_C",var = "Age")
slope_42 <- lstrends(results_age_F[[ 42 ]], "X12.30895251.T.C_C",var = "Age")
slope_43 <- lstrends(results_age_F[[ 43 ]], "X12.30895251.T.C_C",var = "Age")
slope_44 <- lstrends(results_age_F[[ 44 ]], "X12.30895251.T.C_C",var = "Age")
slope_45 <- lstrends(results_age_F[[ 45 ]], "X12.30895251.T.C_C",var = "Age")
slope_46 <- lstrends(results_age_F[[ 46 ]], "X12.30895251.T.C_C",var = "Age")
slope_47 <- lstrends(results_age_F[[ 47 ]], "X17.16081958.C.T_T",var = "Age")
slope_48 <- lstrends(results_age_F[[ 48 ]], "X17.16081958.C.T_T",var = "Age")
slope_49 <- lstrends(results_age_F[[ 49 ]], "X17.16081958.C.T_T",var = "Age")
slope_50 <- lstrends(results_age_F[[ 50 ]], "X17.16081958.C.T_T",var = "Age")
slope_51 <- lstrends(results_age_F[[ 51 ]], "X17.16081958.C.T_T",var = "Age")
slope_52 <- lstrends(results_age_F[[ 52 ]], "X17.16081958.C.T_T",var = "Age")
slope_53 <- lstrends(results_age_F[[ 53 ]], "X17.16081958.C.T_T",var = "Age")
slope_54 <- lstrends(results_age_F[[ 54 ]], "X17.16081958.C.T_T",var = "Age")
slope_55 <- lstrends(results_age_F[[ 55 ]], "X17.16081958.C.T_T",var = "Age")
slope_56 <- lstrends(results_age_F[[ 56 ]], "X2.27677691.G.A_A",var = "Age")
slope_57 <- lstrends(results_age_F[[ 57 ]], "X2.27677691.G.A_A",var = "Age")
slope_58 <- lstrends(results_age_F[[ 58 ]], "X2.27677691.G.A_A",var = "Age")
slope_59 <- lstrends(results_age_F[[ 59 ]], "X2.27677691.G.A_A",var = "Age")
slope_60 <- lstrends(results_age_F[[ 60 ]], "X2.27677691.G.A_A",var = "Age")
slope_61 <- lstrends(results_age_F[[ 61 ]], "X2.27677691.G.A_A",var = "Age")
slope_62 <- lstrends(results_age_F[[ 62 ]], "X2.27677691.G.A_A",var = "Age")
slope_63 <- lstrends(results_age_F[[ 63 ]], "X2.27677691.G.A_A",var = "Age")
slope_64 <- lstrends(results_age_F[[ 64 ]], "X2.27677691.G.A_A",var = "Age")
slope_65 <- lstrends(results_age_F[[ 65 ]], "X2.27677691.G.A_A",var = "Age")
slope_66 <- lstrends(results_age_F[[ 66 ]], "X2.27677691.G.A_A",var = "Age")
slope_67 <- lstrends(results_age_F[[ 67 ]], "X2.27677691.G.A_A",var = "Age")
slope_68 <- lstrends(results_age_F[[ 68 ]], "X2.27677691.G.A_A",var = "Age")
slope_69 <- lstrends(results_age_F[[ 69 ]], "X2.27677691.G.A_A",var = "Age")
slope_70 <- lstrends(results_age_F[[ 70 ]], "X2.27677691.G.A_A",var = "Age")
slope_71 <- lstrends(results_age_F[[ 71 ]], "X2.27677691.G.A_A",var = "Age")
slope_72 <- lstrends(results_age_F[[ 72 ]], "X2.27677691.G.A_A",var = "Age")
slope_73 <- lstrends(results_age_F[[ 73 ]], "X2.27677691.G.A_A",var = "Age")
slope_74 <- lstrends(results_age_F[[ 74 ]], "X2.27677691.G.A_A",var = "Age")
slope_75 <- lstrends(results_age_F[[ 75 ]], "X2.27677691.G.A_A",var = "Age")
slope_76 <- lstrends(results_age_F[[ 76 ]], "X2.27677691.G.A_A",var = "Age")
slope_77 <- lstrends(results_age_F[[ 77 ]], "X2.27677691.G.A_A",var = "Age")
slope_78 <- lstrends(results_age_F[[ 78 ]], "X2.27677691.G.A_A",var = "Age")
slope_79 <- lstrends(results_age_F[[ 79 ]], "X2.27677691.G.A_A",var = "Age")
slope_80 <- lstrends(results_age_F[[ 80 ]], "X2.27677691.G.A_A",var = "Age")
slope_81 <- lstrends(results_age_F[[ 81 ]], "X2.27677691.G.A_A",var = "Age")
slope_82 <- lstrends(results_age_F[[ 82 ]], "X2.27677691.G.A_A",var = "Age")
slope_83 <- lstrends(results_age_F[[ 83 ]], "X2.27741237.T.C_T",var = "Age")
slope_84 <- lstrends(results_age_F[[ 84 ]], "X2.27741237.T.C_T",var = "Age")
slope_85 <- lstrends(results_age_F[[ 85 ]], "X2.27741237.T.C_T",var = "Age")
slope_86 <- lstrends(results_age_F[[ 86 ]], "X2.27741237.T.C_T",var = "Age")
slope_87 <- lstrends(results_age_F[[ 87 ]], "X2.27741237.T.C_T",var = "Age")
slope_88 <- lstrends(results_age_F[[ 88 ]], "X2.27741237.T.C_T",var = "Age")
slope_89 <- lstrends(results_age_F[[ 89 ]], "X2.27741237.T.C_T",var = "Age")
slope_90 <- lstrends(results_age_F[[ 90 ]], "X2.27741237.T.C_T",var = "Age")
slope_91 <- lstrends(results_age_F[[ 91 ]], "X3.100472848.T.C_C",var = "Age")
slope_92 <- lstrends(results_age_F[[ 92 ]], "X3.100472848.T.C_C",var = "Age")
slope_93 <- lstrends(results_age_F[[ 93 ]], "X3.100472848.T.C_C",var = "Age")
slope_94 <- lstrends(results_age_F[[ 94 ]], "X3.100472848.T.C_C",var = "Age")
slope_95 <- lstrends(results_age_F[[ 95 ]], "X3.100472848.T.C_C",var = "Age")
slope_96 <- lstrends(results_age_F[[ 96 ]], "X3.183980416.T.C_T",var = "Age")
slope_97 <- lstrends(results_age_F[[ 97 ]], "X3.183980416.T.C_T",var = "Age")
slope_98 <- lstrends(results_age_F[[ 98 ]], "X3.183980416.T.C_T",var = "Age")
slope_99 <- lstrends(results_age_F[[ 99 ]], "X3.183980416.T.C_T",var = "Age")
slope_100 <- lstrends(results_age_F[[ 100 ]], "X3.183980416.T.C_T",var = "Age")
slope_101 <- lstrends(results_age_F[[ 101 ]], "X3.183980416.T.C_T",var = "Age")
slope_102 <- lstrends(results_age_F[[ 102 ]], "X3.183980416.T.C_T",var = "Age")
slope_103 <- lstrends(results_age_F[[ 103 ]], "X3.183980416.T.C_T",var = "Age")
slope_104 <- lstrends(results_age_F[[ 104 ]], "X3.183980416.T.C_T",var = "Age")
slope_105 <- lstrends(results_age_F[[ 105 ]], "X3.183980416.T.C_T",var = "Age")
slope_106 <- lstrends(results_age_F[[ 106 ]], "X6.119143549.T.C_T",var = "Age")
slope_107 <- lstrends(results_age_F[[ 107 ]], "X6.119143549.T.C_T",var = "Age")
slope_108 <- lstrends(results_age_F[[ 108 ]], "X6.119143549.T.C_T",var = "Age")
slope_109 <- lstrends(results_age_F[[ 109 ]], "X6.32282854.A.G_G",var = "Age")
slope_110 <- lstrends(results_age_F[[ 110 ]], "X6.32282854.A.G_G",var = "Age")
slope_111 <- lstrends(results_age_F[[ 111 ]], "X6.32282854.A.G_G",var = "Age")
slope_112 <- lstrends(results_age_F[[ 112 ]], "X6.32577380.A.G_G",var = "Age")
slope_113 <- lstrends(results_age_F[[ 113 ]], "X6.32577380.A.G_G",var = "Age")
slope_114 <- lstrends(results_age_F[[ 114 ]], "X6.32577380.A.G_G",var = "Age")
slope_115 <- lstrends(results_age_F[[ 115 ]], "X6.32577380.A.G_G",var = "Age")
slope_116 <- lstrends(results_age_F[[ 116 ]], "X6.32577380.A.G_G",var = "Age")
slope_117 <- lstrends(results_age_F[[ 117 ]], "X6.32577380.A.G_G",var = "Age")


Slopes_F <- list(slope_1,slope_2,slope_3,slope_4,slope_5,slope_6,slope_7,slope_8,slope_9,slope_10,
                 slope_11,slope_12,slope_13,slope_14,slope_15,slope_16,slope_17,slope_18,slope_19,slope_20,
                 slope_21,slope_22,slope_23,slope_24,slope_25,slope_26,slope_27,slope_28,slope_29,slope_30,
                 slope_31,slope_32,slope_33,slope_34,slope_35,slope_36,slope_37,slope_38,slope_39,slope_40,
                 slope_41,slope_42,slope_43,slope_44,slope_45,slope_46,slope_47,slope_48,slope_49,slope_50,
                 slope_51,slope_52,slope_53,slope_54,slope_55,slope_56,slope_57,slope_58,slope_59,slope_60,
                 slope_61,slope_62,slope_63,slope_64,slope_65,slope_66,slope_67,slope_68,slope_69,slope_70,
                 slope_71,slope_72,slope_73,slope_74,slope_75,slope_76,slope_77,slope_78,slope_79,slope_80,
                 slope_81,slope_82,slope_83,slope_84,slope_85,slope_86,slope_87,slope_88,slope_89,slope_90,
                 slope_91,slope_92,slope_93,slope_94,slope_95,slope_96,slope_97,slope_98,slope_99,slope_100,
                 slope_101,slope_102,slope_103,slope_104,slope_105,slope_106,slope_107,slope_108,slope_109,
                 slope_110,slope_111,slope_112,slope_113,slope_114,slope_115,slope_116,slope_117)


### Compare slopes
results_age_F_F <- lapply(Slopes_F, pairs)
results_pvalues <- lapply(results_age_F_F, summary)

### Extract p-values
p_val_F <- lapply(results_pvalues, function(x)(x$p.value))
names(p_val_F) <- sex_genes$formula_age

p_values_F <- as.data.frame(do.call(rbind, p_val_F))
p_values_F$sig <- ifelse((p_values_F$V1 < 0.05|p_values_F$V2 < 0.05|p_values_F$V3 < 0.05), "Significant", "Not Significant")

##### In males

### Obtain slopes

slope_1_M <- lstrends(results_age_M[[ 1 ]], "X1.205652472.G.A_A",var = "Age")
slope_2_M <- lstrends(results_age_M[[ 2 ]], "X1.205652472.G.A_A",var = "Age")
slope_3_M <- lstrends(results_age_M[[ 3 ]], "X1.205652472.G.A_A",var = "Age")
slope_4_M <- lstrends(results_age_M[[ 4 ]], "X1.205652472.G.A_A",var = "Age")
slope_5_M <- lstrends(results_age_M[[ 5 ]], "X1.205652472.G.A_A",var = "Age")
slope_6_M <- lstrends(results_age_M[[ 6 ]], "X1.205652472.G.A_A",var = "Age")
slope_7_M <- lstrends(results_age_M[[ 7 ]], "X1.205735891.G.A_A",var = "Age")
slope_8_M <- lstrends(results_age_M[[ 8 ]], "X1.205735891.G.A_A",var = "Age")
slope_9_M <- lstrends(results_age_M[[ 9 ]], "X1.205735891.G.A_A",var = "Age")
slope_10_M <- lstrends(results_age_M[[ 10 ]], "X1.205735891.G.A_A",var = "Age")
slope_11_M <- lstrends(results_age_M[[ 11 ]], "X1.205735891.G.A_A",var = "Age")
slope_12_M <- lstrends(results_age_M[[ 12 ]], "X1.205735891.G.A_A",var = "Age")
slope_13_M <- lstrends(results_age_M[[ 13 ]], "X1.205739266.C.T_C",var = "Age")
slope_14_M <- lstrends(results_age_M[[ 14 ]], "X1.205739266.C.T_C",var = "Age")
slope_15_M <- lstrends(results_age_M[[ 15 ]], "X1.205739266.C.T_C",var = "Age")
slope_16_M <- lstrends(results_age_M[[ 16 ]], "X1.205739266.C.T_C",var = "Age")
slope_17_M <- lstrends(results_age_M[[ 17 ]], "X1.205739266.C.T_C",var = "Age")
slope_18_M <- lstrends(results_age_M[[ 18 ]], "X1.205739266.C.T_C",var = "Age")
slope_19_M <- lstrends(results_age_M[[ 19 ]], "X1.243674772.A.G_G",var = "Age")
slope_20_M <- lstrends(results_age_M[[ 20 ]], "X1.243674772.A.G_G",var = "Age")
slope_21_M <- lstrends(results_age_M[[ 21 ]], "X1.243674772.A.G_G",var = "Age")
slope_22_M <- lstrends(results_age_M[[ 22 ]], "X1.243717682.A.G_A",var = "Age")
slope_23_M <- lstrends(results_age_M[[ 23 ]], "X1.243717682.A.G_A",var = "Age")
slope_24_M <- lstrends(results_age_M[[ 24 ]], "X1.243717682.A.G_A",var = "Age")
slope_25_M <- lstrends(results_age_M[[ 25 ]], "X1.243717682.A.G_A",var = "Age")
slope_26_M <- lstrends(results_age_M[[ 26 ]], "X1.243717682.A.G_A",var = "Age")
slope_27_M <- lstrends(results_age_M[[ 27 ]], "X1.243717682.A.G_A",var = "Age")
slope_28_M <- lstrends(results_age_M[[ 28 ]], "X10.121462047.G.A_A",var = "Age")
slope_29_M <- lstrends(results_age_M[[ 29 ]], "X10.121462047.G.A_A",var = "Age")
slope_30_M <- lstrends(results_age_M[[ 30 ]], "X10.121708929.T.C_C",var = "Age")
slope_31_M <- lstrends(results_age_M[[ 31 ]], "X10.121708929.T.C_C",var = "Age")
slope_32_M <- lstrends(results_age_M[[ 32 ]], "X10.121708929.T.C_C",var = "Age")
slope_33_M <- lstrends(results_age_M[[ 33 ]], "X10.121708929.T.C_C",var = "Age")
slope_34_M <- lstrends(results_age_M[[ 34 ]], "X12.123404966.C.A_A",var = "Age")
slope_35_M <- lstrends(results_age_M[[ 35 ]], "X12.123404966.C.A_A",var = "Age")
slope_36_M <- lstrends(results_age_M[[ 36 ]], "X12.123404966.C.A_A",var = "Age")
slope_37_M <- lstrends(results_age_M[[ 37 ]], "X12.123404966.C.A_A",var = "Age")
slope_38_M <- lstrends(results_age_M[[ 38 ]], "X12.123404966.C.A_A",var = "Age")
slope_39_M <- lstrends(results_age_M[[ 39 ]], "X12.123404966.C.A_A",var = "Age")
slope_40_M <- lstrends(results_age_M[[ 40 ]], "X12.123404966.C.A_A",var = "Age")
slope_41_M <- lstrends(results_age_M[[ 41 ]], "X12.30895251.T.C_C",var = "Age")
slope_42_M <- lstrends(results_age_M[[ 42 ]], "X12.30895251.T.C_C",var = "Age")
slope_43_M <- lstrends(results_age_M[[ 43 ]], "X12.30895251.T.C_C",var = "Age")
slope_44_M <- lstrends(results_age_M[[ 44 ]], "X12.30895251.T.C_C",var = "Age")
slope_45_M <- lstrends(results_age_M[[ 45 ]], "X12.30895251.T.C_C",var = "Age")
slope_46_M <- lstrends(results_age_M[[ 46 ]], "X12.30895251.T.C_C",var = "Age")
slope_47_M <- lstrends(results_age_M[[ 47 ]], "X17.16081958.C.T_T",var = "Age")
slope_48_M <- lstrends(results_age_M[[ 48 ]], "X17.16081958.C.T_T",var = "Age")
slope_49_M <- lstrends(results_age_M[[ 49 ]], "X17.16081958.C.T_T",var = "Age")
slope_50_M <- lstrends(results_age_M[[ 50 ]], "X17.16081958.C.T_T",var = "Age")
slope_51_M <- lstrends(results_age_M[[ 51 ]], "X17.16081958.C.T_T",var = "Age")
slope_52_M <- lstrends(results_age_M[[ 52 ]], "X17.16081958.C.T_T",var = "Age")
slope_53_M <- lstrends(results_age_M[[ 53 ]], "X17.16081958.C.T_T",var = "Age")
slope_54_M <- lstrends(results_age_M[[ 54 ]], "X17.16081958.C.T_T",var = "Age")
slope_55_M <- lstrends(results_age_M[[ 55 ]], "X17.16081958.C.T_T",var = "Age")
slope_56_M <- lstrends(results_age_M[[ 56 ]], "X2.27677691.G.A_A",var = "Age")
slope_57_M <- lstrends(results_age_M[[ 57 ]], "X2.27677691.G.A_A",var = "Age")
slope_58_M <- lstrends(results_age_M[[ 58 ]], "X2.27677691.G.A_A",var = "Age")
slope_59_M <- lstrends(results_age_M[[ 59 ]], "X2.27677691.G.A_A",var = "Age")
slope_60_M <- lstrends(results_age_M[[ 60 ]], "X2.27677691.G.A_A",var = "Age")
slope_61_M <- lstrends(results_age_M[[ 61 ]], "X2.27677691.G.A_A",var = "Age")
slope_62_M <- lstrends(results_age_M[[ 62 ]], "X2.27677691.G.A_A",var = "Age")
slope_63_M <- lstrends(results_age_M[[ 63 ]], "X2.27677691.G.A_A",var = "Age")
slope_64_M <- lstrends(results_age_M[[ 64 ]], "X2.27677691.G.A_A",var = "Age")
slope_65_M <- lstrends(results_age_M[[ 65 ]], "X2.27677691.G.A_A",var = "Age")
slope_66_M <- lstrends(results_age_M[[ 66 ]], "X2.27677691.G.A_A",var = "Age")
slope_67_M <- lstrends(results_age_M[[ 67 ]], "X2.27677691.G.A_A",var = "Age")
slope_68_M <- lstrends(results_age_M[[ 68 ]], "X2.27677691.G.A_A",var = "Age")
slope_69_M <- lstrends(results_age_M[[ 69 ]], "X2.27677691.G.A_A",var = "Age")
slope_70_M <- lstrends(results_age_M[[ 70 ]], "X2.27677691.G.A_A",var = "Age")
slope_71_M <- lstrends(results_age_M[[ 71 ]], "X2.27677691.G.A_A",var = "Age")
slope_72_M <- lstrends(results_age_M[[ 72 ]], "X2.27677691.G.A_A",var = "Age")
slope_73_M <- lstrends(results_age_M[[ 73 ]], "X2.27677691.G.A_A",var = "Age")
slope_74_M <- lstrends(results_age_M[[ 74 ]], "X2.27677691.G.A_A",var = "Age")
slope_75_M <- lstrends(results_age_M[[ 75 ]], "X2.27677691.G.A_A",var = "Age")
slope_76_M <- lstrends(results_age_M[[ 76 ]], "X2.27677691.G.A_A",var = "Age")
slope_77_M <- lstrends(results_age_M[[ 77 ]], "X2.27677691.G.A_A",var = "Age")
slope_78_M <- lstrends(results_age_M[[ 78 ]], "X2.27677691.G.A_A",var = "Age")
slope_79_M <- lstrends(results_age_M[[ 79 ]], "X2.27677691.G.A_A",var = "Age")
slope_80_M <- lstrends(results_age_M[[ 80 ]], "X2.27677691.G.A_A",var = "Age")
slope_81_M <- lstrends(results_age_M[[ 81 ]], "X2.27677691.G.A_A",var = "Age")
slope_82_M <- lstrends(results_age_M[[ 82 ]], "X2.27677691.G.A_A",var = "Age")
slope_83_M <- lstrends(results_age_M[[ 83 ]], "X2.27741237.T.C_T",var = "Age")
slope_84_M <- lstrends(results_age_M[[ 84 ]], "X2.27741237.T.C_T",var = "Age")
slope_85_M <- lstrends(results_age_M[[ 85 ]], "X2.27741237.T.C_T",var = "Age")
slope_86_M <- lstrends(results_age_M[[ 86 ]], "X2.27741237.T.C_T",var = "Age")
slope_87_M <- lstrends(results_age_M[[ 87 ]], "X2.27741237.T.C_T",var = "Age")
slope_88_M <- lstrends(results_age_M[[ 88 ]], "X2.27741237.T.C_T",var = "Age")
slope_89_M <- lstrends(results_age_M[[ 89 ]], "X2.27741237.T.C_T",var = "Age")
slope_90_M <- lstrends(results_age_M[[ 90 ]], "X2.27741237.T.C_T",var = "Age")
slope_91_M <- lstrends(results_age_M[[ 91 ]], "X3.100472848.T.C_C",var = "Age")
slope_92_M <- lstrends(results_age_M[[ 92 ]], "X3.100472848.T.C_C",var = "Age")
slope_93_M <- lstrends(results_age_M[[ 93 ]], "X3.100472848.T.C_C",var = "Age")
slope_94_M <- lstrends(results_age_M[[ 94 ]], "X3.100472848.T.C_C",var = "Age")
slope_95_M <- lstrends(results_age_M[[ 95 ]], "X3.100472848.T.C_C",var = "Age")
slope_96_M <- lstrends(results_age_M[[ 96 ]], "X3.183980416.T.C_T",var = "Age")
slope_97_M <- lstrends(results_age_M[[ 97 ]], "X3.183980416.T.C_T",var = "Age")
slope_98_M <- lstrends(results_age_M[[ 98 ]], "X3.183980416.T.C_T",var = "Age")
slope_99_M <- lstrends(results_age_M[[ 99 ]], "X3.183980416.T.C_T",var = "Age")
slope_100_M <- lstrends(results_age_M[[ 100 ]], "X3.183980416.T.C_T",var = "Age")
slope_101_M <- lstrends(results_age_M[[ 101 ]], "X3.183980416.T.C_T",var = "Age")
slope_102_M <- lstrends(results_age_M[[ 102 ]], "X3.183980416.T.C_T",var = "Age")
slope_103_M <- lstrends(results_age_M[[ 103 ]], "X3.183980416.T.C_T",var = "Age")
slope_104_M <- lstrends(results_age_M[[ 104 ]], "X3.183980416.T.C_T",var = "Age")
slope_105_M <- lstrends(results_age_M[[ 105 ]], "X3.183980416.T.C_T",var = "Age")
slope_106_M <- lstrends(results_age_M[[ 106 ]], "X6.119143549.T.C_T",var = "Age")
slope_107_M <- lstrends(results_age_M[[ 107 ]], "X6.119143549.T.C_T",var = "Age")
slope_108_M <- lstrends(results_age_M[[ 108 ]], "X6.119143549.T.C_T",var = "Age")
slope_109_M <- lstrends(results_age_M[[ 109 ]], "X6.32282854.A.G_G",var = "Age")
slope_110_M <- lstrends(results_age_M[[ 110 ]], "X6.32282854.A.G_G",var = "Age")
slope_111_M <- lstrends(results_age_M[[ 111 ]], "X6.32282854.A.G_G",var = "Age")
slope_112_M <- lstrends(results_age_M[[ 112 ]], "X6.32577380.A.G_G",var = "Age")
slope_113_M <- lstrends(results_age_M[[ 113 ]], "X6.32577380.A.G_G",var = "Age")
slope_114_M <- lstrends(results_age_M[[ 114 ]], "X6.32577380.A.G_G",var = "Age")
slope_115_M <- lstrends(results_age_M[[ 115 ]], "X6.32577380.A.G_G",var = "Age")
slope_116_M <- lstrends(results_age_M[[ 116 ]], "X6.32577380.A.G_G",var = "Age")
slope_117_M <- lstrends(results_age_M[[ 117 ]], "X6.32577380.A.G_G",var = "Age")

Slopes_M <- list(slope_1_M,slope_2_M,slope_3_M,slope_4_M,slope_5_M,slope_6_M,slope_7_M,slope_8_M,slope_9_M,slope_10_M,
                 slope_11_M,slope_12_M,slope_13_M,slope_14_M,slope_15_M,slope_16_M,slope_17_M,slope_18_M,slope_19_M,
                 slope_20_M,slope_21_M,slope_22_M,slope_23_M,slope_24_M,slope_25_M,slope_26_M,slope_27_M,slope_28_M,
                 slope_29_M,slope_30_M,slope_31_M,slope_32_M,slope_33_M,slope_34_M,slope_35_M,slope_36_M,slope_37_M,
                 slope_38_M,slope_39_M,slope_40_M,slope_41_M,slope_42_M,slope_43_M,slope_44_M,slope_45_M,slope_46_M,
                 slope_47_M,slope_48_M,slope_49_M,slope_50_M,slope_51_M,slope_52_M,slope_53_M,slope_54_M,slope_55_M,
                 slope_56_M,slope_57_M,slope_58_M,slope_59_M,slope_60_M,slope_61_M,slope_62_M,slope_63_M,slope_64_M,
                 slope_65_M,slope_66_M,slope_67_M,slope_68_M,slope_69_M,slope_70_M,slope_71_M,slope_72_M,slope_73_M,
                 slope_74_M,slope_75_M,slope_76_M,slope_77_M,slope_78_M,slope_79_M,slope_80_M,slope_81_M,slope_82_M,
                 slope_83_M,slope_84_M,slope_85_M,slope_86_M,slope_87_M,slope_88_M,slope_89_M,slope_90_M,slope_91_M,
                 slope_92_M,slope_93_M,slope_94_M,slope_95_M,slope_96_M,slope_97_M,slope_98_M,slope_99_M,slope_100_M,
                 slope_101_M,slope_102_M,slope_103_M,slope_104_M,slope_105_M,slope_106_M,slope_107_M,slope_108_M,slope_109_M,
                 slope_110_M,slope_111_M,slope_112_M,slope_113_M,slope_114_M,slope_115_M,slope_116_M,slope_117_M)

### Compare slopes
results_age_M_F <- lapply(Slopes_M, pairs)
results_pvalues_M <- lapply(results_age_M_F, summary)

### Extract p-values
p_val_M <- lapply(results_pvalues_M, function(x)(x$p.value))
names(p_val_M) <- sex_genes$formula

p_values_M <- as.data.frame(do.call(rbind, p_val_M))
p_values_M$sig <- ifelse((p_values_M$V1 < 0.05|p_values_M$V2 < 0.05|p_values_M$V3 < 0.05), "Significant", "Not Significant")

Age_results <- cbind(p_values_F,p_values_M)
colnames(Age_results) <- c("0vs1_F", "0vs2_F","1vs2_F", "Sig_F", "0vs1_M", "0vs2_M","1vs2_M", "Sig_M")

sex_genes <- cbind(sex_genes,Age_results)

write.csv(as.data.frame(sex_genes), file="/Conj_Geno_Age.csv",  row.names = TRUE)

