### ### ### ### ### ### ### ### ### ### ### 
###   Differential expression analysis  ###
### ### ### ### ### ### ### ### ### ### ### 


### Packages required:
library("DESeq2")
library("tibble")
library("DEGreport")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("tidyr")
library("ensembldb")


#### In neuropathologically healthy controls (Braak LB stage 0)

load("/dds_Ctrl.Rda")
coldata_Ctrl <- as.data.frame(colData(dds_Ctrl))


# Check the design formula
design(dds_Ctrl)
# ~Age_death + qSV1 + qSV2 + Sex

dds_Ctrl$Sex <- relevel(dds_Ctrl$Sex, ref = "M")


#########################
# Wald test #
#########################

dds_Ctrl <- DESeq(dds_Ctrl)

### Save results
### Extract results for specific comparisons

res_Ctrl <- results(dds_Ctrl)
summary(res_Ctrl)
sum(res_Ctrl$padj < 0.1, na.rm = TRUE) 

### subset significant genes 
res_Ctrl_0.05 <- results(dds_Ctrl, alpha=0.05)
summary(res_Ctrl_0.05)
sum(res_Ctrl_0.05$padj < 0.05, na.rm = TRUE)


### Create a tibble for  results
res_Ctrl_tb <- res_Ctrl %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

### Save file with all genes
write.csv(as.data.frame(res_Ctrl_tb), file="/res_Ctrl_genes.csv",  row.names = FALSE)

 
### Subset to return genes with padj < 0.1
 res_Ctrl_0.1 <- res_Ctrl_tb %>% 
   dplyr::filter(padj < 0.1)
 
 
### Create a tibble for LRT results
res_Ctrl_0.05_tb <- res_Ctrl_0.05 %>%
   data.frame() %>%
   rownames_to_column(var="gene") %>% 
   as_tibble()
 
 
### Subset to return genes with padj < 0.05
res_Ctrl_0.05 <- res_Ctrl_0.05_tb %>% 
   dplyr::filter(padj < 0.05)
 
### Reorder based on adjusted p-value
res_Ctrl_0.05 <- res_Ctrl_0.05[order(res_Ctrl_0.05$padj),]
res_Ctrl_0.1 <- res_Ctrl_0.1[order(res_Ctrl_0.1$padj),]
 
### Save file with significant genes
write.csv(as.data.frame(res_Ctrl_0.05), file="/res_Ctrl_0.05.csv",  row.names = FALSE)
write.csv(as.data.frame(res_Ctrl_0.1), file="/res_Ctrl_0.1.csv",  row.names = FALSE)
 

#### In PD patients (Braak LB stage 5)

load("/dds_5.Rda")
coldata_PD <- as.data.frame(colData(dds_PD))

# Check the desing formula
design(dds_PD)  <- ~Age_death + qSV1 + qSV2 + qSV3 + qSV4 + Sex
# ~Age_death + qSV1 + qSV2 + qSV3 + qSV4 + Sex

dds_PD$Sex <- relevel(dds_PD$Sex, ref = "M")

#########################
# Wald test #
#########################

dds_PD <- DESeq(dds_PD)

### Save results
### Extract results for specific comparisons

res_PD <- results(dds_PD)
summary(res_PD)
sum(res_PD$padj < 0.1, na.rm = TRUE)

### subset significant genes 
res_PD_0.05 <- results(dds_PD, alpha=0.05)
summary(res_PD_0.05)
sum(res_PD_0.05$padj < 0.05, na.rm = TRUE) 

### Create a tibble for  results
res_PD_tb <- res_PD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

### Save file with all genes
write.csv(as.data.frame(res_PD_tb), file="/res_5_genes.csv",  row.names = FALSE)


### Subset to return genes with padj < 0.1
res_PD_0.1 <- res_PD_tb %>% 
  dplyr::filter(padj < 0.1)


### Create a tibbl  e for LRT results
res_PD_0.05_tb <- res_PD_0.05 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


### Subset to return genes with padj < 0.05
res_PD_0.05 <- res_PD_0.05_tb %>% 
  dplyr::filter(padj < 0.05)

### Reorder based on adjusted p-value
res_PD_0.05 <- res_PD_0.05[order(res_PD_0.05$padj),]
res_PD_0.1 <- res_PD_0.1[order(res_PD_0.1$padj),]

### Save file with significant genes
write.csv(as.data.frame(res_PD_0.05), file="/res_5_0.05.csv",  row.names = FALSE)
write.csv(as.data.frame(res_PD_0.1), file="/res_5_0.1.csv",  row.names = FALSE)

