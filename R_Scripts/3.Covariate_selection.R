###   Covariate Selection


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
library("ComplexHeatmap")

### Load dds
load("./dds_84.Rda")
coldata <- colData(dds)

# Subset control group
coldata_Ctrl <- coldata[coldata$Neuropath_diagnosis == "CONTR", ]
keep_Ctrl <- rownames(coldata_Ctrl)

dds_Ctrl <- dds[,keep_Ctrl]
dds_Ctrl$Neuropath_diagnosis <- NULL
coldata_Ctrl <- as.data.frame(coldata_Ctrl)
coldata_Ctrl$Sex <- as.factor(coldata_Ctrl$Sex)
dds_Ctrl$Sex <- coldata_Ctrl$Sex

design(dds_Ctrl) <- ~ Age_death + Sex



#####################################################################################################################
### qSVA: correct for RNA quality
#####################################################################################################################

### Generate degradation matrix

# Import degradation matrix to get the sample names from the first column
degMat <- read.delim("/degradation.matrix_84.txt") 
degMat_Ctrl <- degMat[degMat$X %in% rownames(coldata_Ctrl),]

write.table(degMat_Ctrl, file = "/degradation.matrix_Ctrl.txt", sep = "\t",  row.names = FALSE, quote = FALSE)


SampleNames_Ctrl <- degMat_Ctrl$X
SampleNames_Ctrl

# Save path to degradation matrix
covFile_Ctrl <- "/degradation.matrix_Ctrl.txt"  

# Import table with number of reads per sample and create a vector
Mreads <- read.csv("/Volumes/Data/RNA-Seq/Results/Final_09.2021/Braak_Stage_Groups/Additional_Files/MReads_84.csv")
Mreads_Ctrl <- Mreads[Mreads$Sample.Name %in% rownames(coldata_Ctrl),]
Reads_Ctrl <- Mreads_Ctrl$M.Seqs

# Generate degradation matrix
degCovAdj2_Ctrl = read.degradation.matrix(
  covFiles = covFile_Ctrl, # coverage file(s) for degradation regions 		-> path to matrix saved as a txt (tab separated file)
  sampleNames = SampleNames_Ctrl, # sample names; creates column names of degradation matrix -> rownames of degradation matrix
  readLength = 150, # read length in base pairs
  totalMapped = Reads_Ctrl, # how many reads per sample (library size normalization)  -> you can get the total number of reads either from the first line of the flagstats files or with the following command: samtools view -c /Volumes/Data/RNA-Seq/Results/HISAT2/LINUX/sample_X.bam
  normFactor = 80,  # common library size to normalize to, 80M as default
  type="region_matrix_all") # whether input are individual 'bwtool' output, 'region_matrix' run on individual samples, or 'region_matrix' run on all samples together


### Estimate qSVs
qSVs_Ctrl <- qsva(degCovAdj2_Ctrl)
coldata_qSVs_Ctrl <- merge(coldata_Ctrl, qSVs_Ctrl, by.x = "row.names", by.y = "row.names")

### add column with qSVs to coldata
dds_Ctrl$qSV1 <- coldata_qSVs_Ctrl$PC1
dds_Ctrl$qSV2 <- coldata_qSVs_Ctrl$PC2

design(dds_Ctrl) <- ~ Age_death + qSV1 + qSV2  + Sex
coldata_Ctrl <- colData(dds_Ctrl)

save(dds_Ctrl, file = "/dds_Ctrl.Rda")
 

## Load cell composition

load("/cells_Ctrl.Rda")

coldata_Ctrl$Neurons <- cell_1$Neurons
coldata_Ctrl$Astrocytes <- cell_1$Astrocytes
coldata_Ctrl$VLMC <- cell_1$VLMC
coldata_Ctrl$Tanycytes <- cell_1$Tanycytes
coldata_Ctrl$OPC <- cell_1$OPC
coldata_Ctrl$Oligodendrocytes <- cell_1$Oligodendrocytes
coldata_Ctrl$Ependymal <- cell_1$Ependymal
coldata_Ctrl$Unknown <- cell_1$Unknown
coldata_Ctrl$Endothelial <- cell_1$Endothelial
coldata_Ctrl$Microglia <- cell_1$Microglia
coldata_Ctrl$NFO <- cell_1$NFO


coldata_Ctrl <- as.data.frame(coldata_Ctrl)
coldata_Ctrl[,c(5:8)] <- NULL
coldata_Ctrl[,c(5:8)] <- NULL

coldata_Ctrl <- coldata_Ctrl %>% mutate(Sex =
                          case_when(Sex == "F" ~ "1",
                                    Sex == "M" ~ "2")
)

coldata_Ctrl <- mutate_all(coldata_Ctrl, function(x) as.numeric(as.character(x)))


colnames(coldata_Ctrl)
coldata_Ctrl <- coldata_Ctrl[, c( "Sex","Age_death","PMD_min","RIN",
                             "Neurons","Astrocytes"  , "VLMC" ,  "Tanycytes" , "OPC" ,"Oligodendrocytes","Ependymal",   "Unknown" ,  "Endothelial" , "Microglia", "NFO",
                             "qSV1" ,  "qSV2"  )]


######Pearson correlation
Pearson_matrix_Ctrl <-  cor_mat(coldata_Ctrl, method = "pearson")

# get the p.values
Pearson_matrix_p_Ctrl <- Pearson_matrix_Ctrl %>% cor_get_pval()

Pearson_matrix_Ctrl <- as.data.frame(Pearson_matrix_Ctrl)
rownames(Pearson_matrix_Ctrl) <- Pearson_matrix_Ctrl$rowname
Pearson_matrix_Ctrl$rowname <- NULL
Pearson_matrix_Ctrl <- mutate_all(Pearson_matrix_Ctrl, function(x) as.numeric(as.character(x)))
Pearson_matrix_Ctrl <- as.matrix(Pearson_matrix_Ctrl)

Pearson_matrix_p_Ctrl <- as.data.frame(Pearson_matrix_p_Ctrl)
rownames(Pearson_matrix_p_Ctrl) <- Pearson_matrix_p_Ctrl$rowname
Pearson_matrix_p_Ctrl$rowname <- NULL
Pearson_matrix_p_Ctrl <- mutate_all(Pearson_matrix_p_Ctrl, function(x) as.numeric(as.character(x)))
Pearson_matrix_p_Ctrl <- as.matrix(Pearson_matrix_p_Ctrl)

### Plot tau coefficients
# Initialize file path

cor_plot_Ctrl <- { 
  corrplot(Pearson_matrix_Ctrl, p.mat = Pearson_matrix_p_Ctrl, sig.level = 0.05, type = "lower", diag = FALSE, 
           tl.col = "black", tl.srt = 45, insig='blank', tl.cex = 1)
  recordPlot()
}

cor_plot_Ctrl

### Plot p-values
Pearson_matrix_p_log_Ctrl <- -log10(Pearson_matrix_p_Ctrl)
Pearson_matrix_p_log_Ctrl

### Plot known variables against qSVs
Pearson_matrix_p_log_plot_Ctrl <- Pearson_matrix_p_log_Ctrl[c(4,3,5:15),16:17]


max(Pearson_matrix_p_log_plot_Ctrl)
col_fun = colorRamp2(c(0, 1.3, 4), c("grey", "white", "red2"))

png("/P-value_cor_qSVs_Ctrl.png", width=5,height= 6,units="in",res=600)
Hm_Ctrl <- Heatmap(Pearson_matrix_p_log_plot_Ctrl, col = col_fun, name = " ", cluster_columns = FALSE, cluster_rows = FALSE, column_names_rot = 0, column_names_centered = 1, row_names_side = "left", rect_gp = gpar(col = "white", lwd = 1), width = unit(3, "in"), 
              height = unit(5, "in"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(Pearson_matrix_p_log_plot_Ctrl[i, j] > 1.29)
                  grid.text(sprintf("%.1f", Pearson_matrix_p_log_plot_Ctrl[i, j]), x, y, gp = gpar(fontsize = 12))
              })
draw(Hm_Ctrl)
Hm_plot_Ctrl <- recordPlot()
dev.off()




# Subset PD group (Braak Lewy body stage 5)
coldata_PD <- coldata[coldata$Braak_aSyn_stage == "5", ]
keep_PD <- rownames(coldata_PD)

dds_PD <- dds[,keep_PD]
coldata_PD <- as.data.frame(coldata_PD)
coldata_PD$Sex <- as.factor(coldata_PD$Sex)
dds_PD$Sex <- coldata_PD$Sex
design(dds_PD) <- ~ Age_death + Sex



#####################################################################################################################
### qSVA: correct for RNA quality
#####################################################################################################################

### Generate degradation matrix

# Import degradation matrix to get the sample names from the first column
degMat_PD <- degMat[degMat$X %in% rownames(coldata_PD),]

write.table(degMat_PD, file = "/degradation.matrix_5.txt", sep = "\t", row.names = FALSE, quote = FALSE)

SampleNames_PD <- degMat_PD$X


# Save path to degradation matrix
covFile_PD <- "/degradation.matrix_5.txt"  

# Import table with number of reads per sample and create a vector
Mreads_PD <- Mreads[Mreads$Sample.Name %in% rownames(coldata_PD),]
Reads_PD <- Mreads_PD$M.Seqs

# Generate degradation matrix
degCovAdj2_PD = read.degradation.matrix(
  covFiles = covFile_PD, # coverage file(s) for degradation regions 		-> path to matrix saved as a txt (tab separated file)
  sampleNames = SampleNames_PD, # sample names; creates column names of degradation matrix -> rownames of degradation matrix
  readLength = 150, # read length in base pairs
  totalMapped = Reads_PD, # how many reads per sample (library size normalization)  -> you can get the total number of reads either from the first line of the flagstats files or with the following command: samtools view -c /Volumes/Data/RNA-Seq/Results/HISAT2/LINUX/sample_X.bam
  normFactor = 80,  # common library size to normalize to, 80M as default
  type="region_matrix_all") # whether input are individual 'bwtool' output, 'region_PDatrix' run on individual samples, or 'region_PDatrix' run on all samples together


### Estimate qSVs
qSVs_PD <- qsva(degCovAdj2_PD)
coldata_qSVs_PD <- merge(coldata_PD, qSVs_PD, by.x = "row.names", by.y = "row.names")


### add colums with qSVs to coldata
dds_PD$qSV1 <- coldata_qSVs_PD$PC1
dds_PD$qSV2 <- coldata_qSVs_PD$PC2
dds_PD$qSV3 <- coldata_qSVs_PD$PC3
dds_PD$qSV4 <- coldata_qSVs_PD$PC4


coldata_PD <- colData(dds_PD)

load("/cells_5.Rda")

coldata_PD$Neurons <- cell_1$Neurons
coldata_PD$Astrocytes <- cell_1$Astrocytes
coldata_PD$VLMC <- cell_1$VLMC
coldata_PD$Tanycytes <- cell_1$Tanycytes
coldata_PD$OPC <- cell_1$OPC
coldata_PD$Oligodendrocytes <- cell_1$Oligodendrocytes
coldata_PD$Ependymal <- cell_1$Ependymal
coldata_PD$Unknown <- cell_1$Unknown
coldata_PD$Endothelial <- cell_1$Endothelial
coldata_PD$Microglia <- cell_1$Microglia
coldata_PD$NFO <- cell_1$NFO


coldata_PD <- as.data.frame(coldata_PD)
coldata_PD[,c(5:8)] <- NULL
coldata_PD[,c(6:9)] <- NULL

coldata_PD <- coldata_PD %>% mutate(Sex =
                                      case_when(Sex == "F" ~ "1",
                                                Sex == "M" ~ "2")
)

coldata_PD <- coldata_PD %>% mutate(Diagnosis =
                                      case_when(Neuropath_diagnosis == "iLBD" ~ "1",
                                                Neuropath_diagnosis == "PD" ~ "2",
                                                Neuropath_diagnosis == "PDD" ~ "3")
)

dds_PD$Diagnosis <- coldata_PD$Diagnosis
save(dds_PD, file = "/dds_5.Rda")

coldata_PD$Neuropath_diagnosis <- NULL
coldata_PD <- mutate_all(coldata_PD, function(x) as.numeric(as.character(x)))


colnames(coldata_PD)
coldata_PD <- coldata_PD[, c( "Sex","Diagnosis","Age_death","PMD_min","RIN",
                              "Neurons","Astrocytes"  , "VLMC" ,  "Tanycytes" , "OPC" ,"Oligodendrocytes","Ependymal",   "Unknown" ,  "Endothelial" , "Microglia", "NFO",
                              "qSV1" ,  "qSV2", "qSV3" , "qSV4" )]


######Pearson correlation
Pearson_matrix_PD <-  cor_mat(coldata_PD, method = "pearson")

# get the p.values
Pearson_matrix_p_PD <- Pearson_matrix_PD %>% cor_get_pval()

Pearson_matrix_PD <- as.data.frame(Pearson_matrix_PD)
rownames(Pearson_matrix_PD) <- Pearson_matrix_PD$rowname
Pearson_matrix_PD$rowname <- NULL
Pearson_matrix_PD <- mutate_all(Pearson_matrix_PD, function(x) as.numeric(as.character(x)))
Pearson_matrix_PD <- as.matrix(Pearson_matrix_PD)

Pearson_matrix_p_PD <- as.data.frame(Pearson_matrix_p_PD)
rownames(Pearson_matrix_p_PD) <- Pearson_matrix_p_PD$rowname
Pearson_matrix_p_PD$rowname <- NULL
Pearson_matrix_p_PD <- mutate_all(Pearson_matrix_p_PD, function(x) as.numeric(as.character(x)))
Pearson_matrix_p_PD <- as.matrix(Pearson_matrix_p_PD)

### Plot tau coefficients
# Initialize file path

cor_plot_PD <- { 
  corrplot(Pearson_matrix_PD, p.mat = Pearson_matrix_p_PD, sig.level = 0.05, type = "lower", diag = FALSE, 
           tl.col = "black", tl.srt = 45, insig='blank', tl.cex = 1)
  recordPlot()
}

cor_plot_PD

### Plot p-values
Pearson_matrix_p_log_PD <- -log10(Pearson_matrix_p_PD)
Pearson_matrix_p_log_PD

### Plot known variables against qSVs
Pearson_matrix_p_log_plot_PD <- Pearson_matrix_p_log_PD[c(2,4,3,5:16),17:20]


max(Pearson_matrix_p_log_plot_PD)
col_fun = colorRamp2(c(0, 1.3, 2), c("grey", "white", "red2"))

png("/P-value_cor_qSVs_5.png", width=5,height= 6,units="in",res=600)
Hm_PD <- Heatmap(Pearson_matrix_p_log_plot_PD, col = col_fun, name = " ", cluster_columns = FALSE, cluster_rows = FALSE, column_names_rot = 0, column_names_centered = 1, row_names_side = "left", rect_gp = gpar(col = "white", lwd = 1), width = unit(3, "in"), 
                 height = unit(5, "in"),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(Pearson_matrix_p_log_plot_PD[i, j] > 1.29)
                     grid.text(sprintf("%.1f", Pearson_matrix_p_log_plot_PD[i, j]), x, y, gp = gpar(fontsize = 12))
                 })
draw(Hm_PD)
Hm_plot_PD <- recordPlot()
dev.off()
