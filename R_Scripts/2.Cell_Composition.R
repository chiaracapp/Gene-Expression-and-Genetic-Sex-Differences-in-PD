###   Cell Composition

### Packages required:
library("rstatix")
library("dplyr")
library("reshape2") 
library("DESeq2")
library("ggplot2")
library("biomaRt")
library("markerGeneProfile")
library("homologene")
library("ggpubr")
library("corrplot")

### Load dds
load("./dds_84.Rda")
coldata <- colData(dds)


#####################################################################################################################
### Cell deconvolution using Scaden 
#####################################################################################################################


### Steps in Scaden:
###   1. Process: Scaden performs pre-processing on your training data, making sure it has the same genes as your prediction (bulk) data and performing some data transformations to make the data suitable for machine learning.
###               For this step you need the training data downloaded (mouse brain) and the bulk data in the form of counts normalized by library size.
###               Scaden will create a new file for training which only contains the intersection of genes between the training and the prediction data. Furthermore, the training data will be log2-transformed and scaled to the range [0,1].
###   2. Train:  Scaden consists of three deep neural network models. By default, each of them will be trained for 5,000 steps, which is the recommended number of training steps.
###   3. Predict: Scaden predicts the cell composition


### Comparison in cell composition between groups

cell <- read.delim("/scaden_predictions.txt", header = TRUE, sep = "\t", dec = ".")

# add columns with braak stage and neuropath
cell <- merge(cell, coldata, by.x = "X", by.y = "row.names")
cell <- as.data.frame(cell)
# remove unwanted columns
cell <- cell[,-c(17:20,22,25)]
rownames(cell) <- cell[,"X"]
cell[,"X"] <- NULL

cell <- cell %>% mutate(Sex =
                                case_when(Sex == "F" ~ "1", 
                                          Sex == "M" ~ "2")
)

cell <- cell %>% mutate(Diagnosis =
                                case_when(Neuropath_diagnosis == "CONTR" ~ "0", 
                                          Neuropath_diagnosis == "iLBD" ~ "1",
                                          Neuropath_diagnosis == "PD" ~ "2",
                                          Neuropath_diagnosis == "PDD" ~ "3")
)

cell$Neuropath_diagnosis <- NULL

cell <- mutate_all(cell, function(x) as.numeric(as.character(x)))
str(cell)

cell <- cell %>% 
      rename(
            Age = Age_death,
            PMD = PMD_min,
            Group = Braak_aSyn_groups,
            Braak_Stage = Braak_aSyn_stage
  )

cell <- cell[, c("Group","Diagnosis", "Braak_Stage", "Sex", "Age", "PMD","RIN",  "qSV1" ,  "qSV2" , "qSV3" , "qSV4", "qSV5"  ,
                       "Neurons","Astrocytes"  , "VLMC" ,  "Tanycytes" , "OPC" ,"Oligodendrocytes","Ependymal",   "Unknown" ,  "Endothelial" , "Microglia", "NFO")]

str(cell)

cell <- mutate_all(cell, function(x) as.numeric(as.character(x)))

###############
#### Ctrl group

cell_1 <- cell[,-c(1,3, 8:12)]
cell_1 <- cell_1[!(cell_1$Diagnosis == "2" | cell_1$Diagnosis == "1" | cell_1$Diagnosis == "3"),]
cell_1$Diagnosis <- NULL

save(cell_1, file = "/cells_Ctrl.Rda")

### Pearson correlation
Pearson_matrix <-  cor_mat(cell_1, method = "pearson")

# get the p.values
Pearson_matrix_p <- Pearson_matrix %>% cor_get_pval()

# Arrange for plotting
Pearson_matrix <- as.data.frame(Pearson_matrix)
rownames(Pearson_matrix) <- Pearson_matrix$rowname
Pearson_matrix$rowname <- NULL
Pearson_matrix <- mutate_all(Pearson_matrix, function(x) as.numeric(as.character(x)))
Pearson_matrix <- as.matrix(Pearson_matrix)

Pearson_matrix_p <- as.data.frame(Pearson_matrix_p)
rownames(Pearson_matrix_p) <- Pearson_matrix_p$rowname
Pearson_matrix_p$rowname <- NULL
Pearson_matrix_p <- mutate_all(Pearson_matrix_p, function(x) as.numeric(as.character(x)))
Pearson_matrix_p <- as.matrix(Pearson_matrix_p)


### Plot tau coefficients
# Initialize file path
png("/Cell_pearson_Ctrl.png",width=4.5,height=6,units="in",res=600)

  corrplot(Pearson_matrix, p.mat = Pearson_matrix_p, sig.level = 0.05, type = "lower", diag = FALSE, 
           tl.col = "black", tl.srt = 45, insig='blank', tl.cex = 0.7, mar = c(0.5,0.5,0.5,0), cl.cex = 0.7)
Plot <- recordPlot()
print(Plot)
dev.off()


### Prepare dataset for ancova

cell_2 <- cell_1
cell_2$Sex <- as.factor(cell_2$Sex)

######## ANCOVA #############

res.Neurons <- cell_2 %>% anova_test(Neurons ~ Age + PMD + RIN + Sex)
get_anova_table(res.Neurons)

res.Astrocytes <- cell_2 %>% anova_test(Astrocytes ~ Age + PMD + RIN + Sex)
get_anova_table(res.Astrocytes)

res.VLMC <- cell_2 %>% anova_test(VLMC ~ Age + PMD + RIN + Sex)
get_anova_table(res.VLMC)

res.Tanycytes <- cell_2 %>% anova_test(Tanycytes ~ Age + PMD + RIN + Sex)
get_anova_table(res.Tanycytes)

res.OPC <- cell_2 %>% anova_test(OPC ~ Age + PMD + RIN + Sex)
get_anova_table(res.OPC)

res.Oligodendrocytes <- cell_2 %>% anova_test(Oligodendrocytes ~ Age + PMD + RIN + Sex)
get_anova_table(res.Oligodendrocytes)

res.Ependymal <- cell_2 %>% anova_test(Ependymal ~ Age + PMD + RIN + Sex)
get_anova_table(res.Ependymal)

res.Unknown <- cell_2 %>% anova_test(Unknown ~ Age + PMD + RIN + Sex)
get_anova_table(res.Unknown)

res.Endothelial <- cell_2 %>% anova_test(Endothelial ~ Age + PMD + RIN + Sex)
get_anova_table(res.Endothelial)

res.Microglia <- cell_2 %>% anova_test(Microglia ~ Age + PMD + RIN + Sex)
get_anova_table(res.Microglia)


p <- c(0.014, 0.011, 1, 1, 0.005, 0.011, 1, 0.003, 1, 1)
padj <- p.adjust(p, n= length(p))




### Plot cell composition 
# reshape the dataset
cell_3 <- cell_2
str(cell_3)
cell_3$Age <- NULL
cell_3$PMD <- NULL
cell_3$RIN <- NULL

cell_plot <- aggregate(. ~ Sex, cell_3, function(x) list(mean = round(mean(x),3), sd =round(sd(x),3)))
cell_plot <- do.call("data.frame", cell_plot) # flatten

Neurons <- cell_plot[,1:3]
Neurons$cell_type <- c(rep("Neurons", 2))
names(Neurons) <- c("Sex", "Mean", "SD", "Cell_Type")
Astrocytes <- cell_plot[,c(1,4,5)]
Astrocytes$cell_type <- c(rep("Astrocytes", 2))
names(Astrocytes) <- c("Sex", "Mean", "SD", "Cell_Type")
VLMC<- cell_plot[,c(1,6,7)]
VLMC$cell_type <- c(rep("VLMC", 2))
names(VLMC) <- c("Sex", "Mean", "SD", "Cell_Type")
Tanycytes <- cell_plot[,c(1,8,9)]
Tanycytes$cell_type <- c(rep("Tanycytes", 2))
names(Tanycytes) <- c("Sex", "Mean", "SD", "Cell_Type")
OPC <- cell_plot[,c(1,10,11)]
OPC$cell_type <- c(rep("OPC", 2))
names(OPC) <- c("Sex", "Mean", "SD", "Cell_Type")
Oligodendrocytes <- cell_plot[,c(1,12,13)]
Oligodendrocytes$cell_type <- c(rep("Oligodendrocytes", 2))
names(Oligodendrocytes) <- c("Sex", "Mean", "SD", "Cell_Type")
Ependymal <- cell_plot[,c(1,14,15)]
Ependymal$cell_type <- c(rep("Ependymal", 2))
names(Ependymal) <- c("Sex", "Mean", "SD", "Cell_Type")
Unknown <- cell_plot[,c(1,16,17)]
Unknown$cell_type <- c(rep("Unknown", 2))
names(Unknown) <- c("Sex", "Mean", "SD", "Cell_Type")
Endothelial <- cell_plot[,c(1,18,19)]
Endothelial$cell_type <- c(rep("Endothelial", 2))
names(Endothelial) <- c("Sex", "Mean", "SD", "Cell_Type")
Microglia <- cell_plot[,c(1,20,21)]
Microglia$cell_type <- c(rep("Microglia", 2))
names(Microglia) <- c("Sex", "Mean", "SD", "Cell_Type")
NFO <- cell_plot[,c(1,20,21)]
NFO$cell_type <- c(rep("NFO", 2))
names(NFO) <- c("Sex", "Mean", "SD", "Cell_Type")


cellplot <- rbind(Neurons,Astrocytes,VLMC,Tanycytes,OPC,Oligodendrocytes,Ependymal,Unknown,Endothelial,Microglia,NFO)

cellplot <- cellplot %>% mutate(Sex =
                                case_when(Sex == "1" ~ "F", 
                                          Sex == "2" ~ "M" )
)

cellplot$Mean <- as.numeric(cellplot$Mean)
cellplot$SD <- as.numeric(cellplot$SD)


my_colors_ID <- c("F" = "#c3cece",  "M" = "#98c5df")


Plot1 <- ggplot(cellplot, aes(x= reorder(Cell_Type, Mean),y=Mean,fill=factor(Sex)))+
            geom_bar(stat="identity",position="dodge",color = "black") +
            theme_minimal() +
            scale_fill_manual(values = my_colors_ID) + 
            scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
            geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,position=position_dodge(.9)) +
            coord_flip() +
            labs( x = NULL, y = "Cell-type proportions", fill = "Sex") +
            theme(plot.title = element_text(size=12, margin=margin(0,0,0,0)))

print(Plot1)

ggsave("/Cell_Composition_Ctrl.png", plot = Plot1, device = "png", width = 92, height = 140, units = 'mm', dpi =300)


###############
#### Braak Lewy body stage 5 group
cell_5_1 <- cell[,-c(1,8:12)]
cell_5_1 <- cell_1[(cell_5_1$Braak_Stage == "5"),]
cell_5_1$Braak_Stage <- NULL

save(cell_5_1, file = "/cells_5.Rda")

# Pearson correlation
Pearson_matrix_5 <-  cor_mat(cell_5_1, method = "pearson")

# get the p.values
Pearson_matrix_p_5 <- Pearson_matrix_5 %>% cor_get_pval()

Pearson_matrix_5 <- as.data.frame(Pearson_matrix_5)
rownames(Pearson_matrix_5) <- Pearson_matrix_5$rowname
Pearson_matrix_5$rowname <- NULL
Pearson_matrix_5 <- mutate_all(Pearson_matrix_5, function(x) as.numeric(as.character(x)))
Pearson_matrix <- as.matrix(Pearson_matrix_5)

Pearson_matrix_p_5 <- as.data.frame(Pearson_matrix_p_5)
rownames(Pearson_matrix_p_5) <- Pearson_matrix_p_5$rowname
Pearson_matrix_p_5$rowname <- NULL
Pearson_matrix_p_5 <- mutate_all(Pearson_matrix_p_5, function(x) as.numeric(as.character(x)))
Pearson_matrix_p_5 <- as.matrix(Pearson_matrix_p_5)


### Plot tau coefficients
# Initialize file path
png("/Cell_pearson_5.png",width=4.5,height=6,units="in",res=600)

corrplot(Pearson_matrix_5, p.mat = Pearson_matrix_p_5, sig.level = 0.05, type = "lower", diag = FALSE, 
         tl.col = "black", tl.srt = 45, insig='blank', tl.cex = 0.7, mar = c(0.5,0.5,0.5,0), cl.cex = 0.7)
Plot <- recordPlot()
print(Plot)
dev.off()


#### Prepare dataset for ancova

cell_2_5 <- cell_5_1
cell_2_5$Sex <- as.factor(cell_2_5$Sex)

######## ANCOVA #############
res.Neurons_5 <- cell_2_5 %>% anova_test(Neurons ~ Age + PMD + RIN + Sex)
get_anova_table(res.Neurons_5)

res.Astrocytes_5 <- cell_2_5 %>% anova_test(Astrocytes ~ Age + PMD + RIN + Sex)
get_anova_table(res.Astrocytes_5)

res.VLMC_5 <- cell_2_5 %>% anova_test(VLMC ~ Age + PMD + RIN + Sex)
get_anova_table(res.VLMC_5)

res.Tanycytes_5 <- cell_2_5 %>% anova_test(Tanycytes ~ Age + PMD + RIN + Sex)
get_anova_table(res.Tanycytes_5)

res.OPC_5 <- cell_2_5 %>% anova_test(OPC ~ Age + PMD + RIN + Sex)
get_anova_table(res.OPC_5)

res.Oligodendrocytes_5 <- cell_2_5 %>% anova_test(Oligodendrocytes ~ Age + PMD + RIN + Sex)
get_anova_table(res.Oligodendrocytes_5)

res.Ependymal_5 <- cell_2_5 %>% anova_test(Ependymal ~ Age + PMD + RIN + Sex)
get_anova_table(res.Ependymal_5)

res.Unknown_5 <- cell_2_5 %>% anova_test(Unknown ~ Age + PMD + RIN + Sex)
get_anova_table(res.Unknown_5)

res.Endothelial_5 <- cell_2_5 %>% anova_test(Endothelial ~ Age + PMD + RIN + Sex)
get_anova_table(res.Endothelial_5)

res.Microglia_5 <- cell_2_5 %>% anova_test(Microglia ~ Age + PMD + RIN + Sex)
get_anova_table(res.Microglia_5)


p_5 <- c(0.131, 0.333, 0.115, 0.149, 0.041, 0.041, 0.138, 0.084, 0.248,0.279)
padj_5 <- p.adjust(p_5, n= length(p_5))



### Plot cell composition (non-adjusted)

# reshape the dataset
cell_3_5 <- cell_2_5
str(cell_3_5)
cell_3_5$Diagnosis <- NULL
cell_3_5$Age <- NULL
cell_3_5$PMD <- NULL
cell_3_5$RIN <- NULL


cell_plot_5 <- aggregate(. ~ Sex, cell_3_5, function(x) list(mean = round(mean(x),3), sd =round(sd(x),3)))
cell_plot_5 <- do.call("data.frame", cell_plot_5) # flatten

Neurons_5 <- cell_plot_5[,1:3]
Neurons_5$cell_type <- c(rep("Neurons", 2))
names(Neurons_5) <- c("Sex", "Mean", "SD", "Cell_Type")
Astrocytes_5 <- cell_plot_5[,c(1,4,5)]
Astrocytes_5$cell_type <- c(rep("Astrocytes", 2))
names(Astrocytes_5) <- c("Sex", "Mean", "SD", "Cell_Type")
VLMC_5<- cell_plot_5[,c(1,6,7)]
VLMC_5$cell_type <- c(rep("VLMC", 2))
names(VLMC_5) <- c("Sex", "Mean", "SD", "Cell_Type")
Tanycytes_5 <- cell_plot_5[,c(1,8,9)]
Tanycytes_5$cell_type <- c(rep("Tanycytes", 2))
names(Tanycytes_5) <- c("Sex", "Mean", "SD", "Cell_Type")
OPC_5 <- cell_plot_5[,c(1,10,11)]
OPC_5$cell_type <- c(rep("OPC", 2))
names(OPC_5) <- c("Sex", "Mean", "SD", "Cell_Type")
Oligodendrocytes_5 <- cell_plot_5[,c(1,12,13)]
Oligodendrocytes_5$cell_type <- c(rep("Oligodendrocytes", 2))
names(Oligodendrocytes_5) <- c("Sex", "Mean", "SD", "Cell_Type")
Ependymal_5 <- cell_plot_5[,c(1,14,15)]
Ependymal_5$cell_type <- c(rep("Ependymal", 2))
names(Ependymal_5) <- c("Sex", "Mean", "SD", "Cell_Type")
Unknown_5 <- cell_plot_5[,c(1,16,17)]
Unknown_5$cell_type <- c(rep("Unknown", 2))
names(Unknown_5) <- c("Sex", "Mean", "SD", "Cell_Type")
Endothelial_5 <- cell_plot_5[,c(1,18,19)]
Endothelial_5$cell_type <- c(rep("Endothelial", 2))
names(Endothelial_5) <- c("Sex", "Mean", "SD", "Cell_Type")
Microglia_5 <- cell_plot_5[,c(1,20,21)]
Microglia_5$cell_type <- c(rep("Microglia", 2))
names(Microglia_5) <- c("Sex", "Mean", "SD", "Cell_Type")
NFO_5 <- cell_plot_5[,c(1,20,21)]
NFO_5$cell_type <- c(rep("NFO", 2))
names(NFO_5) <- c("Sex", "Mean", "SD", "Cell_Type")


cellplot_5 <- rbind(Neurons_5,Astrocytes_5,VLMC_5,Tanycytes_5,OPC_5,Oligodendrocytes_5,Ependymal_5,Unknown_5,Endothelial_5,Microglia_5,NFO_5)

cellplot_5 <- cellplot_5 %>% mutate(Sex =
                                  case_when(Sex == "1" ~ "F", 
                                            Sex == "2" ~ "M" )
)

cellplot_5$Mean <- as.numeric(cellplot_5$Mean)
cellplot_5$SD <- as.numeric(cellplot_5$SD)


my_colors_ID_5 <- c("F" = "#202828",  "M" = "#6a8a9c")

# make stacked bar chart
Plot1_5 <- ggplot(cellplot_5, aes(x= reorder(Cell_Type, Mean),y=Mean,fill=factor(Sex)))+
  geom_bar(stat="identity",position="dodge",color = "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID_5) + 
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,position=position_dodge(.9)) +
  coord_flip() +
  labs( x = NULL, y = "Cell-type proportions", fill = "Sex") +
  theme(plot.title = element_text(size=12, margin=margin(0,0,0,0)))

print(Plot1_5)

ggsave("/Cell_Composition_5.png", plot = Plot1, device = "png", width = 92, height = 140, units = 'mm', dpi =300)

