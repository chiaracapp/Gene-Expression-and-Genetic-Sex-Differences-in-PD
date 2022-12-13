# Gene Expression and Genetic Sex Differences in PD

## Background

This repository contains analysis code for the publication found [here](add link). In this pubblication we performed RNA-sequencing after ribosomal RNA removal of 84 frontal cortex samples. Samples were collected postmortem from 23 non-neurological individuals and 61 individuals with varying degree of α-synuclein-related pathology. The latter were samples collected from donors with iLDB, PD or PDD. Cell composition was assessed for all the samples. Differential expression analysis was performed between males and females both in the non-neurological individual and in patients at Braak Lewy body stage 5. 83 samples weere used for sex-specific expression quantitative trait locus (eQTL) analysis and for sex-specific age-related gene expression analysis.

## Overview

R_Scripts folder contains:
- 1.Transcript_to_Gene_&_Pre-filtering.R
- 2.Cell_Composition.R
- 3.Covariate_selection.R
- 4.DE_analysis.R
- 5.Sex-specific_eQTL &_Age_related_expression_PD.R
- 6.Sex-specific_eQTL &_Age_related_expression_11newloci.R


Files folder contains:
- Coldata_Final_84.csv
- Conjsnps.txt
- Genotypes_PDsnps.txt
- Sample_ID_New.csv
- degradation.matrix.84.txt
- MReads_84.csv
- Expressed_genes_Conj.csv
- Gene_ID.csv
- Nalls_2019_RiskSNP.csv
- scaden_predictions.txt

## Citation

If you use any of the code or data from this repository, please cite our paper.
