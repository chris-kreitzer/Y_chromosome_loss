## GTEx Portal; Genotype-Tissue Expression project:
## tissue-specific gene expression and regulation: samples were collected
## from 54 non-diseased tissue sites across nearly 1000 individuals
##
## I will use these data, to highlight which Y-chromosome genes are ubiquitinously expressed in the
## body - in non-diseased tissue
## 
## Start analysis: 10/23/2020
## Start revision: 03/25/2021
## 


rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/Y_chromosome_loss/')


## Libraries:
library(vroom)
library(ComplexHeatmap)
library(grid)
library(grDevices)
library(circlize)


## Input:
#' TPM values for male specific tissue on Y-chromosome
GTEX.TPM.raw_df = vroom::vroom(file = '~/Documents/MSKCC/00_Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', delim = '\t', skip = 2)
GTEX.Sex.Annotation_df = vroom::vroom(file = '~/Documents/MSKCC/00_Data/GTEx_SexAnnotation.txt', delim = '\t')
GTEX.annotation_df = vroom::vroom(file = '~/Documents/MSKCC/00_Data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', delim = '\t')

#' Gene data from BIOMART
GenesY_df = read.csv(file = '~/Documents/MSKCC/00_Data/Gene.Model.Biomart.txt', sep = '\t')

## Functions:


## Processing:
#' selcet only male samples: (1 = male)
GTEX.TPM.annotation_df = GTEX.TPM.raw_df[, c(1, 2)]
male_samples = GTEX.Sex.Annotation_df[which(GTEX.Sex.Annotation_df$SEX == 1), 'SUBJID']
male_samples = as.character(male_samples$SUBJID)

#' grep male columns
index_colnames = c()
for(i in 1:length(male_samples)){
  ii = grep(pattern = male_samples[i], x = colnames(GTEX.TPM.raw_df))
  index_colnames = c(index_colnames, ii)
}

GTEX.TPM.raw_df = GTEX.TPM.raw_df[, index_colnames]
GTEX.TPM.raw_df = cbind(GTEX.TPM.annotation_df, GTEX.TPM.raw_df)

#' replace SampleIDs with tissue type from annotation file
#' column SMTS == Tissue Type, area from which the tissue sample was taken. Parent value to SMTSD
#' column SMTSD == Tissue Type, more specific detail of tissue type

#' in this case I will use the SMTS column (parent); ending up in 29 tissue sites
for(i in 3:ncol(GTEX.TPM.raw_df)){
  if(colnames(GTEX.TPM.raw_df)[i] %in% GTEX.annotation_df$SAMPID){
    colnames(GTEX.TPM.raw_df)[i] = GTEX.annotation_df$SMTS[which(colnames(GTEX.TPM.raw_df)[i] == GTEX.annotation_df$SAMPID)]
  }
}


#' Y-chromosome genes to keep
#' Target gene list from Hanisha | Subhi (299) (under input section):
GenesY_df = GenesY_df[GenesY_df$Chromosome.scaffold.name == 'Y', ]
GenesY_df = GenesY_df[which(GenesY_df$Gene.type == 'protein_coding' | 
                              GenesY_df$Gene.type == 'processed_pseudogene' | 
                              GenesY_df$Gene.type == 'lncRNA' |
                              GenesY_df$Gene.type == 'snRNA' | 
                              GenesY_df$Gene.type == 'snoRNA'), ]
GenesY_df = GenesY_df[!duplicated(GenesY_df$Gene.stable.ID.version), ]
GenesY_df = GenesY_df[!duplicated(GenesY_df$Gene.name),, drop = F]
GenesY_df = GenesY_df[, c('Gene.name', 'Gene.start..bp.', 'Gene.end..bp.')]
colnames(GenesY_df) = c('Gene', 'Start', 'End')

GTExY.expression = GTEX.TPM.raw_df[which(GTEX.TPM.raw_df$Description %in% GenesY_df$Gene),, drop = F]
GTExY.expression = GTExY.expression[!duplicated(GTExY.expression$Description),, drop = F]

#' remove tissue sites < 10 samples
table(colnames(GTExY.expression)[3:ncol(GTExY.expression)]) # non to remove



## Analysis:
#' calculate median expression of genes
GTEx.analysis = GTExY.expression[, 3:ncol(GTExY.expression)]
colnames(GTEx.analysis) = gsub('\\..*', x = colnames(GTEx.analysis), replacement = '')
row.names(GTEx.analysis) = GTExY.expression$Description

# calculate median expression across different tissues
GTEX_expression.median.summary = t(apply(GTEx.analysis, 1, function(x) tapply(x, colnames(GTEx.analysis), median)))
GTEX_expression.median.summary = as.data.frame(GTEX_expression.median.summary)

#' Log transformation (for visualization)
GTEX_expression.median.summary = log2(GTEX_expression.median.summary[] + 1)

#' remove genes with no expression 
GTEX_expression.median.summary = GTEX_expression.median.summary[rowSums(GTEX_expression.median.summary) > 0, ]
#' remove genes with median expression == 0
gene_median = apply(GTEX_expression.median.summary, 1, median)
gene_keep = names(which(apply(GTEX_expression.median.summary, 1, median) != 0))

#' subset parent data-frame
GTEX_expression.median.summary = GTEX_expression.median.summary[which(row.names(GTEX_expression.median.summary) %in% gene_keep),, drop = F]

## Visualization:
#' ComplexHeatmap: Heatmap
#' add sample number 
Freq.tissue = as.data.frame(table(colnames(GTEx.analysis)))
colnames(GTEX_expression.median.summary) = paste(colnames(GTEX_expression.median.summary), 
                                                 paste0('(n=', Freq.tissue$Freq[which(Freq.tissue$Var1 == colnames(GTEX_expression.median.summary))], ')'))


#' Column split in male and NON-male tissues
meta.information = data.frame(Tissue = colnames(GTEX_expression.median.summary))
meta.information$grouping = rep('non male specific tissue', nrow(meta.information))
meta.information$grouping[grep('^Prostate.*', meta.information$Tissue)] = 'male specific tissue'
meta.information$grouping[grep('^Testis.*', meta.information$Tissue)] = 'male specific tissue'
row.names(meta.information) = meta.information$Tissue
meta.information$Tissue = NULL

#' Heatmap body
y_chromosome = ComplexHeatmap::Heatmap(matrix = GTEX_expression.median.summary, 
                            #row_dend_width = unit(2, "cm"),
                            rect_gp = gpar(col = "white", lwd = 2),
                            col = colorRamp2(c(min(GTEX_expression.median.summary), 
                                               #max(GTEX_expression.median.summary) / 2, 
                                               max(GTEX_expression.median.summary)), 
                                             c("#DFF4F7", "#006395")),
                            column_names_rot = 45,
                            column_names_gp = gpar(fontsize = 8),
                            #row_km = 3,
                            #row_km_repeats = 100,
                            row_names_gp = gpar(fontsize = 10),
                            cluster_columns = FALSE,
                            column_split = meta.information$grouping,
                            row_dend_reorder = TRUE,
                            border = T,
                            column_gap = unit(4, "mm"),
                            heatmap_legend_param = list(title = "Expression\n(TPM)", 
                                                        at = c(0.04702180, 9), 
                                                        labels = c("low", 'high'), 
                                                        border = "black",
                                                        legend_height = unit(4, 'cm'),
                                                        grid_width = unit(3, 'mm'),
                                                        direction = 'horizontal',
                                                        title_positon = 'topcenter',
                                                        annotation_legend_side = 'bottom'),
                            
                            width = unit(15, "cm"), 
                            height = unit(15, "cm"))
y_chromosome

# save as 10 x 10
