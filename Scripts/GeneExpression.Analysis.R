## GTEx Portal; Genotype-Tissue Expression project:
## tissue-specific gene expression and regulation: samples were collected
## from 54 non-diseased tissue sites across nearly 1000 individuals
##
## I will use these data, to highlight which Y-chromosome genes are ubiquitously expressed
## and which are testis specfic
## 
## Start analysis: 10/23/2020
## Start revision: 03/25/2021
## Start revision: 02/16/2022
## Start revision: 06/27/2022


clean()
setup(working.path = '~/Documents/GitHub/Y_chromosome_loss/')
dataStorage = '~/Documents/MSKCC/00_Data/'

## Libraries:
library(vroom)
library(ComplexHeatmap)
library(grid)
library(grDevices)
library(circlize)


## Input:
#' TPM values for male specific tissue on Y-chromosome
#' Median gene-level TPM by tissue
GTEX.TPM.raw_df = vroom::vroom(file = paste0(dataStorage, 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'), 
                               delim = '\t', 
                               skip = 2)
GTEX.Sex.Annotation_df = vroom::vroom(file = paste0(dataStorage, 'GTEx_SexAnnotation.txt'), 
                                      delim = '\t')
GTEX.annotation_df = vroom::vroom(file = paste0(dataStorage, 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'), 
                                  delim = '\t')
NCBI = read.csv(file = 'NCBI_Y.chromosome.txt', sep = '\t')
NCBI = NCBI[!duplicated(NCBI$GeneID), ]
write.table('~/Documents/MSKCC/10_MasterThesis/Data/NCBI_GeneTableY.txt', sep = '\t', row.names = F, quote = F)


#' Gene data from BIOMART
GenesY_df = read.csv(file = paste0(dataStorage, 'Gene.Model.Biomart.txt'), 
                                   sep = '\t')
GenesY_df = GenesY_df[which(GenesY_df$Chromosome.scaffold.name == 'Y'), ]
GenesY_df = GenesY_df[!duplicated(GenesY_df$Gene.stable.ID.version), ]



## Data Processing:
male_samples = GTEX.Sex.Annotation_df[which(GTEX.Sex.Annotation_df$SEX == 1), 'SUBJID']
male_samples = as.character(male_samples$SUBJID)

#' fetch male samples from expression table
index_colnames = c()
for(i in 1:length(male_samples)){
  ii = grep(pattern = male_samples[i], x = colnames(GTEX.TPM.raw_df))
  index_colnames = c(index_colnames, ii)
}

male_expression = GTEX.TPM.raw_df[, c(1, 2, index_colnames)]


#' replace sample IDs with tissue type from annotation file
#' column SMTS == Tissue Type, area from which the tissue sample was taken. Parent value to SMTSD
#' column SMTSD == Tissue Type, more specific detail of tissue type
#' in this case I will use the SMTS column (parent); ending up in 29 tissue sites
for(i in 3:ncol(male_expression)){
  if(colnames(male_expression)[i] %in% GTEX.annotation_df$SAMPID){
    colnames(male_expression)[i] = GTEX.annotation_df$SMTS[which(colnames(male_expression)[i] == GTEX.annotation_df$SAMPID)]
  }
}


#' Y-chromosome genes to keep
GenesY_nonCoding = GenesY_df[which(GenesY_df$Gene.type == 'processed_pseudogene' | 
                            GenesY_df$Gene.type == 'lncRNA' |
                            GenesY_df$Gene.type == 'snRNA' | 
                            GenesY_df$Gene.type == 'snoRNA'), ]
GenesY_coding = GenesY_df[which(GenesY_df$Gene.type == 'protein_coding'), ]

GenesY_df = GenesY_df[!duplicated(GenesY_df$Gene.stable.ID.version), ]
GenesY_df = GenesY_df[!duplicated(GenesY_df$Gene.name),, drop = F]
GenesY_df = GenesY_df[, c('Gene.name', 'Gene.start..bp.', 'Gene.end..bp.')]
colnames(GenesY_df) = c('Gene', 'Start', 'End')

GTExY.expression = male_expression[which(male_expression$Description %in% GenesY_df$Gene),, drop = F]
GTExY.expression = GTExY.expression[!duplicated(GTExY.expression$Description),, drop = F]



###############################################################################
###############################################################################
#' subset gene table to those which are protein coding genes from the 
#' Y-chromosome
Y_protein_coding = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/proteinCodingY.txt', sep = '\t', skip = 2)
Y_protein_coding = Y_protein_coding$Gene.name

GTEX_proteinCoding = GTExY.expression[which(GTExY.expression$Description %in% Y_protein_coding), ]


#' Analysis:
#' calculate median expression of genes
GTEx.analysis = GTEX_proteinCoding[, 3:ncol(GTEX_proteinCoding)]
colnames(GTEx.analysis) = gsub('\\..*', x = colnames(GTEx.analysis), replacement = '')
row.names(GTEx.analysis) = GTEX_proteinCoding$Description

# calculate median expression across different tissues
GTEX_expression.median.summary = t(apply(GTEx.analysis, 1, function(x) tapply(x, colnames(GTEx.analysis), median)))
GTEX_expression.median.summary = as.data.frame(GTEX_expression.median.summary)

#' Log transformation (for visualization)
GTEX_expression.median.summary = log2(GTEX_expression.median.summary[] + 1)

#' #' remove genes with no expression 
#' GTEX_expression.median.summary = GTEX_expression.median.summary[rowSums(GTEX_expression.median.summary) > 0, ]
#' #' remove genes with median expression == 0
#' gene_median = apply(GTEX_expression.median.summary, 1, median)
#' gene_keep = names(which(apply(GTEX_expression.median.summary, 1, median) != 0))
#' 
#' #' subset parent data-frame
#' GTEX_expression.median.summary = GTEX_expression.median.summary[which(row.names(GTEX_expression.median.summary) %in% gene_keep),, drop = F]


##############################
## Visualization:
#' ComplexHeatmap:
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

#' row split in X-transposed, X-degenerate and Ampliconic
x_transposed = c('TGIF2LY',
                 'PCDH11Y')
x_degenerate = c('SRY',
                 'RPS4Y1',
                 'ZFY',
                 'AMELY',
                 'TBL1Y',
                 'USP9Y',
                 'DDX3Y',
                 'UTY',
                 'TMSB4Y',
                 'NLGN4Y',
                 'KDM5D',
                 'EIF1AY',
                 'RPS4Y2')
ampliconic = setdiff(Y_protein_coding, union(x_transposed, x_degenerate))

meta.row = data.frame(gene = unique(GTEX_proteinCoding$Description))
meta.row$region = NA
for(i in 1:nrow(meta.row)){
  if(meta.row$gene[i] %in% x_transposed){
    meta.row$region[i] = 'x_transposed'
  } else if(meta.row$gene[i] %in% x_degenerate){
    meta.row$region[i] = 'x_degenerate'
  } else {
    meta.row$region[i] = 'ampliconic'
  }
}
row.names(meta.row) = meta.row$gene
meta.row$gene = NULL


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
                            row_names_gp = gpar(fontsize = 10),
                            cluster_columns = FALSE,
                            column_split = meta.information$grouping,
                            row_split = meta.row$region,
                            row_dend_reorder = F,
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


## Look at the expression of X-homologues
X_homologues = c('SOX3',
                 'RPS4X',
                 'ZFX',
                 'AMELX',
                 'TBL1X',
                 'USP9X',
                   'KDM6A',
                 'TMSB4X',
                 'NLGN4X',
                   'EIF1AX',
                 'RPS4X')

GTEX_X_homologue = male_expression[which(male_expression$Description %in% X_homologues), ]
GTEX_X = GTEX_X_homologue[, 3:ncol(GTEX_X_homologue)]
colnames(GTEX_X) = gsub('\\..*', x = colnames(GTEX_X), replacement = '')
row.names(GTEX_X) = GTEX_X_homologue$Description

# calculate median expression across different tissues
GTEX_X_median = t(apply(GTEX_X, 1, function(x) tapply(x, colnames(GTEX_X), median)))
GTEX_X.median.summary = as.data.frame(GTEX_X_median)

#' Log transformation (for visualization)
GTEX_X.median.summary = log2(GTEX_X.median.summary[] + 1)


#' Visualization
x_chromosome = ComplexHeatmap::Heatmap(matrix = GTEX_X.median.summary, 
                                       rect_gp = gpar(col = "white", lwd = 2),
                                       col = colorRamp2(c(min(GTEX_X.median.summary), 
                                                          max(GTEX_X.median.summary)), 
                                                        c("#DFF4F7", "#006395")),
                                       column_names_rot = 45,
                                       column_names_gp = gpar(fontsize = 8),
                                       row_names_gp = gpar(fontsize = 10),
                                       cluster_columns = FALSE,
                                       column_split = meta.information$grouping,
                                       #row_split = meta.row$region,
                                       #row_dend_reorder = F,
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
x_chromosome


###### X-Y expression differences
install.packages('plotrix')
library(plotrix)

Y_protein_coding = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/proteinCodingY.txt', sep = '\t', skip = 2)
XY_homologue = Y_protein_coding[which(Y_protein_coding$X.1 != ''), ]
XY_homologue = XY_homologue[which(XY_homologue$X.1 != '-'), ]
XY_homologue = XY_homologue[-13, ]
XY_homologue[9, 'X.1'] = 'KDM6A'

#' just look into prostate tissue
Prostate_tissue = male_expression[,c(1,2, which(colnames(male_expression) == 'Prostate'))]

all_out = data.frame()
for(i in 1:nrow(XY_homologue)){
  Y_expression = as.numeric(as.data.frame(Prostate_tissue[which(Prostate_tissue$Description == XY_homologue$Gene.name[i]), c(3:ncol(Prostate_tissue))]))
  Y_expression_median = median(as.numeric(Y_expression))
  X_expression = as.numeric(as.data.frame(Prostate_tissue[which(Prostate_tissue$Description == XY_homologue$X.1[i]), c(3:ncol(Prostate_tissue))])) 
  X_expression_median = median(as.numeric(X_expression))
  test = t.test(Y_expression, X_expression, paired = T)$p.value
  
  out = data.frame(Y_gene = XY_homologue$Gene.name[i],
                   X_gene = XY_homologue$X.1[i],
                   Y_median = Y_expression_median,
                   X_median = X_expression_median,
                   statistic = test)

  all_out = rbind(all_out, out)
}


a = all_out[, c(1, 3, 5)]
colnames(a) = c('name', 'value', 'statistic')
b = all_out[, c(2, 4, 5)]
colnames(b) = c('name', 'value', 'statistic')

plot_out = rbind(a,b)
plot_out$group = c(rep('Y-linked', 13), rep('X-linked', 13))
plot_out$subject = c(seq(1,13, 1), seq(1, 13, 1))


p = ggplot(plot_out, aes(x = group, y = value, group = subject)) +
  geom_line(size = 0.4) +
  geom_point(aes(col = group)) +
  geom_text(aes(label = name), size = 2, vjust = 0, hjust = 0) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 1)
p + annotation_logticks(sides = 'l')




library(scales)
library(MASS)
