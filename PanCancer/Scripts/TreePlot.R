##----------------+
## CancerTypes; cohort
## overview;
##
## start: 12/13/2022
## revision: 01/30/2023
## revision: 06/07/2023
## 
## chris-kreitzer
##----------------+

clean()
gc()
setwd('~/Documents/MSKCC/10_MasterThesis/')
#install.packages("treemapify")
library(treemapify)
library(webr)
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')

## TODO:
#' map the final samples; rather than full cohort
#' go with GI tract (Bastien); 
#' rather than cancer type detailed


impact = readRDS('Data/00_CohortData/Cohort_071322.rds')
impact = impact[which(impact$Study_include == 'yes'), ]
# oncoMuts = read.csv('Data/signedOut/data_mutations_extended.oncokb.txt.gz', sep = '\t')
impact_group = impact[,c('CANCER_TYPE', 'CANCER_TYPE_ritika', 'CANCER_TYPE_DETAILED_ritika', 'ONCOTREE_CODE', 'Age_Sequencing', 'SAMPLE_TYPE')]
load('Data/00_CohortData/0.0.esthetics.Rdata')
esthetics$names



##------- Assign hypermutator phenotype
##  hypermutated (HM) tumors (COAD):
##  - oncogenic POLE mutation
##  - MSIsensor > 10
##  - (Niu et al., 2014)
# oncoPOLE = oncoMuts[which(oncoMuts$Hugo_Symbol == 'POLE'), ]
# oncoPOLE = oncoPOLE[which(oncoPOLE$ONCOGENIC %in% c('Likely Oncogenic', 'Oncogenic')), ]
# impact = merge(impact, oncoPOLE[, c('Tumor_Sample_Barcode', 'Hugo_Symbol')], by.x = 'SAMPLE_ID', by.y = 'Tumor_Sample_Barcode', all.x = T)
# impact$hypermutated = ifelse(!is.na(impact$Hugo_Symbol) | impact$MSI_SCORE >= 10, 'yes', 'no')
# impact$CANCER_TYPE[which(impact$CANCER_TYPE == 'Colorectal Cancer' & impact$hypermutated == 'yes')] = paste0('Colorectal Hypermutated')
# impact$CANCER_TYPE[which(impact$CANCER_TYPE == 'Colorectal Cancer' & impact$hypermutated == 'no')] = paste0('Colorectal MSS')
# 
# colnames(impact)[40] = 'POLE_mut'
# impact$POLE_mut = ifelse(impact$POLE_mut == 'POLE', 'yes', 'no')
# 
# cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
# cohort = merge(cohort, impact[,c('SAMPLE_ID', 'POLE_mut', 'hypermutated')], by.x = 'SAMPLE_ID', by.y = 'SAMPLE_ID', all.x = T)
# 
# cohort$CANCER_TYPE[which(cohort$CANCER_TYPE == 'Colorectal Cancer' & cohort$hypermutated == 'yes')] = paste0('Colorectal Hypermutated')
# cohort$CANCER_TYPE[which(cohort$CANCER_TYPE == 'Colorectal Cancer' & cohort$hypermutated == 'no')]  = paste0('Colorectal MSS')
# cohort$CANCER_TYPE[which(cohort$CANCER_TYPE == 'Colorectal Cancer' & is.na(cohort$hypermutated))] = paste0('Colorectal MSS')
# saveRDS(object = cohort, file = 'Data/00_CohortData/Cohort_071322.rds')



##----------------+
## Cancer type summary stats
##----------------+
all_types = data.frame()
for(i in unique(impact_group$CANCER_TYPE)){
  main_name = i
  sub_data = impact_group[which(impact_group$CANCER_TYPE == i), ]
  for(j in unique(sub_data$ONCOTREE_CODE)){
    sub_name = j
    cancer_type_detailed = unique(sub_data$CANCER_TYPE_DETAILED_ritika[which(sub_data$ONCOTREE_CODE == j)])
    sub_length = length(sub_data$ONCOTREE_CODE[which(sub_data$ONCOTREE_CODE == j)])
    mean_age = mean(sub_data$Age_Sequencing[which(sub_data$ONCOTREE_CODE == j)], na.rm = T)
    if(any('Primary' %in% sub_data$SAMPLE_TYPE[which(sub_data$ONCOTREE_CODE == j)])){
      mean_primary = table(sub_data$SAMPLE_TYPE[which(sub_data$ONCOTREE_CODE == j)])[['Primary']] / sum(table(sub_data$SAMPLE_TYPE[which(sub_data$ONCOTREE_CODE == j)]))
      mean_primary = mean_primary * 100
      mean_mets = 100 - mean_primary
    } else {
      mean_primary = NA
      mean_mets = NA
    }
    
    out = data.frame(CancerType = i,
                     CancerTypeDetailed = cancer_type_detailed,
                     ONCOTREE = sub_name,
                     n = sub_length,
                     mean_Age = mean_age,
                     primary = mean_primary,
                     metastasis = mean_mets)
    
    all_types = rbind(all_types, out)
  }
}

write.table(all_types, file = 'Data/00_CohortData/CancerTypes_Summary.txt', sep = '\t', row.names = F, quote = F)


##-------
## Treeplot visualization
##-------
TreePlot = ggplot(all_types, aes(area = n, fill = CancerType, 
                      label = CancerType, 
                      subgroup = ONCOTREE)) +
  geom_treemap() +
  geom_treemap_text(color = 'black') +
  geom_treemap_subgroup_text(color = 'black') +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c('Mesothelioma' = '#682767',
                             'Bladder Cancer' = '#f17d81',
                             'Non-Small Cell Lung Cancer' = '#B368D9',
                             'Esophagogastric Cancer' = '#7aa8d9',
                             'Germ Cell Tumor' = '#7e1918',
                             'Prostate Cancer' = '#be1e2d',
                             "Hepatobiliary Cancer" = '#4474ba',
                             'Thyroid Cancer' = '#cccc33',
                             'Anal Cancer' = '#2e3a46',
                             'Bone Cancer' = '#00a89b',
                             'Soft Tissue Sarcoma' = '#005c54',
                             'Gastrointestinal Stromal Tumor' = '#3d546c',
                             'Salivary Gland Cancer' = '#026631',
                             'Pancreatic Cancer' = '#355a91',
                             'Skin Cancer, Non-Melanoma' = '#9acd32',
                             'Glioma' = '#bf9fc5',
                             'Head and Neck Cancer' = '#049347',
                             'Colorectal MSS' = '#007eb5',
                             'Colorectal Hypermutated' = '#9eddf9',
                             'Melanoma' = '#4ebe03',
                             'Renal Cell Carcinoma' = '#f8afb3',
                             'Small Bowel Cancer' = '#006490',
                             'Penile Cancer' = '#7e1918',
                             'Breast Cancer' = '#e6308e',
                             'Embryonal Tumor' = '#55934d',
                             'Thymic Tumor' = '#ba9495',
                             'Nerve Sheath Tumor' = 'grey85',
                             'Sellar Tumor' = 'grey85',
                             'Appendiceal Cancer' = '#2e3a46',
                             'Miscellaneous Neuroepithelial Tumor' = 'grey85',
                             'Gastrointestinal Neuroendocrine Tumor' = '#557597',
                             'Miscellaneous Brain Tumor' = 'grey85',
                             'Ampullary Cancer' = '#c6d5ea',
                             'Sex Cord Stromal Tumor' = '#7e1918',
                             'Peripheral Nervous System' = 'grey85',
                             'Pheochromocytoma' = '#026631',
                             'Wilms Tumor' = '#f8afb3',
                             'Choroid Plexus Tumor' = 'grey85',
                             'Primary CNS Melanocytic Tumors' = 'grey85',
                             'Tubular Adenoma of the Colon' = '#9eddf9',
                             'Gastrointestinal Neuroendocrine Tumors of the Esophagus/Stomach' = '#557597',
                             'Parathyroid Cancer' = '#cccc33',
                             'Melanocytoma' = 'grey85', 
                             'Malignant Glomus Tumor' = '#00dbca'))
                             
ggsave_golden(filename = 'Figures_original/TreePlot.pdf', plot = TreePlot, width = 9)




##----------------+
## Cohort overview: 
## PieDonut
##----------------+
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
cancer = sort(table(cohort$CANCER_TYPE), decreasing = T)
cancer = cancer[1:19]

bigcohort = cohort[which(cohort$CANCER_TYPE %in% names(cancer)), c('CANCER_TYPE', 'CANCER_TYPE_DETAILED_ritika', 'ONCOTREE_CODE')]
cc = data.frame()
for(i in unique(bigcohort$CANCER_TYPE)){
  data_sub = sort(table(bigcohort$CANCER_TYPE_DETAILED_ritika[which(bigcohort$CANCER_TYPE == i)]), decreasing = T)
  data_sub = data_sub[1:2]
  out = data.frame(cancer = i,
                   cancer_detailed = names(data_sub))
  cc = rbind(cc, out)
}

pie_cohort = cohort
pie_cohort$CANCER_TYPE[!pie_cohort$CANCER_TYPE %in% cc$cancer] = 'Other'
uu = data.frame()
for(i in unique(pie_cohort$CANCER_TYPE)){
  data_sub = pie_cohort[which(pie_cohort$CANCER_TYPE == i), ]
  data_sub$CANCER_TYPE_DETAILED_ritika[!data_sub$CANCER_TYPE_DETAILED_ritika %in% cc$cancer_detailed[which(cc$cancer == i)]] = 'Other'
  uu = rbind(uu, data_sub)
}

cohortPie = PieDonut(uu, 
         aes(CANCER_TYPE, CANCER_TYPE_DETAILED_ritika),
         labelposition = 1, labelpositionThreshold = 0.9,
         showPieName = F,
         pieLabelSize = 4,
         donutLabelSize = 4,
         r0=0.0, r1 = 0.9, start=pi/2,explodeDonut=TRUE)

ggsave_golden(filename = 'Figures_original/PieDonut_Cohort.pdf', plot = cohortPie, width = 12)


##-- Numbers
a_out = data.frame()
for(i in unique(uu$CANCER_TYPE)){
  datasub = uu[which(uu$CANCER_TYPE == i), ]
  n = nrow(datasub)
  for(j in unique(datasub$CANCER_TYPE_DETAILED_ritika)){
    datasub2 = datasub[which(datasub$CANCER_TYPE_DETAILED_ritika == j), ]
    n2 = nrow(datasub2)
    out = data.frame(cancer_type = i,
                     cancer_type_detailed = j,
                     n = n2,
                     Frequency = round((n2 / n) * 100, 2))
    a_out = rbind(a_out, out)
  }
}


write.xlsx(x = a_out, file = 'Data/00_CohortData/CancerTypes.xlsx')
write.xlsx(x = aa, file = 'Data/00_CohortData/CancerTypes_Summary.xlsx')

#' out