##----------------+
## CancerTypes; cohort
## overview;
##
## start: 12/13/2022
## revision: 01/30/2023
## 
## chris-kreitzer
##----------------+

clean()
gc()
setwd('~/Documents/MSKCC/10_MasterThesis/')
#install.packages("treemapify")
library(treemapify)

## TODO:
#' map the final samples; rather than full cohort
#' go with GI tract (Bastien); 
#' rather than cancer type detailed



impact = readRDS('Data/00_CohortData/Cohort_071322.rds')
impact = impact[which(impact$Study_include == 'yes'), ]
impact_group = impact[,c('CANCER_TYPE', 'CANCER_TYPE_ritika', 'CANCER_TYPE_DETAILED_ritika', 'ONCOTREE_CODE', 'Age_Sequencing')]
load('Data/00_CohortData/0.0.esthetics.Rdata')
esthetics$names


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
    out = data.frame(CancerType = i,
                     CancerTypeDetailed = cancer_type_detailed,
                     ONCOTREE = sub_name,
                     n = sub_length,
                     mean_Age = mean_age)
    
    all_types = rbind(all_types, out)
  }
}

all_types$MAIN = sub("\\s+[^ ]+$", "", all_types$meta_main)
ggplot(all_types, aes(area = n, fill = CancerType, 
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
                             'Colorectal Cancer' = '#9eddf9',
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
                             


#' out