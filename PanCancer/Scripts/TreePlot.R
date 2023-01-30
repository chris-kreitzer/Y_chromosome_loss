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
  geom_treemap_subgroup_text(color = 'black')
  

scale_fill_manual(values = c("Peritoneum" = '#7aa8d9',
                               "Bladder/Urinary Tract" = '#f17d81',
                               "Pleura" = '#c3b2d6',
                               "Esophagus/Stomach" = '#7aa8d9',
                               "Testis" = '#f17d81',
                               'Lung' = '#c3b2d6',
                               "CNS/Brain" = 'grey35',
                               "Biliary Tract" = '#4474ba',
                               'Prostate' = '#f17d81',
                               'Soft Tissue' = '#00a89b',
                               "Liver" = '#4474ba',
                               "Head and Neck" = '#049347',
                               "Thyroid" = '#e4e07a',
                               'Bowel' = '#7aa8d9',
                               "Bone" = 'orange',
                               "Skin" = '#bad540', 
                               "Kidney" = '#f17d81',
                               "Breast" = '#e6308e',
                               "Pancreas" = '#4474ba',
                               "Other" = 'white',
                               "Peripheral Nervous System" = 'grey35',
                               "Thymus" = '#e4e07a',
                               "Adrenal Gland" = '#e4e07a',
                               "Eye" = '#049347',
                               "Penis" = '#f17d81', 
                               "Ampulla of Vater" = '#4474ba',
                               "Vulva/Vagina" = 'white',
                               "Ovary/Fallopian Tube" = '#f8b74a'))

#' out