##----------------+
## Data selection for 
## Y-chromosome loss determination
##----------------+
##
## Read in clinical data from cBioPortal: 07/11/2022
## This should be at the patient level  
## Need to select one sample per patient 
## 
## start: 07/11/2022
## revision: 08/24/2022
## revision: 12/02/2022
## revision: 12/07/2022
## chris-kreitzer


clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/')

library(dplyr)
library(data.table)
library(usefun)


##----------------+
## load clinical data and 
## merge data_clinical_oncoKB with data_from the portal
##----------------+
portal_data = fread('Data/signedOut/mskimpact_clinical_data_07112022.tsv', stringsAsFactors = F, data.table = F)
colnames(portal_data) = toupper(colnames(portal_data))
colnames(portal_data) = gsub(pattern = ' ', replacement = '_', colnames(portal_data), fixed = T )
rownames(portal_data) = portal_data$SAMPLE_ID
FacetsPaths = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/facets_annotated.cohort.txt.gz', sep = '\t')
FacetsPaths = FacetsPaths[!duplicated(FacetsPaths$tumor_sample), ]

#' exclude heme patients;
#' and any lymphoid; myeloid cancers
portal_data$STUDY_ID = NULL
portal_data$DIAGNOSIS_AGE = NULL
portal_data = portal_data[!portal_data$GENE_PANEL %in% c("ACCESS129", "IMPACT-HEME-400", "IMPACT-HEME-468"), ]
portal_data$`AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)` = as.character(as.numeric(portal_data$`AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)`))
portal_data$TUMOR_PURITY = as.numeric(as.character(portal_data$TUMOR_PURITY))


tropism = readxl::read_excel('~/Documents/MSKCC/10_MasterThesis/Data/00_CohortData/Tropism_Bastien.xlsx', sheet = 'Table S1B', skip = 2)
tropism = tropism[,c('sample_id', 'curated_organ_system', 'cancer_type', 'oncotree_code', 'curated_subtype_display', 'curated_subtype_abbr')]
tropism = tropism[!duplicated(tropism$oncotree_code), ]

portal_data$curated_organ_system = NA
portal_data$curated_subtype_display = NA
for(i in 1:nrow(portal_data)){
  print(i)
  if(portal_data$ONCOTREE_CODE[i] %in% tropism$oncotree_code){
    portal_data$curated_organ_system[i] = tropism$curated_organ_system[which(tropism$oncotree_code == portal_data$ONCOTREE_CODE[i])]
    portal_data$curated_subtype_display[i] = tropism$curated_subtype_display[which(tropism$oncotree_code == portal_data$ONCOTREE_CODE[i])]
  } else {
    portal_data$curated_organ_system[i] = NA
    portal_data$curated_subtype_display[i] = NA
  }
}

portal_data = portal_data[, c('SAMPLE_ID','PATIENT_ID', 'SEX', 'RACE_CATEGORY', 'curated_organ_system',
                              'CANCER_TYPE', 'ONCOTREE_CODE', 'curated_subtype_display', 'CANCER_TYPE_DETAILED',
                              'PRIMARY_TUMOR_SITE', 'SAMPLE_TYPE', 'METASTATIC_SITE', 'GENE_PANEL',
                              'SAMPLE_COVERAGE', 'TUMOR_PURITY', 'MSI_TYPE', 'MSI_SCORE', 
                              'FRACTION_GENOME_ALTERED', 'IMPACT_TMB_SCORE', 'MUTATION_COUNT',
                              'AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)',
                              'OVERALL_SURVIVAL_STATUS', 'OVERALL_SURVIVAL_(MONTHS)')]

##----------------+
## Oncotree Annotations;
##----------------+

# oncotree_table = read.csv('Data/00_CohortData/oncotree.json', sep = '\t')
# oncotree_table$metanci = NULL
# oncotree_table$metaumls = NULL
# oncotree_table$history = NULL
# oncotree_table$level_7 = NULL
# 
# ne = data.frame()
# for(i in 1:nrow(oncotree_table)){
#   print(i)
#   new = data.frame(parent = paste(oncotree_table[i,c(1,2,3,4,5,6)], sep = ';', collapse = ';'))
#   ne = rbind(ne, new)
# }
# ne$parent = gsub(pattern = ';+', ';', ne$parent)
# ne$parent = substr(x = ne$parent, start = 1, stop = nchar(ne$parent) - 1)
# ne$last = NA
# ne$cancertype = NA
# for(i in 1:nrow(ne)){
#   all_types = unlist(strsplit(x = ne$parent[i], split = ';', fixed = T))
#   all_types_length = length(all_types)
#   ne$last[i] = all_types[all_types_length]
#   ne$cancertype[i] = all_types[1]
# }
# 
# ne$code = gsub("[\\(\\)]", "", regmatches(ne$last,gregexpr("\\(.*?\\)", ne$last)))
# ne = ne[,c(3,2,4,1)]
# colnames(ne) = c('CancerType', 'cancer_type_detail', 'oncotree_code', 'name_pileup')
# write.table(ne, file = 'Data/00_CohortData/Oncotree_Code_formated.txt', sep = '\t', row.names = F, quote = F)

oncotree = read.csv('Data/00_CohortData/Oncotree_Code_formated.txt', sep = '\t', header = F)
colnames(oncotree) = c('CancerType', 'cancer_type_detail', 'oncotree_code', 'name_pileup')

portal_data$CANCER_TYPE_ritika = NA
portal_data$CANCER_TYPE_DETAILED_ritika = NA

for(i in 1:nrow(portal_data)){
  print(i)
  if(portal_data$ONCOTREE_CODE[i] %in% oncotree$oncotree_code){
    portal_data$CANCER_TYPE_ritika[i] = oncotree$CancerType[which(oncotree$oncotree_code == portal_data$ONCOTREE_CODE[i])]
    portal_data$CANCER_TYPE_DETAILED_ritika[i] = oncotree$cancer_type_detail[which(oncotree$oncotree_code == portal_data$ONCOTREE_CODE[i])]
  } else {
    portal_data$CANCER_TYPE_ritika[i] = NA
    portal_data$CANCER_TYPE_DETAILED_ritika[i] = NA
  }
}

portal_data = portal_data[!is.na(portal_data$CANCER_TYPE), ]
portal_data = portal_data[!portal_data$CANCER_TYPE_ritika %in% c('Myeloid (MYELOID)', 
                                                                 'Lymphoid (LYMPH)'), ] 


##----------------+
## only male samples;
##----------------+
portal_data = portal_data[which(portal_data$SEX == 'Male'), ]


##----------------+
## one sample per patient;
##----------------+
## if just one sample per patient --> use it (regardless of primary or metastasis)
## if multiple samples, priority is given to i) primaries, ii) tumor purity and iii) coverage

sample_selection = function(data){
  data.analysis = data.frame()
  data.in = as.data.frame(data)
  
  for(i in unique(data.in$PATIENT_ID)){
    try({
      print(i)
      # check how many samples per patient
      data.sub = data.in[which(data.in$PATIENT_ID == i), ]
      
      # one sample; go for it
      if(nrow(data.sub) == 1){ 
        chosen.sample = data.sub$SAMPLE_ID
        sample.type = data.sub$SAMPLE_TYPE
        metastatic.site = data.sub$METASTATIC_SITE
        panel = data.sub$GENE_PANEL[which(data.sub$SAMPLE_ID == chosen.sample)]
        sample_number = nrow(data.sub)
        
        #' if multiple samples; chose highest purity otherwise highest coverage
      } else {
        sample_pool = data.sub[!is.na(data.sub$TUMOR_PURITY), ]
        chosen.sample = ifelse(nrow(sample_pool) == 0,
                               data.sub$SAMPLE_ID[which.max(data.sub$SAMPLE_COVERAGE)],
                               ifelse(nrow(sample_pool) > 1 & length(unique(sample_pool$TUMOR_PURITY)) == 1,
                                      sample_pool$SAMPLE_ID[1], sample_pool$SAMPLE_ID[which.max(sample_pool$TUMOR_PURITY)]))
        
        sample.type = data.sub$SAMPLE_TYPE[which(data.sub$SAMPLE_ID == chosen.sample)]
        metastatic.site = data.sub$METASTATIC_SITE[which(data.sub$SAMPLE_ID == chosen.sample)]
        panel = data.sub$GENE_PANEL[which(data.sub$SAMPLE_ID == chosen.sample)]
        sample_number = nrow(data.sub)
      }
      
      data.out = data.frame(PATIENT_ID = i,
                            SAMPLE_ID = chosen.sample,
                            PANEL = panel,
                            SAMPLE_TYPE = sample.type,
                            METASTATIC_SITE = metastatic.site,
                            SAMPLE_NUMBER = sample_number)
      
      data.analysis = rbind(data.analysis, data.out)
    })
  }
  
  # clean up
  rm(data.sub)
  rm(chosen.sample)
  rm(data.out)
  return(data.analysis)
}

MSK_one_pts_sample = sample_selection(data = portal_data)


##-----------------
## Get count-matrices
##-----------------
MSK_one_pts_sample = merge(MSK_one_pts_sample, 
                           FacetsPaths[, c('tumor_sample', 'counts_file')],
                           by.x = 'SAMPLE_ID', 
                           by.y = 'tumor_sample', all.x = T)

sh = MSK_one_pts_sample[is.na(MSK_one_pts_sample$counts_file), ]
mi = read.csv('Data/00_CohortData/Countmatrices_Facets.txt', sep = '\t')
sh_mi = merge(sh, mi, by.x = 'SAMPLE_ID', by.y = 'id', all.x = T)
sh_mi$counts_file = NULL
colnames(sh_mi)[ncol(sh_mi)] = 'counts_file'

MSK_f = MSK_one_pts_sample[!MSK_one_pts_sample$SAMPLE_ID %in% sh$SAMPLE_ID, ]
MSK_out = rbind(MSK_f, sh_mi)
MSK_out = MSK_out[!is.na(MSK_out$counts_file), ]


##----------------+
## Get bam-file paths;
##----------------+
keys = read.csv('~/Documents/MSKCC/dmp-2021/Genomics/keys.txt', sep = '\t')
MSK_out$counts_file = gsub(pattern = '//', replacement = '/', MSK_out$counts_file)
MSK_out$Normal_file = NA
for(i in 1:nrow(MSK_out)){
  print(i)
  MSK_out$Normal_file[i] = substr(x = strsplit(MSK_out$counts_file[i], split = '/')[[1]][11], 
                                  start = 19, stop = 35)
}

path = '/juno/res/dmpcollab/dmpshare/share/irb12_245/'
MSK_out$Normal_bam = NA
for(i in 1:nrow(MSK_out)){
  print(i)
  if(MSK_out$Normal_file[i] %in% keys$DMPID){
    MSK_out$Normal_bam[i] = paste0(path, substr(x = keys$AnonymizedID[which(keys$DMPID == MSK_out$Normal_file[i])], start = 1, stop = 1),
                               '/', substr(x = keys$AnonymizedID[which(keys$DMPID == MSK_out$Normal_file[i])], start = 2, stop = 2),
                               '/', keys$AnonymizedID[which(keys$DMPID == MSK_out$Normal_file[i])], '.bam')
  } else {
    MSK_out$Normal_bam[i] = NA
  }
}


##----------------+
## Clinical variables
##----------------+
Clinical_annotation = portal_data
MSK_out$PANEL = NULL
MSK_out$SAMPLE_TYPE = NULL
MSK_out$PATIENT_ID = NULL
MSK_out$METASTATIC_SITE = NULL
MSK_out$SAMPLE_NUMBER = NULL

MSK_out = merge(MSK_out, Clinical_annotation, by.x = 'SAMPLE_ID', by.y = 'SAMPLE_ID', all.x = T)
colnames(MSK_out)[25] = 'OS_Status'
colnames(MSK_out)[26] = 'OS_months'
colnames(MSK_out)[24] = 'Age_Sequencing'

MSK_out = MSK_out[,c(1,5, 49, 50, 6,7,8,9,10,11, 12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,2,3,4)]
write.table(MSK_out, 
            "~/Documents/MSKCC/10_MasterThesis/Data/00_CohortData/IMPACT_dataFreeze_07.13.22.txt", 
            sep = "\t", row.names = F, quote = F)

Cohort071322 = list(MSK_IMPACT_cohort = MSK_out)
saveRDS(Cohort071322, file = 'Data/00_CohortData/Cohort_071322.rds')




##----------------+
## Mosaic annotation: 
## 08/24/2022
##----------------+
cohort = readRDS('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')
clinical = cohort$IMPACT_clinicalAnnotation
IMPACT_LOY = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/03_Mosaicism/IMPACT_LOY.txt', sep = '\t')


all_out = data.frame()
for(i in 1:nrow(clinical)){
  if(clinical$SAMPLE_ID[i] %in% IMPACT_LOY$sample.id){
    sample = clinical$SAMPLE_ID[i]
    mLRR_Y = IMPACT_LOY$corrected.CN[which(IMPACT_LOY$sample.id == clinical$SAMPLE_ID[i])]
    LOY = IMPACT_LOY$LOY[which(IMPACT_LOY$sample.id == clinical$SAMPLE_ID[i])]
  } else {
    sample = clinical$SAMPLE_ID[i]
    mLRR_Y = 'N/A'
    LOY = 'N/A'
  }
  out = data.frame(sample = sample,
                   mLRR_Y = mLRR_Y,
                   LOY = LOY)
  all_out = rbind(all_out, out)
}

all_out = all_out[!duplicated(all_out$sample), ]
all_out = merge(all_out, clinical[,c('SAMPLE_ID', 'MSI_TYPE', 'AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)',
                                     'OVERALL_SURVIVAL_STATUS', 'OVERALL_SURVIVAL_(MONTHS)')],
                by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
colnames(all_out)[5] = 'Age_Sequencing'
colnames(all_out)[6] = 'OS_Status'
colnames(all_out)[7] = 'OS_months'

#' fetch and save
cohort = list(IMPACT_cohort = cohort$IMPACT_cohort,
              IMPACT_clinicalAnnotation = cohort$IMPACT_clinicalAnnotation,
              IMPACT_LOY = all_out)
saveRDS(cohort, file = '~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')




##-----------------
## Chromosome Y loss
##-----------------
cohort = readRDS('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')
IMPACT = cohort$IMPACT_cohort
LOY = cohort$IMPACT_LOY
IMPACT = merge(IMPACT, LOY[,c('sample', 'LOY')], by.x = 'SAMPLE_ID', by.y = 'sample', all.x = T)

cohort = list(IMPACT_cohort = IMPACT,
              IMPACT_clinicalAnnotation = cohort$IMPACT_clinicalAnnotation,
              IMPACT_LOY = cohort$IMPACT_LOY)
saveRDS(cohort, file = '~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')


##-----------------
## data merge
##-----------------
cohort = readRDS(file = 'Data/signedOut/Cohort_07132022.rds')
impact = cohort$IMPACT_cohort
clinical = cohort$IMPACT_clinicalAnnotation

impact = merge(impact, clinical[, c('SAMPLE_ID', 'CANCER_TYPE', 'CANCER_TYPE_DETAILED')], by = 'SAMPLE_ID', all.x = T)

LOY = cohort$IMPACT_binaryY_call
impact = merge(impact, LOY[,c('sample.id', 'Y_call')], by.x = 'SAMPLE_ID', by.y = 'sample.id', all.x = T)

Cohort = list(IMPACT_cohort = impact,
              IMPACT_clinicalAnnotation = cohort$IMPACT_clinicalAnnotation,
              IMPACT_LOY = cohort$IMPACT_LOY,
              IMPACT_binaryY_call = cohort$IMPACT_binaryY_call)

saveRDS(Cohort, file = '~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')


##-----------------
## precise Y-status classification
##-----------------
cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
QC_true = cohort$IMPACT_binaryY_call$sample.id
Y_CNA = read.csv('Data/04_Loss/Categorial_Classification_Y.txt', sep = '\t')
Y_CNA = Y_CNA[which(Y_CNA$sample %in% QC_true), ]

Cohort = list(IMPACT_cohort = cohort$IMPACT_cohort,
              IMPACT_clinicalAnnotation = cohort$IMPACT_clinicalAnnotation,
              IMPACT_LOY = cohort$IMPACT_LOY,
              IMPACT_binaryY_call = cohort$IMPACT_binaryY_call,
              IMPACT_Y_classification_final = Y_CNA)

saveRDS(Cohort, file = '~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')


##----------------+
## Genome-alteration annotation
##----------------+
cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
Cohort = list(IMPACT_cohort = cohort$IMPACT_cohort,
              IMPACT_clinicalAnnotation = cohort$IMPACT_clinicalAnnotation,
              IMPACT_LOY = cohort$IMPACT_LOY,
              IMPACT_binaryY_call = cohort$IMPACT_binaryY_call,
              IMPACT_Y_classification_final = cohort$IMPACT_Y_classification_final,
              IMPACT_ARM_level_changes = arm_changes_full_cohort)

saveRDS(Cohort, file = '~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')



#######################################
#######################################
#' prepare TCGA data
TCGA_clinical_annotation = read.csv('Data_in/TCGA_Clinical_Annotation.txt', sep = '\t')
TCGA_filepaths = read.csv('Data_in/TCGA_filepaths.txt', sep = '\t', header = F)
TCGA_filepaths$sample_id = substr(TCGA_filepaths$V1, start = 1, stop = 12)
TCGA_male = TCGA_clinical_annotation[which(TCGA_clinical_annotation$gender == 'MALE'), c('bcr_patient_barcode', 'type', 'age_at_initial_pathologic_diagnosis', 'gender') ]
TCGA_cohort = TCGA_male[which(TCGA_male$bcr_patient_barcode %in% TCGA_filepaths$sample_id), ]

#' use the first TCGA sample; regardless of number of samples
#' discard all additional samples
TCGA_filepaths$keep = NA
for(i in unique(TCGA_filepaths$sample_id)){
  if(length(TCGA_filepaths$sample_id[which(TCGA_filepaths$sample_id == i)]) > 1){
    TCGA_filepaths[which(TCGA_filepaths$sample_id == i), 'keep'][1] = 'keep'
    TCGA_filepaths[which(TCGA_filepaths$sample_id == i), 'keep'][2:length(TCGA_filepaths$sample_id[which(TCGA_filepaths$sample_id == i)])] = 'discard'
  } else {
    TCGA_filepaths[which(TCGA_filepaths$sample_id == i), 'keep'] = 'keep'
  }
}

TCGA_filepaths = TCGA_filepaths[which(TCGA_filepaths$keep %in% 'keep'), ]
TCGA_cohort = merge(TCGA_cohort, TCGA_filepaths, by.x = 'bcr_patient_barcode', by.y = 'sample_id', all.x = T)


#' create processing Y-chromosome loss pipeline; n = 4,875 male samples TCGA
TCGA_cohort$path = paste0('/juno/work/ccs/tcga_facets/snp-pileup/', TCGA_cohort$V1)
TCGA_cohort$V1 = NULL
write.table(TCGA_cohort, file = 'Data_in/TCGA_paths.txt', sep = '\t', row.names = F, quote = F)


##' create common data-container for distribution and sharing
TCGA_data = read.csv('Data_in/TCGA_paths.txt', sep = '\t')
cohortData = readRDS('Data_out/cohort_data.rds')
cohortData = append(cohortData, list(TCGA_data))
names(cohortData)[5] = 'male_TCGA_cohort'
cohortData$male_TCGA_cohort
saveRDS(cohortData, file = 'Data_out/cohort_data.rds')


#' binary calls in TCGA
TCGA_binary = read.csv('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/TCGA/TCGA.binaryLoss_out.txt', sep = '\t')
TCGA = a$male_TCGA_cohort

TCGA_full = merge(TCGA, TCGA_binary[,c('Y_call', 'sample.id', 'ploidy', 'purity')],
                  by.x = 'bcr_patient_barcode', by.y = 'sample.id', all.x = T)

cohortData = readRDS('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/cohort_data.rds')
cohortData = append(cohortData, list(TCGA_full))
names(cohortData)[6] = 'TCGA.cohort'
cohortData$male_TCGA_cohort = NULL
saveRDS(cohortData, file = '~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/cohort_data.rds')


#' add mosaicism to cohort
#' MSK_WES_GermlineCN = read.csv('Data_out/WES/MSK_WES_GermlineCN.txt', sep = '\t')
sample_summary$id = substr(sample_summary$sample, start = 3, stop = 17)
MSK_WES = sample_summary[which(sample_summary$target == 24), c('target', 'corrected.CN', 'id', 'Mosaic')]
colnames(MSK_WES)[2] = 'corrected_GermlineCN'
colnames(MSK_WES)[4] = 'Mosaic'
MSK_WES$Mosaic[which(MSK_WES$Mosaic == 'intakt')] = 'no'
MSK_WES$Mosaic[which(MSK_WES$Mosaic == 'mosaic')] = 'yes'

#' merge with current data cohort
WES_cohort = readRDS('Data_out/cohort_data.rds')$WES.cohort
WES_cohort$Y_mosaic = NA

for(i in 1:nrow(WES_cohort)){
  if(WES_cohort$CMO_Sample_ID[i] %in% MSK_WES$id){
    WES_cohort$Y_mosaic[i] = MSK_WES$Mosaic[which(MSK_WES$id == WES_cohort$CMO_Sample_ID[i])]
  } else {
    WES_cohort$Y_mosaic[i] = NA
  }
}

#' add age annotation to MSK-WES;
cohort = readRDS('Data_out/cohort_data.rds')
age = cohort$IMPACT.cohort[,c('SAMPLE_ID', 'AGE_AT_SEQ_REPORTED_YEARS')]
WES_cohort = merge(WES_cohort, age, by.x = 'DMP_Sample_ID', by.y = 'SAMPLE_ID', all.x = T)


cohort$WES.cohort = WES_cohort
saveRDS(cohort, file = '~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/cohort_data.rds')


#' clinical associations:
clinical = read.csv('Data_out/sub_clinical.tsv', sep = '\t')
clinical = clinical[, c('Sample.ID', 'Overall.Survival..Months.', 'Overall.Survival.Status')]
colnames(clinical) = c('Sample_ID', 'OS', 'OS_Status')

MSK_WES = merge(MSK_WES, clinical, by.x = 'DMP_Sample_ID', by.y = 'Sample_ID', all.x = T)
cohort$WES.cohort = MSK_WES
saveRDS(cohort, file = '~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/cohort_data.rds')


#' add IMPACT MOSAIC data to cohort
cohort = append(cohort, list(IMPACT_age))
names(cohort)[6] = 'IMPACT.mosaic'


#' add IMPACT CLINICAL data to cohort
cohort = append(cohort, list(clinical))
names(cohort)[7] = 'IMPACT.clinical'
saveRDS(cohort, file = '~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/cohort_data.rds')

