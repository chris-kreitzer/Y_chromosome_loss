## Data selection for Y-chromosome loss determination
## 
## Read in clinical data from cBioPortal: 07/11/2022
## This should be at the patient level  
## Need to select one sample per patient 
# 
# use_for_patient_level (T/F): column indicating which of the tumor samples to use for patient-level analyses. Relevant when a patient has multiple biopsies sequenced
# use_for_sample_level (T/F): column indicating whether sample passed basic purity check.
##
## start: 07/11/2022
## chris-kreitzer

clean()
gc()
setwd('~/Documents/MSKCC/10_MasterThesis/')

library(dplyr)
library(data.table)
library(usefun)

# load clinical and merge data_clinical_oncoKB with data_from the portal
portal_data = fread('Data/mskimpact_clinical_data_07112022.tsv', stringsAsFactors = F, data.table = F)
colnames(portal_data) = toupper(colnames(portal_data))
colnames(portal_data) = gsub(pattern = ' ', replacement = '_', colnames(portal_data), fixed = T )
rownames(portal_data) = portal_data$SAMPLE_ID
FacetsPaths = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/facets_annotated.cohort.txt.gz', sep = '\t')
FacetsPaths = FacetsPaths[!duplicated(FacetsPaths$tumor_sample), ]

# exclude heme patients
portal_data$STUDY_ID = NULL
portal_data$DIAGNOSIS_AGE = NULL
portal_data = portal_data[!portal_data$GENE_PANEL %in% c("ACCESS129", "IMPACT-HEME-400", "IMPACT-HEME-468"), ]
portal_data$`AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)` = as.character(as.numeric(portal_data$`AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)`))
portal_data$TUMOR_PURITY = as.numeric(as.character(portal_data$TUMOR_PURITY))
# exclude patients with MSI-I or TMB >= 20 #
# portal_data = portal_data[!portal_data$MSI_TYPE %in% 'Instable', ]
# portal_data = portal_data[which(portal_data$IMPACT_TMB_SCORE <= 20),]

# concentrate on males:
portal_data = portal_data[which(portal_data$SEX == 'Male'), ]

##-----------------
# one sample per patient:

## function which selects ONE sample per patient
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
  
  # clean up the mess
  rm(data.sub)
  rm(chosen.sample)
  rm(data.out)
  
  return(data.analysis)
}

MSK_one_pts_sample = sample_selection(data = portal_data)


##-----------------
## Get countFiles
MSK_one_pts_sample = merge(MSK_one_pts_sample, FacetsPaths[, c('tumor_sample', 'counts_file')],
      by.x = 'SAMPLE_ID', by.y = 'tumor_sample', all.x = T)


##-----------------
## missing countfiles
# missing_countfiles = MSK_one_pts_sample[is.na(MSK_one_pts_sample$counts_file), 'SAMPLE_ID']
# write.table(missing_countfiles, file = 'Data/PanCancerMissingCountfiles.txt', row.names = F, quote = F, sep = '\t')
# 
# all_missing = data.frame()
# for(i in unique(missing_countfiles)){
#   try({
#     common = '/juno/work/ccs/shared/resources/impact/facets/all/'
#     path_high = list.files(path = paste0(common, substr(i, start = 1, stop = 7), '/'))
#     path_lower = path_high[grep(pattern = i, path_high)]
#     path_folder = list.files(path = paste0(common, substr(i, start = 1, stop = 7), '/', path_lower, '/'), full.names = T)
#     path_counts = grep(pattern = 'countsMerged', path_folder, value = T)
#     out = data.frame(SAMPLE_ID = i,
#                      counts_file = path_counts)
#     
#     all_missing = rbind(all_missing, out)
#   })
# }
# 
# missing_paths = read.csv('Data/missing.txt', sep = '\t')
# missing_paths = missing_paths[!duplicated(missing_paths$SAMPLE_ID), ]

MSK_one_pts_sample = merge(MSK_one_pts_sample, missing_paths, by.x = 'SAMPLE_ID', by.y = 'SAMPLE_ID', all.x = T)
MSK_IMPACT_cohort = MSK_one_pts_sample[!is.na(MSK_one_pts_sample$counts_file.x), ]
MSK_IMPACT_cohort$counts_file.y = NULL
colnames(MSK_IMPACT_cohort)[7] = 'counts_file'

##-------------------------------------
write.table(MSK_IMPACT_cohort, "Data/IMPACT_dataFreeze_07.13.22.txt", sep = "\t", row.names = F, quote = F)


Cohort = list(IMPACT_cohort = MSK_IMPACT_cohort)
saveRDS(Cohort, file = 'Data/signedOut/Cohort_07132022.rds')













#' make an automated data selection process
data_freeze1 = read.csv('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/data_freeze_cancer_type_082321.txt', sep = '\t')
data_WES = read.csv('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/WES_cBio_ID.match.txt', sep = '\t')
Facets_annotation = read.csv('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/Facets_annotated.cohort.txt', sep = '\t')
Facets_Countsfile = Facets_annotation[, c('tumor_sample', 'counts_file')]
data_freeze1 = merge(data_freeze1, Facets_Countsfile, by.x = 'SAMPLE_ID', by.y = 'tumor_sample', all.x = T)
data_freeze1 = data_freeze1[-which(duplicated(data_freeze1$SAMPLE_ID)), ]

#' WES countsfile: working with n = 957
WES_data = data_WES[which(data_WES$DMP_Sample_ID %in% data_freeze1$SAMPLE_ID), ]



#' save output; for generell purpose
Y_chromosome_cohort_data = list(IMPACT.cohort = data_freeze1,
                                WES.cohort = WES_data)
saveRDS(object = Y_chromosome_cohort_data, file = '~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/cohort_data.rds')



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

