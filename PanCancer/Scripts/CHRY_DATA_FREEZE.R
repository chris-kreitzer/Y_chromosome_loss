# Read in clinical data
# This should be at the patient level 
# Need to select one sample per patient
# use_for_patient_level (T/F): column indicating which of the tumor samples to use for patient-level analyses. Relevant when a patient has multiple biopsies sequenced
# use_for_sample_level (T/F): column indicating whether sample passed basic purity check.

pat_df = read.delim("~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/IMPACT_clinical_annotation.tsv")
pat_df = pat_df %>% 
  mutate_if(is.factor, as.character) 
samp_df <- read.delim("~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/data_clinical_sample.oncokb.txt")
samp_df <- samp_df %>% 
  mutate_if(is.factor, as.character) %>%
  filter(use_for_patient_level == TRUE)
samp_df 

data_freeze = samp_df %>% 
  dplyr::select(PATIENT_ID, SAMPLE_ID, CANCER_TYPE, facets_qc, TUMOR_PURITY, 
                SAMPLE_TYPE, AGE_AT_SEQ_REPORTED_YEARS, GENE_PANEL, SAMPLE_COVERAGE) %>%
  left_join(select(pat_df, PATIENT_ID, Sex), by = "PATIENT_ID") %>%
  filter(Sex != "", Sex != "Female", facets_qc == "TRUE" ) %>% distinct()

write.table(data_freeze, "~/Downloads/data_freeze_cancer_type_082321.txt", sep = "\t", row.names = F, quote = F)



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






