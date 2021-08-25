# Read in clinical data
# This should be at the patient level 
# Need to select one sample per patient
# use_for_patient_level (T/F): column indicating which of the tumor samples to use for patient-level analyses. Relevant when a patient has multiple biopsies sequenced
# use_for_sample_level (T/F): column indicating whether sample passed basic purity check.

pat_df = read.delim("~/Downloads/mskimpact_clinical_data (65).tsv")
pat_df = pat_df %>% 
  mutate_if(is.factor, as.character) 
samp_df <- read.delim("~/Desktop/Projects/40K_Analysis/data/data_clinical_sample.oncokb.txt")
samp_df <- samp_df %>% 
  mutate_if(is.factor, as.character) %>%
  filter(use_for_patient_level == TRUE)
samp_df 

data_freeze <- samp_df %>% 
  dplyr::select(PATIENT_ID, SAMPLE_ID, CANCER_TYPE, facets_qc,TUMOR_PURITY,SAMPLE_TYPE,AGE_AT_SEQ_REPORTED_YEARS,GENE_PANEL,SAMPLE_COVERAGE) %>%
  left_join(select(pat_df, PATIENT_ID, Sex), by = "PATIENT_ID")%>%
  filter(Sex != "",Sex !="Female", facets_qc == "TRUE" )%>% distinct()

write.table(data_freeze, "~/Downloads/data_freeze_cancer_type_082321.txt", sep = "\t", row.names = F, quote = F)



#' make an automated data selection process
data_freeze1 = read.csv('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/data_freeze_cancer_type_082321.txt', sep = '\t')
data_WES = read.csv('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/WES_cBio_ID.match.txt', sep = '\t')
Facets_annotation = read.csv('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/Facets_annotated.cohort.txt', sep = '\t')
Facets_Countsfile = Facets_annotation[, c('tumor_sample', 'counts_file')]
data_freeze1 = merge(data_freeze1, Facets_Countsfile, by.x = 'SAMPLE_ID', by.y = 'tumor_sample', all.x = T)
data_freeze1 = data_freeze1[-which(duplicated(data_freeze1$SAMPLE_ID)), ]

#' save output; for generell purpose
Y_chromosome_cohort_data = list(IMPACT.cohort = data_freeze1)
saveRDS(object = Y_chromosome_cohort_data, file = '~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/cohort_data.rds')
