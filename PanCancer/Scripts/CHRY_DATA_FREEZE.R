# Read in clinical data
# This should be at the patient level 
# Need to select on sample per patient
# use_for_patient_level (T/F): column indicating which of the tumor samples to use for patient-level analyses. Relevant when a patient has multiple biopsies sequenced
# use_for_sample_level (T/F): column indicating whether sample passed basic purity check.

pat_df <- read.delim("~/Downloads/mskimpact_clinical_data (65).tsv")
pat_df <- pat_df %>% 
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


data_freeze1  <- table(data_freeze$CANCER_TYPE,data_freeze$Sex)

