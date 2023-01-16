##----------------+
## Survival analysis of LOY
## and cancer-types;
## Run Cox model and filter out
## cancer types with significant
## association of deprived OS;
##----------------+
## 
## start: 08/16/2021
## revision: 08/24/2022
## revision: 01/16/2023
## 
## chris-kreitzer
## 
## 
## TODO:
## - QC TRUE/FALSE samples
## - complete_loss as TRUE? what's with partial loss?

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
source('Scripts/UtilityFunctions.R')

library(readxl)
library(data.table)
library(survival)




##-------
## Data rendering;
##-------
Cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
OS_cohort = Cohort[,c('PATIENT_ID', 'QC', 'Age_Sequencing', 'genome_doubled', 'SAMPLE_TYPE', 'MSI_TYPE',
                      'classification', 'OS_Status', 'OS_months')]
OS_cohort = OS_cohort[!is.na(OS_cohort$OS_Status) & !is.na(OS_cohort$OS_months), ]
colnames(OS_cohort)[7] = 'LOY'
OS_cohort$LOY = ifelse(OS_cohort$LOY == 'complete_loss', TRUE, FALSE)
OS_cohort$OS_Status = ifelse(OS_cohort$OS_Status == '1:DECEASED', 'Dead', 'Alive')
OS_cohort$OS_Status_INT = ifelse(OS_cohort$OS_Status == 'Dead', 2L, 1L)


##-------
## COX-proportional hazard model;
## PanCancer: not accounting for co-variates
##-------
loy.pancancer = coxph(formula = Surv(OS_months, as.numeric(factor(OS_Status))) ~ LOY, data = OS_cohort)
summary(loy.pancancer)
coxph(Surv(OS_months, OS_Status_INT) ~ LOY, data = OS_cohort) %>%
  gtsummary::tbl_regression(exp = TRUE)



##-------
## Kaplan-Meier method
##-------
loy.pancancer_KM = survfit(Surv(OS_months, OS_Status_INT) ~ LOY, data = OS_cohort)
p1 = ggsurvplot(loy.pancancer_KM,
                palette = c("#3d4397", "#a22231"),
                size = 0.8,
                legend = 'top',
                xlab = 'Time (months)',
                ylab = 'Overall survival (%)',
                )$plot +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(aspect.ratio = 1,
        legend.position = 'none')
 
ggsave_golden(filename = 'Figures_original/OS_wholeCohort.pdf', plot = p1, width = 6)
summary(loy.pancancer)




##----------------+
## Adjustment of CoxP model
## for age and tumor site;
## 
## we can potentially include
## more variables (later)
##----------------+
loy.pancancer_MV = coxph(formula = Surv(OS_months, as.numeric(factor(OS_Status))) ~ LOY+genome_doubled+SAMPLE_TYPE, data = OS_cohort)
summary(loy.pancancer_MV)


















