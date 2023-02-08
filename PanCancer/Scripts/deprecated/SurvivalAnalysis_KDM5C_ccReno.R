##----------------+
## Survival Analysis based on
## Renal Cell Carcinoma
## KDM5C/LOY status
##----------------+

library(survminer)
library(survival)

cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
RenoCancer = read.csv('Data/05_Association/KDM5C_muts_RenalCancers.txt', sep = '\t')
copyNumbers = read.csv('Data/04_Loss/IMPACT_copynumber_out.txt', sep = '\t')
copyNumbers = copyNumbers[which(copyNumbers$ID %in% RenoCancer$sample), ]

##----------------+
## Assess KDM5C CNA status
##----------------+
KD = data.table(chrom = 23,
                start = 53220503,
                end = 53254604)
key.col = c('chrom', 'start', 'end')

chromosomeX = data.frame()
for(i in unique(copyNumbers$ID)){
  try({
    sample = i
    segs = as.data.table(copyNumbers[which(copyNumbers$ID == i), ])
    setkey(segs, chrom, start, end)
    overlap = foverlaps(KD, segs, by.x = key.col, by.y = key.col, nomatch = 0)
    patient_tcn = overlap$tcn.em
    out = data.frame(sample = sample,
                     chrom = 23,
                     tcn = patient_tcn)
    chromosomeX = rbind(chromosomeX, out)
  })
}

#' there are samples (chromosome X)
#' where two segments are called (delete
#' duplicated segments)
chromosomeX = chromosomeX[-522,]
chromosomeX = chromosomeX[-238,]
chromosomeX = chromosomeX[-140,]
chromosomeX = chromosomeX[-89,]

Reno = left_join(RenoCancer, chromosomeX)
Reno$ploidy = NULL
Reno$TP53 = NULL
Reno$VHL = NULL
Reno$chrom = NULL

RenoLOY = merge(cohort$IMPACT_Y_classification_final, cohort$IMPACT_clinicalAnnotation,
                by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
RenoLOY = RenoLOY[which(RenoLOY$CANCER_TYPE == 'Renal Cell Carcinoma'), ]

RenoKDM5C = merge(RenoLOY, Reno, by.x = 'sample', by.y = 'sample', all.x = T)
RenoKDM5C$Y_call.x = NULL
RenoKDM5C$ploidy = NULL
RenoKDM5C$SAMPLE_COVERAGE = NULL
RenoKDM5C$Y_expected = NULL
RenoKDM5C$PRIMARY_TUMOR_SITE = NULL
colnames(RenoKDM5C)[11] = 'OS_Status'
colnames(RenoKDM5C)[12] = 'OS_months'
colnames(RenoKDM5C)[13] = 'Y_call'
colnames(RenoKDM5C)[length(RenoKDM5C)] = 'tcn_X'

##----------------+
## Assign Status
##----------------+
RenoKDM5C$status = ifelse(RenoKDM5C$Y_call == 1 &
                            RenoKDM5C$KDM5C == 0 & 
                            RenoKDM5C$tcn_X > 0, 'mono',
                          ifelse(RenoKDM5C$Y_call == 0 &
                                   RenoKDM5C$KDM5C == 1 &
                                   RenoKDM5C$tcn_X > 0, 'mono',
                                 ifelse(RenoKDM5C$Y_call == 0 &
                                          RenoKDM5C$KDM5C == 0 &
                                          RenoKDM5C$tcn_X == 0, 'mono',
                                        ifelse(RenoKDM5C$Y_call == 1 &
                                                 RenoKDM5C$KDM5C == 0 &
                                                 RenoKDM5C$tcn_X == 0, 'bi',
                                               ifelse(RenoKDM5C$Y_call == 1 &
                                                        RenoKDM5C$KDM5C == 1 &
                                                        RenoKDM5C$tcn_X > 0, 'bi',
                                                      ifelse(RenoKDM5C$Y_call == 0 &
                                                               RenoKDM5C$KDM5C == 0 &
                                                               RenoKDM5C$tcn_X > 0, 'wt', NA))))))
                                 
RenoKDM5C$OS_Status[which(RenoKDM5C$OS_Status == '1:DECEASED')] = 2L
RenoKDM5C$OS_Status[which(RenoKDM5C$OS_Status == '0:LIVING')] = 1L
RenoKDM5C$OS_Status = as.numeric(as.character(RenoKDM5C$OS_Status))
RenoKDM5C = RenoKDM5C[!is.na(RenoKDM5C$OS_months), ]
RenoKDM5C$Y_call = factor(RenoKDM5C$Y_call, levels = c('0', '1'))
RenoKDM5C$status = factor(RenoKDM5C$status, levels = c('wt', 'mono', 'bi'))


coxph(Surv(OS_months, OS_Status) ~ status, data = RenoKDM5C) %>%
  gtsummary::tbl_regression(exp = TRUE)


#' survival analysis
f1 = survfit(Surv(OS_months, OS_Status) ~ status, data = RenoKDM5C)
IMPACT_LOY_survival = ggsurvplot(f1,
                                 size = 0.8,
                                 censor.shape = '',
                                 pval = T,
                                 font.x = c(12, 'bold'),
                                 font.y = c(12, 'bold'),
                                 font.caption = c(12, 'bold'),
                                 risk.table = T,
                                 risk.table.height = 0.25,
                                 font.tickslab = c(10, 'bold', 'black'),
                                 break.time.by = 12,
                                 ggtheme = theme_bw(),
                                 censor = F,
                                 xlim = c(0, 90),
                                 legend = 'top',
                                 xlab = ('OS [months]'))
IMPACT_LOY_survival
