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

library(readxl)
library(data.table)
library(survival)




##-------
## Data rendering;
##-------
Cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
OS_cohort = Cohort[,c('PATIENT_ID', 'QC', 'Age_Sequencing', 'classification', 'OS_Status', 'OS_months')]
colnames(OS_cohort)[4] = 'LOY'
OS_cohort$LOY = ifelse(OS_cohort$LOY == 'complete_loss', TRUE, FALSE)
OS_cohort$OS_Status = ifelse(OS_cohort$OS_Status == '1:DECEASED', 'Dead', 'Alive')


loy.pancancer = coxph(formula = Surv(OS_months, as.numeric(factor(OS_Status))) ~ LOY, data = OS_cohort)
summary(loy.pancancer)


# Load per-patient WGD and outcome data
os.dat <-
  as.data.table(readxl::read_excel(
    "~/Downloads/NIHMS970114-supplement-1.xlsx",
    skip = 2,
    n_max = 9692
  ))
colnames(os.dat)[7] = 'OS'
colnames(os.dat)[8] = 'Use'
colnames(os.dat)[6] = 'VitalStatus'
model.pancan <-
  coxph(formula = Surv(OS, as.numeric(factor(VitalStatus))) ~ WGD,
        data = os.dat[Use == T])

summary(model.pancan)




Cohort = readRDS('Data/signedOut/Cohort_07132022.rds')


##-----------------
## LOY:Age distribution
##-----------------
LOY = Cohort$IMPACT_LOY
LOY$Age_Sequencing = as.integer(as.character(LOY$Age_Sequencing))

hist(LOY$Age_Sequencing[which(LOY$LOY == 'yes')], nclass = 30,
     xlim = c(20, 78),
     col = 'white',
     yaxt = 'n',
     xlab = '',
     ylab = '',
     main = '')
axis(side = 2, las = 2)
box(lwd = 2)
mtext(text = 'Frequency', side = 2, line = 2)
mtext(text = 'Age-distribution LOY', side = 3, line = 1, adj = 0)


##-----------------
## Survival analysis
##-----------------
library(survminer)
library(survival)

IMPACT_survival = LOY[!is.na(LOY$LOY) & !is.na(LOY$OS_Status), ]
IMPACT_survival = IMPACT_survival[!IMPACT_survival$LOY == 'N/A', ]
IMPACT_survival$OS_Status[which(IMPACT_survival$OS_Status == '1:DECEASED')] = 2L
IMPACT_survival$OS_Status[which(IMPACT_survival$OS_Status == '0:LIVING')] = 1L
IMPACT_survival$OS_Status = as.numeric(as.character(IMPACT_survival$OS_Status))
IMPACT_survival$LOY = factor(IMPACT_survival$LOY, levels = c('no', 'yes'))


coxph(Surv(OS_months, OS_Status) ~ LOY + Age_Sequencing, data = IMPACT_survival) %>%
  gtsummary::tbl_regression(exp = TRUE)


#' survival analysis
f1 = survfit(Surv(OS_months, OS_Status) ~ LOY + Age_Sequencing, data = IMPACT_survival)
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
               palette = 'aaas')
IMPACT_LOY_survival


#' out



##-----------------
## Old cohort: 2021
##-----------------

#' cohort = readRDS('Data_out/cohort_data.rds')
#' MSK_WES = cohort$WES.cohort
#' MSK_WES$AGE_AT_SEQ_REPORTED_YEARS = as.numeric(as.character(MSK_WES$AGE_AT_SEQ_REPORTED_YEARS))
#' MSK_WES = MSK_WES[!is.na(MSK_WES$AGE_AT_SEQ_REPORTED_YEARS), ]
#' 
#' 
#' mosaic_age = data.frame()
#' for(i in seq(5, 90, 5)){
#'   da = MSK_WES[between(MSK_WES$AGE_AT_SEQ_REPORTED_YEARS, i-4, i), ]
#'   
#'   if(length(table(da$Y_mosaic)) == 2){
#'     ratio = (table(da$Y_mosaic)[2] / sum(table(da$Y_mosaic))) * 100
#'     n = sum(table(da$Y_mosaic))
#'     out = data.frame(group = i,
#'                      ratio = ratio,
#'                      n = n)
#'     
#'   } else if (length(table(da$Y_mosaic)) == 1) {
#'     table.ratio = table(da$Y_mosaic)
#'     ratio = ifelse(names(table.ratio) == 'no', 0, 100)
#'     n = table.ratio[[1]]
#'     out = data.frame(group = i,
#'                      ratio = ratio,
#'                      n = n)
#'   } else {
#'     ratio = NA
#'     n = dim(da)[[1]]
#'     out = data.frame(group = i,
#'                      ratio = ratio, 
#'                      n = n)
#'   }
#'   mosaic_age = rbind(mosaic_age, out)
#' }
#'   
#' 
#' #' make a visualization
#' age_mosaic = ggplot(mosaic_age[!mosaic_age$group %in% c(5, 10, 15, 20), ], aes(x = group, y = ratio)) + 
#'   geom_hline(yintercept = seq(0, 50, 25), linetype = 'dashed', size = 0.2, color = 'grey35') +
#'   geom_bar(stat = 'identity') +
#'   stat_smooth(method = 'gam', formula = y ~ s(x, k = 3), size = 2, se = F, color = 'red', alpha = 0.2) +
#'   scale_x_continuous(expand = c(0, 0),
#'                      breaks = seq(25, 90, 5),
#'                      labels = paste0(mosaic_age[!mosaic_age$group %in% c(5, 10, 15, 20), 'group'], '\n',
#'                                      ' n=',
#'                                      mosaic_age[!mosaic_age$group %in% c(5, 10, 15, 20), 'n'])) +
#'   scale_y_continuous(expand = c(0, 0),
#'                      limits = c(0, 60),
#'                      breaks = c(0, 50, 25)) +
#'   
#'   theme_bw(base_size = 22) +
#'   theme(panel.grid = element_blank()) +
#'   labs(x = 'AGE groups', y = 'Percent mLOY')
#'   
#' ggsave_golden(filename = 'Figures/WES_Age_Mosaic.pdf', plot = age_mosaic, width = 12)
#' 
#' 
#' 
#' #' look into overall survival associated with mosaic loss of Y-chromosome
#' WES_mosaic = MSK_WES[!is.na(MSK_WES$Y_mosaic), ]
#' 
#' #' survival analysis:
#' MSK_WES = readRDS('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/cohort_data.rds')$WES.cohort
#' 
#' library(survminer)
#' library(survival)
#' 
#' WES_survival = MSK_WES[!is.na(MSK_WES$Y_mosaic) & !is.na(MSK_WES$OS_Status), ]
#' WES_survival$OS_Status[which(WES_survival$OS_Status == '1:DECEASED')] = 2L
#' WES_survival$OS_Status[which(WES_survival$OS_Status == '0:LIVING')] = 1L
#' WES_survival$OS_Status = as.numeric(as.character(WES_survival$OS_Status))
#' WES_survival$Y_mosaic = factor(WES_survival$Y_mosaic, levels = c('no', 'yes'))
#' 
#' 
#' coxph(Surv(OS, OS_Status) ~ Y_mosaic, data = WES_survival) %>%
#'   gtsummary::tbl_regression(exp = TRUE)
#' 
#' #' survival analysis
#' f1 = survfit(Surv(OS, OS_Status) ~ Y_mosaic, data = WES_survival)
#' WES_survival$age = ifelse(WES_survival$AGE_AT_SEQ_REPORTED_YEARS < 50, 'group1', 
#'                           ifelse(WES_survival$AGE_AT_SEQ_REPORTED_YEARS > 50 & WES_survival$AGE_AT_SEQ_REPORTED_YEARS < 70, 'group2', 'group3')) 
#'                                  
#' f2 = survfit(Surv(OS, OS_Status) ~ Y_mosaic + age, data = WES_survival)
#' 
#' #' make the plot
#' a = ggsurvplot(f1,
#'            size = 0.8,
#'            censor.shape = '',
#'            pval = T,
#'            risk.table = F,
#'            risk.table.height = 0.25,
#'            break.time.by = 12,
#'            ggtheme = theme_light(),
#'            surv.median.line = 'v',
#'            xlim = c(0, 84),
#'            xlab = ('OS [months]'),
#'            palette = 'aaas')
#'   
#' b = a$plot + theme(panel.background = element_blank()) +
#'   theme(aspect.ratio = 0.5,
#'         panel.border = element_rect(size = 1.5, colour = 'black'),
#'         axis.text = element_text(size = 14, color = 'black', face = 'bold'),
#'         legend.position = 'top', 
#'         legend.text = element_text(size = 14)) +
#'   scale_color_discrete(name = "", labels = c("full Y-chromosome", "mosaic Y loss"), 
#'                        type = c('#414A8F', '#CD2703'))
#' 
#' ggsave_golden(filename = 'Figures/OS_Mosaic_YLoss.pdf', plot = b, width = 9)
#' 
#' 
#' 
#' 
#' ###############################################################################
#' #' Survival analysis for IMPACT cases
#' #' modify input data first
#' clinical = read.csv('Data_out/IMPACT/IMPACT_clinical.txt', sep = '\t')
#' IMPACT_mosaic = read.csv('Data_out/IMPACT/GermlineCN.txt', sep = '\t')
#' IMPACT_mosaic = IMPACT_mosaic[which(IMPACT_mosaic$target == 'Y'), ]
#' 
#' OS_Impact = merge(IMPACT_mosaic, clinical[, c('SAMPLE_ID', 'FGA', 'OS_Months', 'OS_Status')], 
#'                   by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
#' 
#' 
#' #' OS analysis
#' library(survminer)
#' library(survival)
#' 
#' OS_Impact = OS_Impact[!is.na(OS_Impact$Y_mosaicism), ]
#' OS_Impact$Y_mosaicism = factor(OS_Impact$Y_mosaicism, levels = c('no', 'yes'))
#' 
#' #' cox-p model
#' summary(coxph(Surv(OS_Months, OS_Status) ~ Y_mosaicism, data = OS_Impact))
#' coxph(Surv(OS_Months, OS_Status) ~ Y_mosaicism, data = OS_Impact) %>%
#'   gtsummary::tbl_regression(exp = TRUE)
#' 
#' model_impact = survfit(Surv(OS_Months, OS_Status) ~ Y_mosaicism, data = OS_Impact)
#' 
#' #' make the plot
#' a = ggsurvplot(model_impact,
#'                size = 0.8,
#'                censor.shape = '',
#'                pval = T,
#'                risk.table = F,
#'                risk.table.height = 0.25,
#'                break.time.by = 12,
#'                ggtheme = theme_light(),
#'                surv.median.line = 'v',
#'                xlim = c(0, 84),
#'                xlab = ('OS [months]'),
#'                palette = 'aaas')
#' 
#' b = a$plot + theme(panel.background = element_blank()) +
#'   theme(aspect.ratio = 0.5,
#'         panel.border = element_rect(size = 1.5, colour = 'black'),
#'         axis.text = element_text(size = 14, color = 'black', face = 'bold'),
#'         legend.position = 'top', 
#'         legend.text = element_text(size = 14)) +
#'   scale_y_continuous(expand = c(0.01, 0.05)) +
#'   scale_x_continuous(expand = c(0.0, 0.0)) +
#'   scale_color_discrete(name = "", labels = c("full Y-chromosome", "mosaic Y loss"), 
#'                        type = c('#414A8F', '#CD2703'))
#' 
#' ggsave_golden(filename = 'Figures/OS_Mosaic_IMPACT.pdf', plot = b, width = 9)