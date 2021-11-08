## Survival Analysis mosaic-Y-loss

setwd('~/Documents/GitHub/Y_chromosome_loss/PanCancer/')
rm(list = ls())
source('Scripts/UtilityFunctions.R')
library(dplyr)

cohort = readRDS('Data_out/cohort_data.rds')
MSK_WES = cohort$WES.cohort
MSK_WES$AGE_AT_SEQ_REPORTED_YEARS = as.numeric(as.character(MSK_WES$AGE_AT_SEQ_REPORTED_YEARS))
MSK_WES = MSK_WES[!is.na(MSK_WES$AGE_AT_SEQ_REPORTED_YEARS), ]


mosaic_age = data.frame()
for(i in seq(5, 90, 5)){
  da = MSK_WES[between(MSK_WES$AGE_AT_SEQ_REPORTED_YEARS, i-4, i), ]
  
  if(length(table(da$Y_mosaic)) == 2){
    ratio = (table(da$Y_mosaic)[2] / sum(table(da$Y_mosaic))) * 100
    n = sum(table(da$Y_mosaic))
    out = data.frame(group = i,
                     ratio = ratio,
                     n = n)
    
  } else if (length(table(da$Y_mosaic)) == 1) {
    table.ratio = table(da$Y_mosaic)
    ratio = ifelse(names(table.ratio) == 'no', 0, 100)
    n = table.ratio[[1]]
    out = data.frame(group = i,
                     ratio = ratio,
                     n = n)
  } else {
    ratio = NA
    n = dim(da)[[1]]
    out = data.frame(group = i,
                     ratio = ratio, 
                     n = n)
  }
  mosaic_age = rbind(mosaic_age, out)
}
  

#' make a visualization
age_mosaic = ggplot(mosaic_age[!mosaic_age$group %in% c(5, 10, 15, 20), ], aes(x = group, y = ratio)) + 
  geom_hline(yintercept = seq(0, 50, 25), linetype = 'dashed', size = 0.2, color = 'grey35') +
  geom_bar(stat = 'identity') +
  stat_smooth(method = 'gam', formula = y ~ s(x, k = 3), size = 2, se = F, color = 'red', alpha = 0.2) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(25, 90, 5),
                     labels = paste0(mosaic_age[!mosaic_age$group %in% c(5, 10, 15, 20), 'group'], '\n',
                                     ' n=',
                                     mosaic_age[!mosaic_age$group %in% c(5, 10, 15, 20), 'n'])) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 60),
                     breaks = c(0, 50, 25)) +
  
  theme_bw(base_size = 22) +
  theme(panel.grid = element_blank()) +
  labs(x = 'AGE groups', y = 'Percent mLOY')
  
ggsave_golden(filename = 'Figures/WES_Age_Mosaic.pdf', plot = age_mosaic, width = 12)



#' look into overall survival associated with mosaic loss of Y-chromosome
WES_mosaic = MSK_WES[!is.na(MSK_WES$Y_mosaic), ]



#' survival analysis:
MSK_WES = readRDS('Data_out/cohort_data.rds')$WES.cohort

library(survminer)
library(survival)

WES_survival = MSK_WES[!is.na(MSK_WES$Y_mosaic) & !is.na(MSK_WES$OS_Status), ]
WES_survival$OS_Status[which(WES_survival$OS_Status == '1:DECEASED')] = 2L
WES_survival$OS_Status[which(WES_survival$OS_Status == '0:LIVING')] = 1L
WES_survival$OS_Status = as.numeric(as.character(WES_survival$OS_Status))
WES_survival$Y_mosaic = factor(WES_survival$Y_mosaic, levels = c('no', 'yes'))


coxph(Surv(OS, OS_Status) ~ Y_mosaic, data = WES_survival) %>%
  gtsummary::tbl_regression(exp = TRUE)

#' survival analysis
f1 = survfit(Surv(OS, OS_Status) ~ Y_mosaic, data = WES_survival)
WES_survival$age = ifelse(WES_survival$AGE_AT_SEQ_REPORTED_YEARS < 50, 'group1', 
                          ifelse(WES_survival$AGE_AT_SEQ_REPORTED_YEARS > 50 & WES_survival$AGE_AT_SEQ_REPORTED_YEARS < 70, 'group2', 'group3')) 
                                 
f2 = survfit(Surv(OS, OS_Status) ~ Y_mosaic + age, data = WES_survival)

#' make the plot
a = ggsurvplot(f1,
           size = 0.8,
           censor.shape = '',
           pval = T,
           risk.table = F,
           risk.table.height = 0.25,
           break.time.by = 12,
           ggtheme = theme_light(),
           surv.median.line = 'v',
           xlim = c(0, 84),
           xlab = ('OS [months]'),
           palette = 'aaas')
  
b = a$plot + theme(panel.background = element_blank()) +
  theme(aspect.ratio = 0.5,
        panel.border = element_rect(size = 1.5, colour = 'black'),
        axis.text = element_text(size = 14, color = 'black', face = 'bold'),
        legend.position = 'top', 
        legend.text = element_text(size = 14)) +
  scale_color_discrete(name = "", labels = c("full Y-chromosome", "mosaic Y loss"), 
                       type = c('#414A8F', '#CD2703'))

ggsave_golden(filename = 'Figures/OS_Mosaic_YLoss.pdf', plot = b, width = 9)




###############################################################################
#' Survival analysis for IMPACT cases
#' modify input data first
clinical = read.csv('Data_out/IMPACT/IMPACT_clinical.txt', sep = '\t')
IMPACT_mosaic = read.csv('Data_out/IMPACT/GermlineCN.txt', sep = '\t')
IMPACT_mosaic = IMPACT_mosaic[which(IMPACT_mosaic$target == 'Y'), ]

OS_Impact = merge(IMPACT_mosaic, clinical[, c('SAMPLE_ID', 'FGA', 'OS_Months', 'OS_Status')], 
                  by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)


#' OS analysis
library(survminer)
library(survival)

OS_Impact = OS_Impact[!is.na(OS_Impact$Y_mosaicism), ]
OS_Impact$Y_mosaicism = factor(OS_Impact$Y_mosaicism, levels = c('no', 'yes'))

#' cox-p model
coxph(Surv(OS_Months, OS_Status) ~ Y_mosaicism, data = OS_Impact) %>%
  gtsummary::tbl_regression(exp = TRUE)

model_impact = survfit(Surv(OS_Months, OS_Status) ~ Y_mosaicism, data = OS_Impact)

#' make the plot
a = ggsurvplot(model_impact,
               size = 0.8,
               censor.shape = '',
               pval = T,
               risk.table = F,
               risk.table.height = 0.25,
               break.time.by = 12,
               ggtheme = theme_light(),
               surv.median.line = 'v',
               xlim = c(0, 84),
               xlab = ('OS [months]'),
               palette = 'aaas')

b = a$plot + theme(panel.background = element_blank()) +
  theme(aspect.ratio = 0.5,
        panel.border = element_rect(size = 1.5, colour = 'black'),
        axis.text = element_text(size = 14, color = 'black', face = 'bold'),
        legend.position = 'top', 
        legend.text = element_text(size = 14)) +
  scale_y_continuous(expand = c(0.01, 0.05)) +
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_color_discrete(name = "", labels = c("full Y-chromosome", "mosaic Y loss"), 
                       type = c('#414A8F', '#CD2703'))

ggsave_golden(filename = 'Figures/OS_Mosaic_IMPACT.pdf', plot = b, width = 9)




















