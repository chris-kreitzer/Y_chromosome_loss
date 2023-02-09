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
## revision: 02/08/2023
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
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')

library(readxl)
library(data.table)
library(survival)
library(survminer)
library(dplyr)



##-------
## Data rendering;
##-------
Cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
Cohort = Cohort[which(Cohort$Study_include == 'yes'), ]
OS_cohort = Cohort[,c('PATIENT_ID', 'SAMPLE_ID', 'CANCER_TYPE', 'QC', 'Age_Sequencing', 'genome_doubled', 'SAMPLE_TYPE', 'MSI_TYPE',
                      'classification', 'OS_Status', 'OS_months')]
OS_cohort = OS_cohort[!is.na(OS_cohort$OS_Status) & !is.na(OS_cohort$OS_months), ]
colnames(OS_cohort)[9] = 'LOY'
OS_cohort$LOY = ifelse(OS_cohort$LOY == 'complete_loss', TRUE, FALSE)
OS_cohort$OS_Status = ifelse(OS_cohort$OS_Status == '1:DECEASED', 'Dead', 'Alive')
OS_cohort$OS_Status_INT = ifelse(OS_cohort$OS_Status == 'Dead', 2L, 1L)
data.gam = readRDS(file = 'Data/05_Mutation/data_gam.rds')
data_gam = data.gam$GAM_Analysis

OS_cohort = merge(OS_cohort, data_gam[, c('sample', 'TP53_mut', 'TP53_Deletion')],
                  by.x = 'SAMPLE_ID', by.y = 'sample', all.x = T)
OS_cohort$TP53_Deletion = as.numeric(as.integer(OS_cohort$TP53_Deletion))

OS_cohort$TP53_Status = ifelse(OS_cohort$TP53_Deletion == 1, 1,
                               ifelse(OS_cohort$TP53_mut == 1, 1, 0))
OS_cohort$TP53_Deletion = NULL
OS_cohort$TP53_mut = NULL
OS_cohort$TP53_Status = as.integer(as.numeric(OS_cohort$TP53_Status))


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




##-------
## COX-proportional hazard model;
## PanCancer: 
## - accounting for TP53 mutational status
## - and cancer type
##-------
pancancer_TP53_ct = coxph(formula = Surv(OS_months, as.numeric(factor(OS_Status))) ~ LOY+TP53_Status+CANCER_TYPE, data = OS_cohort)
summary(pancancer_TP53_ct)
coxph(Surv(OS_months, OS_Status_INT) ~ LOY+TP53_Status+CANCER_TYPE, data = OS_cohort) %>%
  gtsummary::tbl_regression(exp = TRUE)



##-------
## Kaplan-Meier method
##-------
pancancerTP53 = survfit(Surv(OS_months, OS_Status_INT) ~ LOY+TP53_Status, data = OS_cohort)
p1 = ggsurvplot(pancancerTP53,
                palette = c('lightblue', 'pink', "#3d4397", "#a22231"),
                size = 0.8,
                legend = 'top',
                xlab = 'Time (months)',
                ylab = 'Overall survival (%)',pval = T)$plot +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(aspect.ratio = 1)
p1

ggsave_golden(filename = 'Figures_original/OS_wholeCohort_TP53mut.pdf', plot = p1, width = 8)




##----------------+
## Look into individual cancer types
##----------------+
OS_cohort = OS_cohort[which(OS_cohort$CANCER_TYPE %in% ctypes_keep), ]
OS_cancertypes = data.frame()
for(i in unique(OS_cohort$CANCER_TYPE)){
  data_sub = OS_cohort[which(OS_cohort$CANCER_TYPE == i), ]
  cancer = i
  n = nrow(data_sub)
  if(n < 60) next
  else {
    model_raw = coxph(formula = Surv(OS_months, OS_Status_INT) ~ LOY, data = data_sub)
    model_adj = coxph(formula = Surv(OS_months, OS_Status_INT) ~ LOY+TP53_Status, data = data_sub)
    estimate_raw = summary(model_raw)$coefficients[[1]]
    estimate_se_raw = summary(model_raw)$coefficients[[3]]
    significance_raw = summary(model_raw)$coefficients[[5]]
    exp_estimate_raw = summary(model_raw)$conf.int[[1]]
    
    
    ##-- adjusted
    estimate_adj = summary(model_adj)$coefficients[[1]]
    estimate_se_adj = summary(model_adj)$coefficients[[5]]
    significance_adj = summary(model_adj)$coefficients[[9]]
    exp_estimate_adj = summary(model_adj)$conf.int[[1]]
    out = data.frame(cancer = cancer,
                     n = n,
                     estimate = estimate_raw,
                     estimate_se = estimate_se_raw,
                     exp_estimate = exp_estimate_raw,
                     significance = significance_raw,
                     estimate_adj = estimate_adj,
                     estimate_se_adj = estimate_se_adj,
                     exp_estimate_adj = exp_estimate_adj,
                     significance_adj = significance_adj)
    OS_cancertypes = rbind(OS_cancertypes, out)
  }
}

OS_cancertypes
OS_cancertypes$cancer = paste0(OS_cancertypes$cancer, ' (n=', OS_cancertypes$n, ')')
OS_cancertypes$plot_color = ifelse(OS_cancertypes$significance <= 0.05, 'plot', 'notplot')
OS_cancertypes$plot_adjusted = ifelse(OS_cancertypes$significance_adj <= 0.05, 'plot', 'notplot')



##-------
## Visualization
##-------
OS_CancerTypes_plot = ggplot(data = OS_cancertypes, 
                             aes(x = reorder(cancer, estimate_adj), 
                                 y = estimate_adj)) +
  geom_point(aes(colour = plot_adjusted), size = 1.5) +
  geom_pointrange(aes(ymin = estimate_adj - estimate_se_adj,
                      ymax = estimate_adj + estimate_se_adj,
                      colour = plot_adjusted),
                  size = 0.5) +
  coord_flip() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        aspect.ratio = 1,
        axis.title.x.top = element_blank(),
        axis.text = element_text(size = 10, face = 'bold', color = 'black')) +
  scale_color_manual(values = c('notplot' = 'grey',
                                'plot' = 'red'),
                     labels = c('not_significant', 'significant'),
                     name = '') +
  geom_hline(yintercept = 0) +
  scale_y_continuous(sec.axis = dup_axis()) +
  labs(x = '', y = 'Model coefficient (log [HR])')

OS_CancerTypes_plot


ggsave_golden(filename = 'Figures_original/OS_CancerTypes_adj.pdf', plot = OS_CancerTypes_plot, width = 8)










##----------------+
## Adjustment of CoxP model
## for age and tumor site;
## 
## we can potentially include
## more variables (later)
##----------------+
loy.pancancer_MV = coxph(formula = Surv(OS_months, as.numeric(factor(OS_Status))) ~ LOY+genome_doubled+SAMPLE_TYPE, data = OS_cohort)

xx = survfit(Surv(OS_months, OS_Status_INT) ~ LOY+genome_doubled, data = OS_cohort)
ggsurvplot(fit = survfit(loy.pancancer_MV), data = OS_cohort)
ggsurvplot(fit = xx)

##---- Craig Example: WGD
xx = read.csv('~/Desktop/breast.txt', sep = '\t')
model.breast <- coxph(
  formula = Surv(as.numeric(os_days) / 30,
                 as.numeric(factor(VitalStatus))) ~ gd + DetailedTumorType2 + age_at_diagnosis + menopause_status + stage + grade + esr1_mutant,
  data = xx
)
summary(model.breast)

str(xx)
head(xx)
##----------------+
## run CoxPh model on per
## cancer-type level
## additionally adjust for
## co-variates;
##----------------+

OS_cancertypes = data.frame()
for(i in unique(OS_cohort$CANCER_TYPE)){
  data_sub = OS_cohort[which(OS_cohort$CANCER_TYPE == i), ]
  cancer = i
  n = nrow(data_sub)
  if(n < 60) next
  else {
    model = coxph(formula = Surv(OS_months, OS_Status_INT) ~ LOY, data = data_sub)
    estimate = summary(model)$coefficients[[1]]
    estimate_se = summary(model)$coefficients[[3]]
    significance = summary(model)$coefficients[[5]]
    exp_estimate = summary(model)$conf.int[[1]]
    lower = summary(model)$conf.int[[3]]
    upper = summary(model)$conf.int[[4]]
    out = data.frame(cancer = cancer,
                     n = n,
                     estimate = estimate,
                     estimate_se = estimate_se,
                     exp_estimate = exp_estimate,
                     lower = lower,
                     upper = upper,
                     significance = significance)
    OS_cancertypes = rbind(OS_cancertypes, out)
    rm(data_sub, model, estimate, significance)
  }
}


OS_cancertypes
OS_cancertypes$cancer = paste0(OS_cancertypes$cancer, ' (n=', OS_cancertypes$n, ')')
OS_cancertypes$plot_color = ifelse(OS_cancertypes$significance <= 0.01, 'plot', 'notplot')


##-------
## Visualization
##-------
OS_CancerTypes_plot = ggplot(data = OS_cancertypes, aes(x = reorder(cancer, estimate), y = estimate)) +
  geom_point(aes(colour = plot_color), size = 1.5) +
  geom_pointrange(aes(ymin = estimate - estimate_se,
                      ymax = estimate + estimate_se,
                      colour = plot_color),
                  size = 0.5) +
  coord_flip() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        aspect.ratio = 1,
        axis.title.x.top = element_blank(),
        axis.text = element_text(size = 10, face = 'bold', color = 'black')) +
  scale_color_manual(values = c('notplot' = 'grey',
                                'plot' = 'red'),
                     labels = c('not_significant', 'significant'),
                     name = '') +
  geom_hline(yintercept = 0) +
  scale_y_continuous(sec.axis = dup_axis()) +
  labs(x = '', y = 'Model coefficient (log [HR])')

ggsave_golden(filename = 'Figures_original/OS_CancerTypes.pdf', plot = OS_CancerTypes_plot, width = 12)




