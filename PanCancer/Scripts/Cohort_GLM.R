##----------------+
## Run GLM model on MSK-IMPACT cohort
## Which variables influence the LOY 
##----------------+
## start: 11/2021
## revision: 12/2021
## revision: 11/18/2022
## revision: 11/21/2022
## revision: 01/13/2023
## revision: 02/03/23
## 
## chris-kreitzer


clean()
gc()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')


## Libraries and input data
library(data.table)
library(car)
library(aod)
library(rcompanion)
library(ggplot2)
library(cowplot)


##-------
## TODO:
## - QC TRUE and FALSE samples
## - what to do with Gain_Loss segments?


cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]

cohort$classification[which(cohort$classification %in% c('complete_loss', 'partial_loss', 'relative_loss', 'gain_loss'))] = 1
cohort$classification[which(cohort$classification %in% c('gain', 'wt', 'partial_gain'))] = 0
cohort = cohort[!is.na(cohort$classification), ]


# cohort$classification[which(cohort$sample == 'P-0002130-T02-IM3')] = 0
# cohort$classification[which(cohort$sample == 'P-0013100-T01-IM5')] = 0
# cohort$classification[which(cohort$sample == 'P-0028409-T01-IM6')] = 1
# cohort$classification[which(cohort$sample == 'P-0064933-T02-IM7')] = 0

clinical_glm = cohort[,c('SAMPLE_ID', 'classification','ploidy', 'purity', 'CANCER_TYPE', 'MSI_TYPE',
                         'Age_Sequencing', 'SAMPLE_TYPE', 'genome_doubled', 'fraction_cna')]


#' only work with cancer types whose frequency is > 1% in the whole cohort (n = 12,405)
#' also modify the data frame columns, so that they are properly encoded
CancerTypes_considered = names(which(table(clinical_glm$CANCER_TYPE) / sum(table(clinical_glm$CANCER_TYPE)) > 0.01))
clinical_glm = clinical_glm[which(clinical_glm$CANCER_TYPE %in% CancerTypes_considered), ]
clinical_glm$CANCER_TYPE = as.factor(as.character(clinical_glm$CANCER_TYPE))
clinical_glm$Age_Sequencing = as.numeric(as.character(clinical_glm$Age_Sequencing))
clinical_glm$SAMPLE_TYPE = as.factor(as.character(clinical_glm$SAMPLE_TYPE))
clinical_glm$MSI_TYPE = as.factor(as.character(clinical_glm$MSI_TYPE))
clinical_glm$MSI_TYPE = factor(clinical_glm$MSI_TYPE, levels = c('Stable', 'Instable', 'Indeterminate', 'Do not report'))
clinical_glm$genome_doubled = as.integer(clinical_glm$genome_doubled)


##----------------+
## full glm approach 
## without CancerType included
##----------------+
glm_model = clinical_glm[, -which(names(clinical_glm) == "CANCER_TYPE")]
glm_model = glm_model[!is.na(glm_model$purity), ]
glm_model$classification = as.integer(as.character(glm_model$classification))
model_full = glm(classification ~ ploidy+purity + MSI_TYPE + Age_Sequencing + SAMPLE_TYPE + genome_doubled + fraction_cna, data = glm_model, family = binomial)
summary(model_full)

#' check for multicollinearity with car::vif()
model.vif = as.data.frame(car::vif(model_full))
model.vif$Test = model.vif$`GVIF^(1/(2*Df))` ^ 2

##----------------+
## excluding variables
## that show VIF/TEST > 2.5
##----------------+
# glm_model = glm_model[,-which(names(glm_model) == 'TMB')]
# glm_model = glm_model[,-which(names(glm_model) == 'WGD')]
# glm_model = glm_model[,-which(names(glm_model) == 'Mut_Count')]
# glm_model = glm_model[,-which(names(glm_model) == 'MSI_SCORE')]
glm_model = glm_model[,-which(names(glm_model) == 'SAMPLE_ID')]
glm_model$genome_doubled = NULL

model_final = glm(classification ~., data = glm_model, family = binomial(link = 'logit'))
car::vif(model_final)
summary(model_final)


##----------------+
## Wald-test
##----------------+
library(aod)
wald.test(Sigma = vcov(model_final), b = coef(model_final), Terms = 2)


##----------------+
## Extract model variables
## compute 95% CI
##----------------+
model.vars = as.data.table(summary(model_final)$coefficients, keep.rownames = T)
model.vars[, OR := exp(Estimate)]
model.vars[, OR_lower := exp(Estimate - 1.96 * `Std. Error`)]
model.vars[, OR_upper := exp(Estimate + 1.96 * `Std. Error`)]
model.vars = model.vars[order(`Pr(>|z|)`)]
model.vars$adjusted = p.adjust(p = model.vars$`Pr(>|z|)`, method = 'fdr')
model.vars = model.vars[model.vars$adjusted < 0.05, ]

#' exclude intercept and MSI type in plotting
model.vars = model.vars[which(model.vars$rn != '(Intercept)'), ]
model.vars = model.vars[which(model.vars$rn != 'MSI_TYPEIndeterminate'), ]
model.vars = model.vars[which(model.vars$rn != 'MSI_TYPEDo not report'), ]
model.vars$rn = c('FGA', 'Purity', 'Ploidy', 'Age', 'MSI')

#' Visualization
cohort_glm_plot = ggplot(data = model.vars, 
       aes(x = reorder(rn, Estimate), 
           y = Estimate)) +
  geom_pointrange(aes(ymin = Estimate - `Std. Error`,
                      ymax = Estimate + `Std. Error`),
                  size = 1, fatten = 2, linewidth = 1) +
  coord_flip() +
  theme_std(base_size = 14) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        aspect.ratio = 1,
        axis.title.x.top = element_blank(),
        axis.line = element_blank()) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_y_continuous(sec.axis = dup_axis(),
                     limits = c(-2,3.5)) +
  labs(x = '', y = 'log(odds ratio)') +
  panel_border(size = 2, color = 'black')

write.table(x = model.vars, file = 'Data/05_Mutation/Cohort_GLM.txt', sep = '\t', row.names = F)
ggsave_golden(filename = 'Figures_original/Cohort_GLM_all.pdf', plot = cohort_glm_plot, width = 6)





#' out