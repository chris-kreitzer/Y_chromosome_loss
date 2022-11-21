##----------------+
## Run a logistic regression model 
## on the whole Impact (male) cohort;
## Look into which variables influence 
## the Y-chromosome loss probability
##----------------+
## start: 11/2021
## revision: 12/2021
## revision: 11/18/2022
## revision: 11/21/2022
## chris-kreitzer


clean()
gc()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
source('Scripts/plot_theme.R')


## Libraries and input data
library(data.table)
library(car)
library(aod)
library(rcompanion)
library(ggplot2)
library(cowplot)


cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
clinical = cohort$IMPACT_clinicalAnnotation
loy = cohort$IMPACT_Y_classification_final
loy_clinical = cohort$IMPACT_binaryY_call
arm_level = cohort$IMPACT_ARM_level_changes

clinical_glm = merge(loy[,c('sample', 'ploidy', 'classification')], loy_clinical[,c('sample.id', 'purity')],
                     by.x = 'sample', by.y = 'sample.id', all.x = T)

clinical_glm$classification[which(clinical_glm$classification %in% c('loss', 'relative_loss'))] = 1
clinical_glm$classification[which(clinical_glm$classification %in% c('gain', 'wt', 'gain_loss'))] = 0
clinical_glm$classification[which(clinical_glm$sample == 'P-0002130-T02-IM3')] = 0
clinical_glm$classification[which(clinical_glm$sample == 'P-0013100-T01-IM5')] = 0
clinical_glm$classification[which(clinical_glm$sample == 'P-0028409-T01-IM6')] = 1
clinical_glm$classification[which(clinical_glm$sample == 'P-0064933-T02-IM7')] = 0

clinical_glm = merge(clinical_glm, clinical[,c('SAMPLE_ID', 'CANCER_TYPE', 'MSI_TYPE', 
                                               'MSI_SCORE', 'AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)', 
                                               'IMPACT_TMB_SCORE', 'MUTATION_COUNT')],
                     by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
clinical_glm = merge(clinical_glm, cohort$IMPACT_cohort[,c('SAMPLE_ID', 'SAMPLE_TYPE')],
                     by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
clinical_glm = merge(clinical_glm, arm_level[,c('id', 'genome_doubled', 'fraction_cna')],
                     by.x = 'sample', by.y = 'id', all.x = T)
clinical_glm$classification = as.integer(as.character(clinical_glm$classification))
colnames(clinical_glm) = c('sample', 'ploidy', 'Y_call', 'purity', 'CANCER_TYPE', 'MSI_TYPE', 
                           'MSI_SCORE', 'Age', 'TMB', 'Mut_Count', 'SAMPLE_TYPE', 'WGD', 'FGA')

#' only work with cancer types whose frequency is > 1% in the whole cohort (n = 12,405)
#' also modify the data frame columns, so that they are properly encoded
CancerTypes_considered = names(which(table(clinical_glm$CANCER_TYPE) / sum(table(clinical_glm$CANCER_TYPE)) > 0.01))
clinical_glm = clinical_glm[which(clinical_glm$CANCER_TYPE %in% CancerTypes_considered), ]
clinical_glm$CANCER_TYPE = as.factor(as.character(clinical_glm$CANCER_TYPE))
clinical_glm$Age = as.integer(as.character(clinical_glm$Age))
clinical_glm$SAMPLE_TYPE = as.factor(as.character(clinical_glm$SAMPLE_TYPE))
clinical_glm$MSI_TYPE = as.factor(as.character(clinical_glm$MSI_TYPE))
clinical_glm$MSI_TYPE = factor(clinical_glm$MSI_TYPE, levels = c('Stable', 'Instable', 'Indeterminate', 'Do not report'))
clinical_glm$WGD[which(clinical_glm$WGD == TRUE)] = 1
clinical_glm$WGD[which(clinical_glm$WGD == FALSE)] = 0
clinical_glm$WGD = as.factor(clinical_glm$WGD)


##----------------+
## full glm approach 
## without CancerType included
##----------------+
glm_model = clinical_glm[, -which(names(clinical_glm) == "CANCER_TYPE")]
model_full = glm(Y_call ~ ploidy+purity+MSI_TYPE+Age+TMB+SAMPLE_TYPE+WGD+FGA, data = glm_model, family = binomial)
summary(model_full)

#' check for multicollinearity with car::vif()
model.vif = as.data.frame(car::vif(model_full))
model.vif$Test = model.vif$`GVIF^(1/(2*Df))` ^ 2

##----------------+
## excluding variables:
## TMB, 
## Mut_Count, 
## WGD as they show high VIF
##----------------+
glm_model = glm_model[,-which(names(glm_model) == 'TMB')]
glm_model = glm_model[,-which(names(glm_model) == 'WGD')]
glm_model = glm_model[,-which(names(glm_model) == 'Mut_Count')]
glm_model = glm_model[,-which(names(glm_model) == 'MSI_SCORE')]
glm_model = glm_model[,-which(names(glm_model) == 'sample')]

model_final = glm(Y_call ~., data = glm_model, family = binomial(link = 'logit'))
car::vif(model_final)
# any model_final$test > 10 - variance inflation factor
summary(model_final)


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

#' Visualization
ggplot(data = model.vars, aes(x = reorder(rn, Estimate), y = Estimate)) +
  geom_point(size = 1.5) +
  geom_pointrange(aes(ymin = Estimate - `Std. Error`,
                      ymax = Estimate + `Std. Error`),
                  size = 0.5) +
  coord_flip() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        aspect.ratio = 1,
        axis.title.x.top = element_blank(),
        axis.text = element_text(size = 10, face = 'bold', color = 'black')) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(sec.axis = dup_axis()) +
  labs(x = '', y = 'log(odds ratio)')




##----------------+
## Investigate the purity
## category; what there
##----------------+

## Hypothesis:
#' The lower the purity, 
#' the less likely we observe 
#' a Y-chromosome loss call. 
library(dplyr)

purity_out = data.frame()
for(i in seq(0, 0.9, by = 0.1)){
  da = clinical_glm[dplyr::between(x = clinical_glm$purity, left = i, right = i+0.1), ]
  frac = (table(da$Y_call)[[2]] / sum(table(da$Y_call))) *100
  out = data.frame(group = i,
                   n = nrow(da),
                   frac = frac,
                   average_age = mean(da$Age, na.rm = T))
  purity_out = rbind(purity_out, out)
}

purity_out$group = as.factor(as.numeric(purity_out$group))

#' Visualization
n.purity = ggplot(purity_out, aes(x = group, y = n)) +
  geom_bar(stat = 'identity', width = 0.8) +
  labs(x = '', y = 'cases (n)') + theme_std() +
  theme(axis.text.x = element_blank())
frac.purity = ggplot(purity_out, aes(x = group, y = frac)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_x_discrete(labels = c('[0-0.1]', '[0.1-0.2]', '[0.2-0.3]', '[0.3-0.4]', '[0.4-0.5]',
                              '[0.5-0.6]', '[0.6-0.7]', '[0.7-0.8]', '[0.8-0.9]', '[0.9-1]')) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(x = 'purity-group', y = 'Y-chromosome loss [%]') +
  theme_std() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

n.purity / frac.purity

#' #' relative contribution cancer types
#' relCancer = data.frame()
#' for(i in seq(0, 0.9, by = 0.1)){
#'   da = clinical_glm[dplyr::between(x = clinical_glm$purity, left = i, right = i+0.1), ]
#'   freq = data.frame(table(da$CANCER_TYPE))
#'   for(j in 1:nrow(freq)){
#'     type = freq$Var1[j]
#'     proc = freq$Freq[j] / sum(freq$Freq)
#'     out = data.frame(group = factor(i),
#'                      CancerType = type,
#'                      frac = proc)
#'     relCancer = rbind(relCancer, out)
#'   }
#' }
#' 
# purity_Cancer = ggplot(relCancer[order(relCancer$frac, decreasing = T), ], 
#                        aes(x = group, y = frac, fill = CancerType, label = CancerType)) +
#   geom_bar(stat = 'identity', position = 'stack') +
#   scale_fill_viridis_d(begin = 0, option = 'B') +
#   geom_text(size = 2, position = position_stack(vjust = 0.5), color = 'white') +
#   labs(x = 'purity-group', y = 'Fraction') +
#   theme(legend.position = 'none')
# 
# 
# sum = plot_grid(n.purity, frac.purity, purity_Cancer, rel_heights = c(1,1,3), nrow = 3, align = 'h')


##----------------+
## MSI_Type; 
## it seems like that MSI instable 
## tumors show a tendency for 
## retaining the Y chromosome
##----------------+ 
clinical_glm$MSI_TYPE = as.character(as.factor(clinical_glm$MSI_TYPE))
clinical_glm$CANCER_TYPE = as.character(as.factor(clinical_glm$CANCER_TYPE))
MSI_out = data.frame()
for(i in unique(clinical_glm$CANCER_TYPE)){
  try({
    da = clinical_glm[which(clinical_glm$CANCER_TYPE == i), ]
    Y.stable = table(da$Y_call[which(da$MSI_TYPE == 'Stable')])[[2]] / sum(table(da$Y_call[which(da$MSI_TYPE == 'Stable')]))
    Y.instable = table(da$Y_call[which(da$MSI_TYPE == 'Instable')])[[2]] / sum(table(da$Y_call[which(da$MSI_TYPE == 'Instable')]))
    show = ifelse(Y.stable > Y.instable, 'yes', 'no')
    out = data.frame(CancerType = i,
                     Type = c('Stable', 'Instable'),
                     prop = c(Y.stable, Y.instable),
                     n = c(nrow(da[which(da$MSI_TYPE == 'Stable'), ]),
                           nrow(da[which(da$MSI_TYPE == 'Instable'), ])),
                     show = rep(show, 2))
    MSI_out = rbind(MSI_out, out)
  })
}

MSI_out = MSI_out[MSI_out$CancerType != 'Thyroid Cancer', ]
MSI_out = MSI_out[MSI_out$CancerType != 'Soft Tissue Sarcoma', ]
MSI_out = MSI_out[MSI_out$CancerType != 'Pancreatic Cancer', ]

#' Visualization
MSI.plot = ggplot(MSI_out, aes(x = Type, y = prop, label = CancerType)) +
  geom_text_repel(aes(label = CancerType), hjust = 1) +
  geom_hline(yintercept = seq(0.2, 0.6, 0.2), color = 'grey85', size = 0.2) +
  geom_line(aes(group = CancerType), size = 0.35) +
  geom_point() +
  theme_std(base_size = 14) +
  theme(panel.border = element_rect(fill = NA, size = 2, color = 'black'),
        aspect.ratio = 1) +
  labs(x = 'MSI Type', y = 'Fraction LOY')
  
ggsave_golden(filename = 'Figures_original/MSI_Type_Y.loss.pdf', plot = MSI.plot, width = 6)


## MSI_TYPE unstable show generally 
## decreased FGA and hence, 
## Y-chromosome loss can be explained by FGA again; 

#' A) Colorectal Cancer
Colorectal = clinical_glm[which(clinical_glm$CANCER_TYPE == 'Colorectal Cancer'), ]
Colorectal = Colorectal[which(Colorectal$MSI_TYPE %in% c('Stable', 'Instable')), ]


boxplot(Colorectal$FGA ~ Colorectal$MSI_TYPE,
        main = 'Colorectal Cancer',
        xlab = 'MSI-Type',
        ylab = 'FGA')

t.test(Colorectal$FGA ~ Colorectal$MSI_TYPE)


#' out
