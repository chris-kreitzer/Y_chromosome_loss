rm(list = ls())
.rs.restartR()
setwd('~/Documents/GitHub/Y_chromosome_loss/PanCancer/')

library(data.table)
library(car)

clinical = read.csv('Data_out/IMPACT/IMPACT_clinical.txt', sep = '\t')
clinical$Y_call[which(clinical$Y_call == 'intact_Y_chrom')] = 0
clinical$Y_call[which(clinical$Y_call == 'Y_chrom_loss')] = 1
clinical$Y_call = as.integer(as.character(clinical$Y_call))
clinical_glm = clinical[, c('CANCER_TYPE', 'AGE_AT_SEQ_REPORTED_YEARS', 'SAMPLE_TYPE',
                            'Y_call', 'ploidy', 'purity', 'FGA', 'TMB', 'MSI_Score', 'MSI_Type')]


#' if we are modelling a glm model without any predictor variable; just an intercept
summary(glm(Y_call ~ 1, data = clinical_glm, family = binomial(link = 'logit')))

#' including one predictor variable: FGA
summary(glm(Y_call ~ FGA, data = clinical_glm, family = binomial(link = 'logit')))
#log(p/(1-p)) = logit(p) = -0.93740 + 1.05054*FGA
-0.93740 + 1.05054*0
f = -0.93740 + 1.05054*1
s = -0.93740 + 1.05054*0.2
t = -0.93740 + 1.05054*0.3
exp(t - s)
exp(-0.9374) / (1+exp(-0.9374))


model1 = glm(Y_call ~ ., data = clinical_glm, family = binomial(link = 'logit'))
model.vif = as.data.table(car::vif(model1))
model.vif[, Test := (`GVIF^(1/(2*Df))`) ^ 2]
any(model.vif$Test > 10) # FALSE

# Extract model variables
model.vars = as.data.table(summary(model1)$coefficients, keep.rownames = T)

# Compute 95% confidence intervals for odds ratios
model.vars[, OR := exp(Estimate)]
model.vars[, OR_lower := exp(Estimate - 1.96 * `Std. Error`)]
model.vars[, OR_upper := exp(Estimate + 1.96 * `Std. Error`)]
model.vars = model.vars[order(`Pr(>|z|)`)]

model.vars = model.vars[model.vars$`Pr(>|z|)` < 0.001, ]
model.vars$adjusted = p.adjust(p = model.vars$`Pr(>|z|)`, method = 'BH')

## print:
# model.print$n = NA
# for(i in unique(keep2)){
#   print(i)
#   n = table(data.out[, i])[[2]][[1]]
#   print(n)
#   model.print$n[which(model.print$rn == i)] = n
# }


#' Visualization
ggplot(data = model.vars, aes(x = reorder(rn, Estimate), y = Estimate)) +
  geom_point(size = 2) +
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


