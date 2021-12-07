## Run a logistic regression model on the whole Impact cohort;
## Look into which variables influence the Y-chromosome loss probability
## 
## start: 11/2021
## revision: 12/2021
## chris-kreitzer


setup(working.path = '~/Documents/GitHub/Y_chromosome_loss/PanCancer/')
clean()
source('Scripts/Plotting_theme.R')


## Libraries and input data
library(data.table)
library(car)
# install.packages('aod')
library(aod)
library(rcompanion)
library(ggplot2)
library(cowplot)


clinical = read.csv('Data_out/IMPACT/IMPACT_clinical.txt', sep = '\t')
clinical$Y_call[which(clinical$Y_call == 'intact_Y_chrom')] = 0
clinical$Y_call[which(clinical$Y_call == 'Y_chrom_loss')] = 1
clinical$Y_call = as.integer(as.character(clinical$Y_call))
clinical_glm = clinical[, c('CANCER_TYPE', 'AGE_AT_SEQ_REPORTED_YEARS', 'SAMPLE_TYPE',
                            'Y_call', 'ploidy', 'purity', 'FGA', 'TMB', 'MSI_Score', 'MSI_Type')]

#' only work with cancer types whose frequency is > 1% in the whole cohort (n = 12,405)
#' also modify the data frame columns, so that they are properly encoded
CancerTypes_considered = names(which(table(clinical_glm$CANCER_TYPE) / sum(table(clinical_glm$CANCER_TYPE)) > 0.01))
clinical_glm = clinical_glm[which(clinical_glm$CANCER_TYPE %in% CancerTypes_considered), ]

clinical_glm$CANCER_TYPE = as.factor(as.character(clinical_glm$CANCER_TYPE))
clinical_glm$AGE_AT_SEQ_REPORTED_YEARS = as.integer(as.character(clinical_glm$AGE_AT_SEQ_REPORTED_YEARS))
clinical_glm$SAMPLE_TYPE = as.factor(as.character(clinical_glm$SAMPLE_TYPE))
clinical_glm$MSI_Type = as.factor(as.character(clinical_glm$MSI_Type))
clinical_glm$MSI_Type = factor(clinical_glm$MSI_Type, levels = c('Stable', 'Instable', 'Indeterminate', 'Do not report'))


#' if we are modelling a glm model without any predictor variable; just an intercept
model_intercept = glm(Y_call ~ 1, data = clinical_glm, family = binomial(link = 'logit'))
summary(model_intercept)
car::vif(model_intercept) #' not possible with only one variable;


#' including one predictor variable: FGA
model_1predictor = glm(Y_call ~ FGA, data = clinical_glm, family = binomial(link = 'logit'))
summary(model_1predictor)
car::vif(model_1predictor) #' also not possible due to only one predictor

#log(p/(1-p)) = logit(p) = -0.93740 + 1.05054*FGA
f = -0.93740 + 1.05054*1
s = -0.93740 + 1.05054*0.1
t = -0.93740 + 1.05054*0.2
f = -0.93740 + 1.05054*0.3
exp(-0.9374) / (1+exp(-0.9374))


#' full glm approach (without CancerType included)
glm_model = clinical_glm[, -which(names(clinical_glm) == "CANCER_TYPE")]
model1 = glm(Y_call ~., data = glm_model, family = binomial(link = 'logit'))
summary(model1)


#' holistic glm approach (including CancerType and all other variables)
model_full = glm(Y_call ~ ., data = clinical_glm, family = binomial(link = 'logit'))
summary(model_full)

#' what about the overall effect of CANCER_TYPE?
#' we are using a wald-test to assess the significance
wald.test(b = coef(model_full), Sigma = vcov(model_full), Terms = 2:20)

#' comparing two models with each other
rcompanion::compareGLM(model_full, model_intercept)


#' check for multicollinearity with car::vif()
model.vif = as.data.frame(car::vif(model_full))
model.vif$Test = model.vif$`GVIF^(1/(2*Df))` ^ 2

#' we will exclude the MSI_Score predictor variable as it shows a high VIF
data_final = clinical_glm[,-which(names(clinical_glm) == 'MSI_Score')]
model_final = glm(Y_call ~., data = data_final, family = binomial(link = 'logit'))
car::vif(model_final)

#' log-likelihood ratio test whether model_final predicts better than null model
A = logLik(model_intercept)
B = logLik(model_final)

LLR = -2*(as.numeric(A) - as.numeric(B))
pchisq(LLR, 28, lower.tail = FALSE)


# Extract model variables
model.vars = as.data.table(summary(model1)$coefficients, keep.rownames = T)

# Compute 95% confidence intervals for odds ratios
model.vars[, OR := exp(Estimate)]
model.vars[, OR_lower := exp(Estimate - 1.96 * `Std. Error`)]
model.vars[, OR_upper := exp(Estimate + 1.96 * `Std. Error`)]
model.vars = model.vars[order(`Pr(>|z|)`)]

model.vars = model.vars[model.vars$`Pr(>|z|)` < 0.001, ]
model.vars$adjusted = p.adjust(p = model.vars$`Pr(>|z|)`, method = 'BH')


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




#' looking into MSI_Type; it seems like that MSI instable tumors show a tendency for retaining the 
#' Y chromosome; 

MSI_out = data.frame()
for(i in unique(clinical_glm$CANCER_TYPE)){
  try({
    da = clinical_glm[which(clinical_glm$CANCER_TYPE == i), ]
    Y.stable = table(da$Y_call[which(da$MSI_Type == 'Stable')])[[2]] / sum(table(da$Y_call[which(da$MSI_Type == 'Stable')]))
    Y.instable = table(da$Y_call[which(da$MSI_Type == 'Instable')])[[2]] / sum(table(da$Y_call[which(da$MSI_Type == 'Instable')]))
    show = ifelse(Y.stable > Y.instable, 'yes', 'no')
    out = data.frame(CancerType = i,
                     Type = c('Stable', 'Instable'),
                     prop = c(Y.stable, Y.instable),
                     n = c(nrow(da[which(da$MSI_Type == 'Stable'), ]),
                           nrow(da[which(da$MSI_Type == 'Instable'), ])),
                     show = rep(show, 2))
    MSI_out = rbind(MSI_out, out)
  })
}


#' Visualization
MSI.plot = ggplot(MSI_out, aes(x = Type, y = prop, color = CancerType, size = n, label = CancerType)) +
  geom_text(aes(label = CancerType), hjust = 0) +
  geom_hline(yintercept = seq(0.2, 0.6, 0.2), color = 'grey85', size = 0.2) +
  geom_line(aes(group = CancerType), 
            color = ifelse(MSI_out$show == 'yes', 'grey85', 'red'), size = 0.35) +
  
  geom_point() +
  
  scale_colour_viridis_d(direction = -1) +
  scale_size_area(limits = c(1, 1500),
                  breaks = c(100, 250, 500, 750, 1000, 1500)) +
  theme_Y(base_size = 14) +
  labs(x = 'MSI Type', y = 'Y-chromosome loss')
  
ggsave_golden(filename = 'Figures/MSI_Type_Y.loss.pdf', plot = MSI.plot, width = 11)


## MSI_TYPE unstable show generally decreased FGA and hence, Y-chromosome loss can be
## explained by FGA again; 

#' A) Colorectal Cancer
Colorectal = clinical_glm[which(clinical_glm$CANCER_TYPE == 'Colorectal Cancer'), ]
Colorectal = droplevels(Colorectal[which(Colorectal$MSI_Type %in% c('Stable', 'Instable')),])


boxplot(Colorectal$FGA ~ Colorectal$MSI_Type,
        main = 'Colorectal Cancer',
        xlab = 'MSI-Type',
        ylab = 'FGA',
        notch = T,
        asp = 10)

t.test(Colorectal$FGA ~ Colorectal$MSI_Type)


#' B) Esophagogastric Cancer
Esophagus = clinical_glm[which(clinical_glm$CANCER_TYPE == 'Esophagogastric Cancer'), ]
Esophagus = droplevels(Esophagus[which(Esophagus$MSI_Type %in% c('Stable', 'Instable')),])
boxplot(Esophagus$FGA ~ Esophagus$MSI_Type,
        main = 'Esophagogastric Cancer',
        xlab = 'MSI-Type',
        ylab = 'FGA',
        notch = T,
        asp = 10)

t.test(Esophagus$FGA ~ Esophagus$MSI_Type)

#' C) Pancreatic Cancer:
Pancreatic = clinical_glm[which(clinical_glm$CANCER_TYPE == 'Pancreatic Cancer'), ]
Pancreatic = droplevels(Pancreatic[which(Pancreatic$MSI_Type %in% c('Stable', 'Instable')),])
boxplot(Pancreatic$FGA ~ Pancreatic$MSI_Type,
        main = 'Esophagogastric Cancer',
        xlab = 'MSI-Type',
        ylab = 'FGA',
        notch = T,
        asp = 10)

t.test(Pancreatic$FGA ~ Pancreatic$MSI_Type)



#' D) Prostate Cancer:
Prostate = clinical_glm[which(clinical_glm$CANCER_TYPE == 'Prostate Cancer'), ]
Prostate = droplevels(Prostate[which(Prostate$MSI_Type %in% c('Stable', 'Instable')),])
boxplot(Prostate$FGA ~ Prostate$MSI_Type,
        main = 'Esophagogastric Cancer',
        xlab = 'MSI-Type',
        ylab = 'FGA',
        notch = T,
        asp = 10)

t.test(Prostate$FGA ~ Prostate$MSI_Type)
dim(Prostate)


###############################################################################
#' The lower the purity, the less likely we observe a Y-chromosome loss call. 
library(dplyr)

purity_out = data.frame()
for(i in seq(0, 0.9, by = 0.1)){
  da = clinical_glm[between(x = clinical_glm$purity, left = i, right = i+0.1), ]
  frac = (table(da$Y_call)[[2]] / sum(table(da$Y_call))) *100
  out = data.frame(group = i,
                   n = nrow(da),
                   frac = frac)
  purity_out = rbind(purity_out, out)
}

purity_out$group = as.factor(as.numeric(purity_out$group))

#' Visualization
n.purity = ggplot(purity_out, aes(x = group, y = n)) +
  geom_bar(stat = 'identity', width = 0.6) +
  labs(x = '', y = 'cases (n)')
frac.purity = ggplot(purity_out, aes(x = group, y = frac)) +
  geom_bar(stat = 'identity', width = 0.6) +
  labs(x = 'purity-group', y = 'Y-chromosome loss [%]')
library(patchwork)
n.purity / frac.purity



#' relative contribution cancer types
relCancer = data.frame()
for(i in seq(0, 0.9, by = 0.1)){
  da = clinical_glm[between(x = clinical_glm$purity, left = i, right = i+0.1), ]
  freq = data.frame(table(da$CANCER_TYPE))
  for(j in 1:nrow(freq)){
    type = freq$Var1[j]
    proc = freq$Freq[j] / sum(freq$Freq)
    out = data.frame(group = factor(i),
                     CancerType = type,
                     frac = proc)
    relCancer = rbind(relCancer, out)
  }
}
  

purity_Cancer = ggplot(relCancer[order(relCancer$frac, decreasing = T), ], 
       aes(x = group, y = frac, fill = CancerType, label = CancerType)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_viridis_d(begin = 0, option = 'B') +
  geom_text(size = 2, position = position_stack(vjust = 0.5), color = 'white') +
  labs(x = 'purity-group', y = 'Fraction') +
  theme(legend.position = 'none')


sum = plot_grid(n.purity, frac.purity, purity_Cancer, rel_heights = c(1,1,3), nrow = 3, align = 'h')
ggsave_golden(filename = 'Figures/PurityCancer.pdf', plot = sum, width = 12)
