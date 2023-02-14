##----------------+
## Prostate; 
## close-up investigations
## especially OS; with CoxPH
##----------------+
##
## start: 02/13/2023
## revision: 02/14/2023
## chris-kreitzer



clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')


cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
prostate = cohort[which(cohort$CANCER_TYPE == "Prostate Cancer"), ]
prostate = prostate[!is.na(prostate$classification), ]
data_gam = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = data_gam$GAM_Analysis
prad_muts = data_gam[which(data_gam$sample %in% prostate$SAMPLE_ID), ]
copynumber = read.csv('Data/04_Loss/010523/CopyNumberStates.txt', sep = '\t')
copynumber = copynumber[which(copynumber$id %in% prostate$SAMPLE_ID), ]
IGV = read.csv('Data/04_Loss/010523/IGV_out.seg', sep = '\t')

full_prostate = merge(prostate, prad_muts, by.x = 'SAMPLE_ID', by.y = 'sample', all.x = T)

full_prostate$TP53_Deletion

table(prostate$classification)
sum(table(prostate$classification))
dim(prostate)
which(is.na(prostate$classification))
View(prostate)

cbio(prostate$SAMPLE_ID[which(prostate$classification == 'complete_loss')])

pro_out = data.frame()
for(i in c('complete_loss', 'wt')){
  for(j in 42:ncol(full_prostate)){
    print(j)
    alteration = colnames(full_prostate)[j]
    freq = sum(full_prostate[which(full_prostate$classification == i), j], na.rm = T)
    rel_freq = freq / length(full_prostate$SAMPLE_ID[which(full_prostate$classification == i)])
    category = i
    out = data.frame(category = category,
                     alteration = alteration,
                     frequency = rel_freq)
    pro_out = rbind(pro_out, out)
  }
}

pro_summary = data.frame()
for(i in unique(pro_out$alteration)){
  complete_loss = pro_out$frequency[which(pro_out$alteration == i & pro_out$category == 'complete_loss')]
  wt = pro_out$frequency[which(pro_out$alteration == i & pro_out$category == 'wt')]
  alteration = i
  out = data.frame(alteration = alteration,
                   wt = wt,
                   complete_loss = complete_loss)
  pro_summary = rbind(pro_summary, out)
}
pro_summary$display = ifelse(abs(pro_summary$wt - pro_summary$complete_loss) > 0.05, 'plot', 'notplot')

alteration_frequency = ggplot(pro_summary, aes(x = wt, y = complete_loss, color = display, label = alteration)) +
  geom_text_repel(aes(label = ifelse(pro_summary$display == 'plot', pro_summary$alteration, ''))) +
  geom_jitter() +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.45)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 0.45)) +
  scale_color_manual(values = c('plot' = 'darkred',
                                'notplot' = 'grey65',
                                name = '')) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = 'dashed', color = 'grey35') +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1,
        legend.position = 'none') +
  labs(x = 'wild type', y = 'complete loss')

ggsave_golden(filename = 'Figures_original/Prostate_AlterationFreq.pdf', plot = alteration_frequency, width = 6)


##----------------+
## Univariate and multivariate Cox model
##----------------+
gain_loss = full_prostate$SAMPLE_ID[which(full_prostate$classification == 'gain_loss')]
gain_loss_prad = copynumber[which(copynumber$id %in% gain_loss & copynumber$chrom == 24), ]
gain_loss_prad = gain_loss_prad[which(gain_loss_prad$tcn.em == 0), ]


pGL = c('P-0003040-T01-IM5',
'P-0013100-T01-IM5',
'P-0016288-T01-IM6',
'P-0025571-T01-IM6',
'P-0026808-T01-IM6',
'P-0029183-T01-IM6',
'P-0037069-T01-IM6',
'P-0042930-T01-IM6',
'P-0053240-T01-IM6',
'P-0060161-T01-IM7',
'P-0067186-T01-IM7',
'P-0075758-T01-IM7')

partial_loss = full_prostate$SAMPLE_ID[which(full_prostate$classification == 'partial_loss')]
partial_loss = copynumber[which(copynumber$id %in% partial_loss & copynumber$chrom == 24 & copynumber$tcn.em == 0), ]
partial_loss[,c('id', "start", 'end')]


pLoss = c('P-0001698-T01-IM3',
'P-0016450-T01-IM6',
'P-0017964-T01-IM6',
'P-0018217-T01-IM6',
'P-0018798-T01-IM6',
'P-0021792-T01-IM6',
'P-0028784-T01-IM6',
'P-0029339-T01-IM6',
'P-0033767-T01-IM6',
'P-0053666-T01-IM6',
'P-0058347-T01-IM6',
'P-0063708-T01-IM7',
'P-0069465-T01-IM7',
'P-0074727-T01-IM7',
'P-0077871-T01-IM7')

pLoss_density = partial_loss[which(partial_loss$id %in% pLoss), ]
pLoss_density = pLoss_density[-1,]
pLoss_vector = data.frame()
for(i in 1:nrow(pLoss_density)){
  print(i)
  sample = pLoss_density$id[i]
  value = seq(from = pLoss_density$start[i], to = pLoss_density$end[i], by = 1000)
  out = data.frame(sample = sample,
                   vec = value)
  pLoss_vector = rbind(pLoss_vector, out)
}

dim(pLoss_vector)
head(pLoss_vector)
head(pLoss_density)
plot(density(pLoss_vector$vec))



IGV_pLoss = IGV[which(IGV$ID %in% c(pGL) & IGV$chrom == 24), ]
breakpoints = c()
for(i in unique(IGV_pLoss$ID)){
  n = length(IGV_pLoss$seg.mean[which(IGV_pLoss$ID == i)])
  m = n-1
  if(m == 0) next
  loc.end = IGV_pLoss$loc.end[which(IGV_pLoss$ID == i)][seq(1,m,1)]
  breakpoints = c(breakpoints, loc.end)
}

plot(density(breakpoints), 
     ylab = '',
     xlab = '',
     xaxt = 'n',
     yaxt = 'n',
     main = '', lwd = 2)
axis(side = 1, at = c(1e+07, 1.5e+07, 2e+07), labels = c('10Mb', '15Mb', '20Mb'))
mtext(text = 'Genomic Position', side = 1, line = 1.8)
mtext(text = 'Density', side = 2, line = 0.8)
mtext(text = 'Gain/Loss; Chromosome Y', side = 3, line = 0.8, adj = 0, cex = 1.5)
text(x = 1.98e+07, y = 0.2e-06, label = 'KDM5D')
box(lwd = 2)



##----------------+
## univariate model
##----------------+
LOY = c(full_prostate$SAMPLE_ID[which(full_prostate$classification == 'complete_loss')],
        pLoss, pGL)

full_prostate$Y_call = ifelse(full_prostate$SAMPLE_ID %in% LOY, 'loss', 'intact')
os_prostate = full_prostate[,c(colnames(full_prostate)[1:32], 'TP53_mut', 'TP53_Deletion', 'Y_call')]
osp = os_prostate
tmperg = read.csv('Data/05_Mutation/TMPRSS2_ETS.tsv', sep = '\t')
tmperg = as.character(tmperg$Sample.ID)
osp$SAMPLE_TYPE[osp$SAMPLE_TYPE == 'Local Recurrence'] = 'Primary'
osp$TMPERG = NA
for(i in 1:nrow(osp)){
  if(osp$SAMPLE_ID[i] %in% tmperg){
    osp$TMPERG[i] = 1
  } else {
    osp$TMPERG[i] = 0
  }
}
osp$TMPERG = as.integer(as.numeric(osp$TMPERG))
osp$TP53_Deletion = as.numeric(as.integer(osp$TP53_Deletion))

osp$TP53_Status = ifelse(osp$TP53_Deletion == 1, 1,
                               ifelse(osp$TP53_mut == 1, 1, 0))
osp$TP53_Deletion = NULL
osp$TP53_mut = NULL
osp$TP53_Status = as.integer(as.numeric(osp$TP53_Status))
osp$Y_call = as.factor(osp$Y_call)
osp$SAMPLE_TYPE = as.factor(osp$SAMPLE_TYPE)
osp$OS_Status_INT = ifelse(osp$OS_Status == 'Dead', 2L, 1L)



summary(coxph(formula = Surv(as.numeric(OS_months), as.numeric(factor(OS_Status_INT))) ~ Y_call+SAMPLE_TYPE+TP53_Status+Age_Sequencing+TMPERG, 
      data = osp))

summary(coxph(formula = Surv(as.numeric(OS_months), as.numeric(factor(OS_Status_INT))) ~ Y_call, 
              data = osp))

fit.prostate = survfit(Surv(as.numeric(OS_months), as.numeric(factor(OS_Status_INT))) ~ Y_call, data = osp)
prostate_KM = ggsurvplot(
  fit.prostate,
  data = osp,
  size = 1,
  palette = 
    c("#FF0000", "#104E8B"),
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  xlab = "Time (Months)",
  ylab = "Overall survival (%)",
  xlim = c(0,96),
  ylim = c(0, 1.01),
  break.time.by = 24,
  axes.offset = FALSE,
  font.legend =
    c(22),
  font.x =
    c(26, "bold"),
  font.y =
    c(26,"bold"),
  legend =
    c(.07, 0.4),
  legend.labs =
    c("wild type", "loss"),
  legend.title = "Chromosome Y",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  fontsize = 8
)

prostate_KM$plot = prostate_KM$plot +
  ggplot2::annotate("text",
                    x = 1, 
                    y = 0.15,
                    label = "HR = 1.48 (1.28, 1.70)", 
                    hjust=0, 
                    size=10)
prostate_KM$plot = prostate_KM$plot +
  ggplot2::annotate("text",
                    x= 1, 
                    y=0.08,
                    label = "P: 6.39 %*% 10^-8",
                    hjust=0,
                    size=10,
                    parse = TRUE)

ggsave_golden(filename = 'Figures_original/ProstateKM_pdf', plot = as_ggplot(ggplotGrob(prostate_KM)), width = 12)

ggsave(filename = 'Figures_original/ProstateKM_pdf', plot = prostate_KM, device = 'pdf')

##----------------+
## multivariate Cox
##----------------+
osp$Y_call = relevel(osp$Y_call, ref = "intact")
osp$TP53_Status = as.factor(osp$TP53_Status)
osp$TP53_Status = relevel(osp$TP53_Status, ref = '0')
osp$TMPERG = as.factor(osp$TMPERG)
osp$TMPERG = relevel(osp$TMPERG, ref = '0')
osp$SAMPLE_TYPE = as.factor(osp$SAMPLE_TYPE)
osp$SAMPLE_TYPE = relevel(osp$SAMPLE_TYPE, ref = 'Primary')
osp$MSI_TYPE[osp$MSI_TYPE == 'Do not report'] = NA
osp$MSI_TYPE[osp$MSI_TYPE == 'Indeterminate'] = NA
osp$MSI_TYPE = as.factor(osp$MSI_TYPE)
osp$MSI_TYPE = relevel(osp$MSI_TYPE, ref = 'Stable')


median(osp$OS_months, na.rm = T)
panels <- list(
  list(width = 0.03),
  list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
  list(width = 0.1, display = ~level),
  list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.55, item = "forest", hjust = 0.5, heading = "Hazard ratio", linetype = "dashed",
    line_x = 0
  ),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA),
  list(
    width = 0.05,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p"
  ),
  list(width = 0.03)
)


osp_forest = forest_model(coxph(Surv(as.numeric(OS_months), as.numeric(factor(OS_Status_INT))) ~ Y_call+SAMPLE_TYPE+MSI_TYPE+Age_Sequencing+TP53_Status+TMPERG, data = osp), panels)

osp_forest = osp_forest +
  theme(
    axis.text.x = element_text(size = 18,
                               face = "bold",
                               color = "black"),
    aspect.ratio = 0.5
  )

ggsave_golden(filename = 'Figures_original/CoxPH_prostate.pdf', plot = osp_forest, width = 14.5)
summary(coxph(Surv(as.numeric(OS_months), as.numeric(factor(OS_Status_INT))) ~ Y_call+SAMPLE_TYPE+MSI_TYPE+Age_Sequencing+TP53_Status+TMPERG, data = osp))




