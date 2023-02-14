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
library(cowplot)

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



##------- Identify all PRAD samples with KDM5D loss
gain_loss = full_prostate$SAMPLE_ID[which(full_prostate$classification == 'gain_loss')]
gain_loss_prad = copynumber[which(copynumber$id %in% gain_loss & copynumber$chrom == 24 & copynumber$tcn.em == 0), ]
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

IGV_pLoss = IGV[which(IGV$ID %in% c(pGL, pLoss) & IGV$chrom == 24), ]

pvec = data.frame()
for(i in 1:nrow(IGV_pLoss)){
  category = ifelse(IGV_pLoss$seg.mean[i] >= 0.1, 'gain', 
                    ifelse(IGV_pLoss$seg.mean[i] <= -0.1, 'loss', 'wt'))
  breakpoint = (IGV_pLoss$loc.end[i] - IGV_pLoss$loc.start[i]) / 2
  loc.end = IGV_pLoss$loc.end[i]
  out = data.frame(alteration = category,
                   value = breakpoint,
                   loc.end = loc.end)
  pvec = rbind(pvec, out)
}

plot(density(pvec$loc.end[which(pvec$alteration == 'loss')]),
     ylab = '',
     xlab = '',
     xlim = c(5.0e06, 3.0e07),
     ylim = c(0, 2e-07),
     xaxt = 'n',
     yaxt = 'n',
     lwd = 2,
     main = '',
     col = 'darkblue')

lines(density(pvec$loc.end[which(pvec$alteration == 'gain')]),
     ylab = '',
     xlab = '',
     lwd = 2,
     col = 'darkred')

lines(density(pvec$loc.end[which(pvec$alteration == 'wt')]),
      ylab = '',
      xlab = '',
      lwd = 1,
      col = 'grey35')

legend('topleft',
       legend = c("Gain", 
                  "Loss",
                  'wild type'),
       col = c("darkred", "darkblue", 'grey35'), 
       lty = 1, lwd = 2.5, cex = 1,
       box.lty = 0)

axis(side = 1, at = c(5.0e06, 1e+07, 1.5e+07, 2e+07, 2.5e07, 3.0e+07), 
     labels = c('5Mb', '10Mb', '15Mb', '20Mb', '25Mb', '30Mb'))

mtext(text = 'Genomic Position', side = 1, line = 1.8)
mtext(text = 'Density', side = 2, line = 0.8)
mtext(text = 'Gain/Loss; Chromosome Y', side = 3, line = 0.8, adj = 0, cex = 1.5)
text(x = 2.16e+07, y = 0.12e-06, label = 'KDM5D', cex = 1.5, font = 3)
box(lwd = 2)

# breakpoints = c()
# for(i in unique(IGV_pLoss$ID)){
#   n = length(IGV_pLoss$seg.mean[which(IGV_pLoss$ID == i)])
#   m = n-1
#   if(m == 0) next
#   loc.end = IGV_pLoss$loc.end[which(IGV_pLoss$ID == i)][seq(1,m,1)]
#   breakpoints = c(breakpoints, loc.end)
# }
# 
# plot(density(breakpoints), 
#      ylab = '',
#      xlab = '',
#      xaxt = 'n',
#      yaxt = 'n',
#      main = '', lwd = 2)
# axis(side = 1, at = c(1e+07, 1.5e+07, 2e+07), labels = c('10Mb', '15Mb', '20Mb'))
# mtext(text = 'Genomic Position', side = 1, line = 1.8)
# mtext(text = 'Density', side = 2, line = 0.8)
# mtext(text = 'Gain/Loss; Chromosome Y', side = 3, line = 0.8, adj = 0, cex = 1.5)
# text(x = 1.98e+07, y = 0.2e-06, label = 'KDM5D')
# box(lwd = 2)




##------- Alteration frequency in wild type VS complete-loss cases
LOY = c(full_prostate$SAMPLE_ID[which(full_prostate$classification == 'complete_loss')],
        pLoss, pGL)
full_prostate$Y_call = ifelse(full_prostate$SAMPLE_ID %in% LOY, 'loss', 'intact')
full_prostate = full_prostate[,c(colnames(full_prostate)[1:41], 'Y_call', colnames(full_prostate)[42:1050])]

pro_out = data.frame()
for(i in unique(full_prostate$Y_call)){
  for(j in 43:ncol(full_prostate)){
    print(j)
    alteration = colnames(full_prostate)[j]
    freq = sum(full_prostate[which(full_prostate$Y_call == i), j], na.rm = T)
    rel_freq = freq / length(full_prostate$SAMPLE_ID[which(full_prostate$Y_call == i)])
    category = i
    out = data.frame(category = category,
                     alteration = alteration,
                     frequency = rel_freq)
    pro_out = rbind(pro_out, out)
  }
}

pro_summary = data.frame()
for(i in unique(pro_out$alteration)){
  complete_loss = pro_out$frequency[which(pro_out$alteration == i & pro_out$category == 'loss')]
  wt = pro_out$frequency[which(pro_out$alteration == i & pro_out$category == 'intact')]
  alteration = i
  out = data.frame(alteration = alteration,
                   wt = wt,
                   complete_loss = complete_loss)
  pro_summary = rbind(pro_summary, out)
}
pro_summary$display = ifelse(abs(pro_summary$wt - pro_summary$complete_loss) >= 0.05, 'plot', 'notplot')
pro_summary$display[which(pro_summary$alteration == 'MYC_Amplification')] = 'plot'

alteration_frequency = ggplot(pro_summary, aes(x = wt, y = complete_loss, color = display, label = alteration)) +
  geom_text_repel(aes(label = ifelse(pro_summary$display == 'plot', pro_summary$alteration, ''))) +
  geom_jitter(size = 3) +
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

alteration_frequency
ggsave_golden(filename = 'Figures_original/Prostate_AlterationFreq.pdf', plot = alteration_frequency, width = 6)



##------- Gene-level: GLM
full_prostate$TP53_Deletion = as.numeric(as.integer(full_prostate$TP53_Deletion))
full_prostate$TP53_Status = ifelse(full_prostate$TP53_Deletion == 1, 1,
                                   ifelse(full_prostate$TP53_mut == 1, 1, 0))
full_prostate$TP53_Deletion = NULL
full_prostate$TP53_mut = NULL
full_prostate$TP53_Status = as.integer(as.numeric(full_prostate$TP53_Status))

glm_prostate = data.frame()
for(i in 43:1049){
  try({
    print(i)
    loss = 240
    intact = 1557
    n_altered_all = sum(full_prostate[,colnames(full_prostate)[i]], na.rm = T)
    rel_all = (n_altered_all / 1797) * 100
    n_altered_LOY = sum(full_prostate[which(full_prostate$Y_call == 'loss'), colnames(full_prostate)[i]], na.rm = T)
    rel_loy = (n_altered_LOY / loss)*100
    n_altered_wt = sum(full_prostate[which(full_prostate$Y_call == 'intact'), colnames(full_prostate)[i]], na.rm = T)
    rel_wt = (n_altered_wt / intact) * 100
    
    formula = as.formula(paste0(colnames(full_prostate)[i], "~ Y_call + ploidy + fraction_cna + purity + TP53_Status + MSI_TYPE"))
    log_results = glm(formula = formula, data = full_prostate, family = binomial)
    log_results_df = as.data.frame(summary(log_results)$coefficients)
    log_results_df$variable = row.names(log_results_df)
    log_results_df$gene = colnames(full_prostate)[i]
    colnames(log_results_df)[1:4] = c("estimate", "std_err", "z_value", "p_value")
    conf_df = as.data.frame(confint.default(log_results))
    conf_df$variable = row.names(conf_df)
    log_results_df = left_join(log_results_df, conf_df, by = "variable")
    log_results_df = log_results_df[c(5,6,1:4,7,8)]
    log_results_df$rel_all = rel_all
    log_results_df$rel_loy = rel_loy
    log_results_df$rel_wt = rel_wt
    glm_prostate = rbind(glm_prostate, log_results_df)
  })
}
glm_prostate$p_adj = p.adjust(p = glm_prostate$p_value, method = 'fdr')
pro = glm_prostate[which(glm_prostate$variable == 'Y_callloss' & glm_prostate$rel_all > 3), ]


gene_glm_prostate = glm_prostate[which(glm_prostate$p_value <= 0.05 & glm_prostate$variable == 'Y_callloss' & glm_prostate$rel_all >= 3), ]
length_genes = length(unique(gene_glm_prostate$gene))
pos = position_jitter(width = 0.2, seed = 2)

gene_glm_plot = ggplot(gene_glm_prostate, 
                       aes(x = reorder(gene, estimate),
                           y = estimate,
                           label = p_value)) +
  geom_pointrange(aes(ymin = estimate - std_err,
                      ymax = estimate + std_err),
                  size = 0.75,
                  position = pos,
                  fatten = 4) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey35', size = 0.25) +
  geom_vline(xintercept = seq(1.5, length_genes, 1), linetype = 'dashed', color = 'grey35', size = 0.4) +
  scale_y_continuous(expand = c(0.1, 0),
                     limits = c(-2, 2),
                     sec.axis = dup_axis()) +
  geom_text_repel(aes(label = paste0('P=', round(p_value, 3))), position = pos) +
  coord_flip() +
  theme_std(base_size = 14) +
  theme(legend.position = 'none',
        axis.line = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = '', y = 'log ODDS') +
  panel_border(size = 2, color = 'black')

ggsave_golden(filename = 'Figures_original/GeneLevel_Prostate.pdf', plot = gene_glm_plot, width = 7)


##------- Univariate CoxPH MODEL
osp = full_prostate[,c(colnames(full_prostate)[1:32], 'TP53_Status', 'Y_call')]
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
osp$Y_call = as.factor(osp$Y_call)
osp$SAMPLE_TYPE = as.factor(osp$SAMPLE_TYPE)
osp$OS_Status_INT = ifelse(osp$OS_Status == 'Dead', 2L, 1L)


##------- Model
summary(coxph(formula = Surv(as.numeric(OS_months), as.numeric(factor(OS_Status_INT))) ~ Y_call, 
              data = osp))

fit.prostate = survfit(Surv(as.numeric(OS_months), as.numeric(factor(OS_Status_INT))) ~ Y_call, data = osp)
prostate_KM = ggsurvplot(
  fit.prostate,
  data = osp, size = 1, palette = c("#FF0000", "#104E8B"), conf.int = FALSE, pval = FALSE,
  risk.table = TRUE,
  xlab = "Time (Months)",
  ylab = "Overall survival (%)",
  xlim = c(0,96),
  ylim = c(0, 1.01),
  break.time.by = 24,
  axes.offset = FALSE, 
  font.legend = c(22),
  font.x = c(26, "bold"),
  font.y = c(26,"bold"),
  legend = c(.07, 0.4),
  legend.labs = c("wild type", "loss"),
  legend.title = "Chromosome Y",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  fontsize = 8)

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

ggsave(filename = 'Figures_original/ProstateKM_pdf', plot = prostate_KM, device = 'pdf')




##----------------+
## multivariate Cox
##----------------+
osp = merge(osp, full_prostate[,c('SAMPLE_ID', 'RB1_mut', 'RB1_Deletion', 'ATM_mut', 'ATM_Deletion', 'SPOP_mut', 'SPOP_Amplification')],
            by = 'SAMPLE_ID', all.x = T)
osp$RB1_Status = ifelse(osp$RB1_mut == 1, '1',
                        ifelse(osp$RB1_Deletion == 1, '1', 0))
osp$RB1_Status = as.factor(osp$RB1_Status)
osp$ATM_Status = ifelse(osp$ATM_Deletion == 1, '1', 
                        ifelse(osp$ATM_mut == 1, '1', 0))
osp$ATM_Status = as.factor(osp$ATM_Status)
osp$SPOP_Status = ifelse(osp$SPOP_Amplification == 1, '1',
                         ifelse(osp$SPOP_mut == 1, '1', 0))
osp$SPOP_Status = as.factor(osp$SPOP_Status)
osp$CANCER_TYPE_DETAILED_ritika = as.factor(osp$CANCER_TYPE_DETAILED_ritika)

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
osp$RB1_Status = relevel(osp$RB1_Status, ref = '0')
osp$ATM_Status = relevel(osp$RB1_Status, ref = '0')
osp$SPOP_Status = relevel(osp$RB1_Status, ref = '0')
osp$CANCER_TYPE_DETAILED_ritika = relevel(osp$CANCER_TYPE_DETAILED_ritika, ref = 'Prostate Adenocarcinoma (PRAD)')


##-- Model
summary(coxph(formula = Surv(as.numeric(OS_months), as.numeric(factor(OS_Status_INT))) ~ Y_call+SAMPLE_TYPE+MSI_TYPE+CANCER_TYPE_DETAILED_ritika+Age_Sequencing+TP53_Status+TMPERG+RB1_Status+ATM_Status+SPOP_Status, 
              data = osp))

##-- Visualization
median(osp$OS_months, na.rm = T)
panels <- list(
  list(width = 0.01),
  list(width = 0.5, display = ~variable, fontface = "bold", heading = "Variable"),
  list(width = 0.75, display = ~level, heading = 'Level', text_size = 2),
  list(width = 0.15, item = 'vline', hjust = 0.5),
  list(width = 0.005, display = ~n, hjust = 1, heading = "N"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.8, item = "forest", hjust = 0.5, heading = "Hazard ratio", linetype = "dashed", line_x = 0),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.8, display = ~ ifelse(reference, "Reference", sprintf("%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA, heading = 'HR (95% CI)'),
  list(width = 0.05,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
    display_na = T, hjust = 1, heading = "P-value"),
  list(width = 0.03, text_size = 1)
)

osp_forest = forest_model(coxph(formula = Surv(as.numeric(OS_months), as.numeric(factor(OS_Status_INT))) ~ 
                                  Y_call+SAMPLE_TYPE+MSI_TYPE+CANCER_TYPE_DETAILED_ritika+Age_Sequencing+TP53_Status+TMPERG+RB1_Status, 
                                data = osp), panels, recalculate_width = F)

osp_forest = osp_forest + theme(
    axis.text.x = element_text(size = 14,
                               face = "bold",
                               color = "black"),
    aspect.ratio = 0.5
  )

osp_forest

ggsave_golden(filename = 'Figures_original/CoxPH_prostate.pdf', plot = osp_forest, width = 16)



