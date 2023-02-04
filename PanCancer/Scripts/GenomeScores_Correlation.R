##----------------+
## Genome-scores; and
## correlations with 
## LOY
## - Purity (bin;cases)
## - Purity: sample coverage
## - LOY correlation with FGA/#arm losses/Aneuploidy score
## - LOY associated with WGD
## - MSI and LOY investigation
##----------------+
##
## start: 11/01/2022
## revision: 01/10/2023
## revision: 01/12/2023
## revision: 01/13/2023
## revision: 01/20/2023
## revision: 02/04/2023
## 
## chris-kreitzer


clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')

library(ggrepel)
library(cowplot)
library(patchwork)
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')


##----------------+
## both QC TRUE and FALSE
## included;
## 
## TODO
## - define cancer types to include
## - work on QC TRUE samples
##----------------+
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]

##----------------+
## purity distribution
## - LOY enriched in low tumor purity
## - LOY enriched in high tumor purity?
##----------------+
Y_calls = cohort
Y_calls$purity[which(Y_calls$purity == 0)] = 0.0001
Y_calls = Y_calls[!is.na(Y_calls$purity), ]

purity_bin = data.frame()
for(i in seq(0.1, 0.9, 0.1)){
  #print(c(i, i+0.2))
  data_sub = Y_calls[which(Y_calls$purity > i & Y_calls$purity <= i + 0.1), ]
  n_loss = length(data_sub$SAMPLE_ID[which(data_sub$classification %in% c('complete_loss'))])
  n_loss_rel = n_loss / nrow(data_sub)
  out = data.frame(purity_bin = i,
                   total = nrow(data_sub),
                   n_loss_rel = n_loss_rel)
  purity_bin = rbind(purity_bin, out)
}

purity_bin$purity_bin = paste0('[', purity_bin$purity_bin, '-', purity_bin$purity_bin+0.1, ']')

cases_purity = ggplot(purity_bin, aes(x = purity_bin, y = total)) +
  geom_bar(stat = 'identity', color = 'grey25', fill = 'grey25') +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_blank()) +
  labs(x = '', y = 'cases (n)')

procent_purity = ggplot(purity_bin, aes(x = purity_bin, y = n_loss_rel)) +
  geom_bar(stat = 'identity', color = 'grey25', fill = 'grey25') +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'purity-group', y = 'Y-chromosome loss [%]')
procent_purity

purity_compact = cases_purity / procent_purity
ggsave_golden(filename = 'Figures_original/LOY_by_Purity.pdf', plot = purity_compact, width = 10)



##----------------+
## Average purity and
## LOY based on cancer types
##----------------+
purity_out = data.frame()
for(i in unique(cohort$CANCER_TYPE)){
  cancer = i
  #cancer = gsub("\\s*\\([^\\)]+\\)", "", i)
  data.sub = cohort[which(cohort$CANCER_TYPE== i), ]
  data.sub = data.sub[!is.na(data.sub$purity), ]
  n = nrow(data.sub)
  out = data.frame(cancer = cancer,
                   n = n,
                   value = data.sub$purity,
                   median_value = median(data.sub$purity, na.rm = T))
  purity_out = rbind(purity_out, out)
  rm(data.sub)
}
purity_out = purity_out[which(purity_out$cancer %in% ctypes_keep), ]
purity_summary = purity_out[,c('cancer', 'n', 'median_value')]
purity_summary = unique(purity_summary)
purity_summary = purity_summary[order(purity_summary$median_value, decreasing = F), ]
fn = factor(purity_summary$cancer, levels = purity_summary$cancer)

##-------
## Visualization
##-------
Purity.plot = ggplot(purity_out, 
                     aes(x = reorder(cancer, median_value), y = value)) +
  geom_quasirandom(dodge.width = 0.6, cex = 1, nbins = 50, alpha = 0.5) +
  stat_summary(fun = median, 
               geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..), 
               position = position_dodge(width = 0.5), 
               linewidth = 0.8, 
               col = 'red',
               size = 1) +
  scale_x_discrete(breaks = purity_summary$cancer,
                   labels = paste0(purity_summary$cancer, '\n(n=', purity_summary$n, ')')) +
  theme_std(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.25),
                     labels = seq(0, 1, 0.25)) +
  labs(x = '', y = 'Purity')

Purity.plot


##-------
## Add info for LOY
##-------
fLOY = data.frame()
for(i in unique(cohort$CANCER_TYPE)){
  print(i)
  cancer = i
  #cancer = gsub("\\s*\\([^\\)]+\\)", "", i)
  data.sub = cohort[which(cohort$CANCER_TYPE == i), ]
  n = nrow(data.sub)
  n_loss = length(data.sub$SAMPLE_ID[which(data.sub$classification %in% c('complete_loss', 'partial_loss', 'relative_loss'))])
  out = data.frame(cancer = cancer,
                   n = n,
                   value = (n_loss / n) * 100)
  fLOY = rbind(fLOY, out)
  rm(data.sub)
}

fLOY = fLOY[which(fLOY$cancer %in% ctypes_keep), ]
fLOY$cancer = factor(fLOY$cancer, levels = levels(fn))

fLOY_plot = ggplot(fLOY, 
                   aes(x = cancer, y = value)) +
  geom_bar(stat = 'identity', fill = 'black', color = 'black') +
  geom_hline(yintercept = seq(0, 100, 25), linewidth = 0.25, color = 'grey55', linetype = 'dashed') +
  scale_y_continuous(expand = c(0, 0)) +
  theme_std(base_size = 14) +
  theme(axis.text.x = element_blank()) +
  labs(x = '', y = 'LOY [%]')


purity_fLOY = fLOY_plot / Purity.plot
ggsave_golden(filename = 'Figures_original/Purity_fractionLOY_cancer.pdf', plot = purity_fLOY, width = 12)


##----------------+
## correlation LOY and purity
##----------------+
purity_loy = merge(purity_summary, fLOY, by = 'cancer', all.x = T)
colnames(purity_loy) = c('cancer', 'n_purity', 'median_purity', 'n_loy', 'median_loy')

mpurity_mloy = ggplot(purity_loy, aes(x = median_purity,
                       y = median_loy,
                       label = cancer)) +
  geom_jitter() +
  geom_text_repel(size = 1.5) +
  scale_x_continuous(limits = c(0, 1), 
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100),
                     expand = c(0,0)) +
  geom_abline(intercept = 100, slope = -100, linetype = 'dashed', linewidth = 0.35, color = 'darkgrey') +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1) +
  panel_border(color = 'black', size = 1.5) +
  labs(x = 'median purity', y = 'median LOY') +
  annotate(geom = 'text', x = 0.85, y = 85, label = expression(paste(rho, ' = -0.49')), size = 6)

ggsave_golden(filename = 'Figures_original/Purity_LOY_correlation.pdf', plot = mpurity_mloy, width = 8)




##----------------+
## LOY in lower pure samples
## associated with sample coverage?
##----------------+
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
low_purity_samples = cohort$counts_file[which(cohort$purity > 0 & cohort$purity < 0.3)]
write.table(low_purity_samples, file = 'Data/01_Coverage_Depth/Low_purity[0:0.29]_samples.txt', sep = '\t', row.names = F, quote = F, col.names = 'sample')
high_purity_samples = cohort$counts_file[which(cohort$purity > 0.8 & cohort$purity < 1)]
write.table(high_purity_samples, file = 'Data/01_Coverage_Depth/High_purity[0.8:1]_samples.txt', sep = '\t', row.names = F, quote = F, col.names = 'sample')


##----------------+
## Investigate the results
##----------------+
low_purity = read.csv('Data/01_Coverage_Depth/Average_Depth_lowPurity.txt', sep = '\t')
low_purity$group = 'low_purity'
high_purity = read.csv('Data/01_Coverage_Depth/Average_Depth_highPurity.txt', sep = '\t')
high_purity$group = 'high_purity'

purity_coverage = rbind(low_purity, high_purity)

Purity_Coverage_plot = ggplot(purity_coverage, aes(x = group, y = average_depth_TUM)) +
  geom_violin(width = 0.45) +
  geom_jitter(position = position_dodge2(width = 0.15), size = 0.25) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_std(base_size = 14, base_line_size = 0.5) +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, linewidth = 2)) +
  labs(x = '', y = 'Average sequencing depth')
  
ggsave_golden(filename = 'Figures_original/Purity_Coverage.pdf', plot = Purity_Coverage_plot, width = 6)  




##----------------+
## Correlation analysis;
## FGA, arm-losses VS 
## Fraction LOY; 
## purity all
##----------------+
clean()
gc()
.rs.restartR()
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]

genome_scores = data.frame()
for(i in unique(cohort$CANCER_TYPE)){
  type = i
  Aneuploidy_score = mean(x = cohort$AS_score[which(cohort$CANCER_TYPE == i)], na.rm = T)
  Loss_score = mean(x = cohort$losses_n[which(cohort$CANCER_TYPE == i)], na.rm = T)
  FGA_score = mean(x = cohort$fraction_cna[which(cohort$CANCER_TYPE == i)], na.rm = T)
  n_loy = length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE == i & cohort$classification %in% c('complete_loss', 'partial_loss', 'relative_loss'))])
  fraction_LOY = n_loy / length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE == i)])
  out = data.frame(CancerType = type,
                   Aneuploidy_score = Aneuploidy_score,
                   Loss_score = Loss_score,
                   FGA_score = FGA_score,
                   n = length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE == i)]),
                   n_loy = n_loy,
                   fraction_LOY = fraction_LOY*100)
  genome_scores = rbind(genome_scores, out)
}

genome_scores = genome_scores[genome_scores$CancerType %in% ctypes_keep, ]
write.table(genome_scores, file = 'Data/05_Mutation/011823/GenomeScores_correlation.txt', sep = '\t', row.names = F)


##-------
## Visualization;
## - Aneuploidy score
##-------
Aneuploidy_correlation = ggplot(genome_scores, 
                                aes(x = fraction_LOY, 
                                    y = Aneuploidy_score)) +
  geom_point(data = genome_scores, aes(size = n), shape = 20) +
  scale_size(breaks = c(100, 500, 1000, 1500), range = c(0, 8),
             labels = c(100, 500, 1000, '>1500')) +
  scale_x_continuous(expand = c(0.01, 0), 
                     limits = c(0, 80)) +
  geom_text_repel(aes(label = CancerType), 
                  size = as_points(7), 
                  point.padding = as_points(1)) +
  stat_smooth(method = lm)  +
  theme_std(base_size = 14, base_line_size = .5) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank()) +
  labs(y = 'Average # of chromosome arms gained & lost', x = 'LOY [%]', size = 'n patients')

Aneuploidy_correlation = Aneuploidy_correlation + 
  stat_cor(method = "spearman", label.x = 4, label.y = 30, size = 6) +
  panel_border(size = 2, color = 'black')


ggsave_golden(filename = 'Figures_original/ArmsGain_Loss_LOY.pdf', plot = Aneuploidy_correlation, width = 12)



##-------
## Visualization;
## - Loss score
##-------
Loss_correlation = ggplot(genome_scores, 
                          aes(x = fraction_LOY, 
                              y = Loss_score)) +
  geom_point(data = genome_scores, aes(size = n), shape = 20) +
  scale_size(breaks = c(100, 500, 1000, 1500), range = c(0, 8),
             labels = c(100, 500, 1000, '>1500')) +
  scale_x_continuous(expand = c(0.01, 0), 
                     limits = c(0, 80)) +
  geom_text_repel(aes(label = CancerType), 
                  size = as_points(7), 
                  point.padding = as_points(1)) +
  stat_smooth(method = lm)  +
  theme_std(base_size = 14, base_line_size = .5) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank()) +
  labs(y = 'Average # of chromosome arms lost', x = 'LOY [%]', size = 'n patients')

Loss_correlation = Loss_correlation + 
  stat_cor(method = "spearman", label.x = 4, label.y = 30, size = 6) +
  panel_border(size = 2, color = 'black')

ggsave_golden(filename = 'Figures_original/ArmsLost_LOY_Correlation.pdf', plot = Loss_correlation, width = 12)



##-------
## Visualization;
## - FGA score
##-------
FGA_correlation = ggplot(genome_scores, 
                         aes(x = fraction_LOY, 
                             y = FGA_score)) +
  geom_point(data = genome_scores, aes(size = n), shape = 20) +
  scale_size(breaks = c(100, 500, 1000, 1500), range = c(0, 8),
             labels = c(100, 500, 1000, '>1500')) +
  scale_x_continuous(expand = c(0.01, 0), 
                     limits = c(0, 80)) +
  geom_text_repel(aes(label = CancerType), 
                  size = as_points(7), 
                  point.padding = as_points(1)) +
  stat_smooth(method = lm)  +
  theme_std(base_size = 14, base_line_size = .5) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank()) +
  labs(y = 'Fraction Genome Altered', x = 'LOY [%]', size = 'n patients')
  
FGA_correlation = FGA_correlation + 
  stat_cor(method = "spearman", label.x = 4, label.y = 0.7, size = 6) +
  panel_border(size = 2, color = 'black')

GenomeScore_Correlation = Aneuploidy_correlation + theme(legend.position = 'none') + FGA_correlation
ggsave_golden(filename = 'Figures_original/GenomeScore_Correlation.pdf', plot = GenomeScore_Correlation, width = 16)


##----------------+
## Exclude low purity 
## samples; is purity
## driving the curves?
##----------------+


##----------------+
## WGD and LOY association;
##----------------+
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]

WGD_LOY = data.frame()
for(i in unique(cohort$genome_doubled)){
  if(i %in% c('TRUE', 'FALSE')){
    n = length(unique(cohort$SAMPLE_ID[which(cohort$genome_doubled == i)]))
    cc = cohort[which(cohort$genome_doubled == i), ]
    cc_df = as.data.frame(table(cc$classification))
    cc_df$group = i
    cc_df$rel_freq = cc_df$Freq / n
    WGD_LOY = rbind(WGD_LOY, cc_df)
  } else next
}
WGD_LOY$group[which(WGD_LOY$group == 'TRUE')] = '1x'
WGD_LOY$group[which(WGD_LOY$group == 'FALSE')] = 'none'
WGD_LOY$group = factor(WGD_LOY$group, levels = c('none', '1x'))
WGD_LOY$Var1 = factor(WGD_LOY$Var1, levels = rev(c('complete_loss',
                                                           'relative_loss',
                                                           'partial_loss',
                                                           'gain',
                                                           'partial_gain',
                                                           'gain_loss',
                                                           'wt')))

##-------
## Visualization
##-------
WGD_LOY_plot = ggplot(WGD_LOY, 
                      aes(x = group, 
                          y = rel_freq, 
                          fill = Var1)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'complete_loss' = '#0E3F7C',
                               'partial_loss' = '#00AEC8',
                               'relative_loss' = '#474089',
                               'gain' = '#D53833',
                               'partial_gain' = '#E3CC98',
                               'gain_loss' = '#f9f8d5'),
                    name = '') +
  scale_y_continuous(expand = c(0.01,0)) +
  #scale_x_discrete(labels = c('none', '1x')) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(aspect.ratio = 2) +
  labs(x = '# genome doublings', y = 'Fraction of tumors')

ggsave_golden(filename = 'Figures_original/Fraction_Y_loss_WGD.pdf', plot = WGD_LOY_plot, width = 4)



##----------------+
## MSI_Type; 
## it seems like that MSI
## tumors show a tendency for 
## retaining the Y chromosome
## - work with QC TRUE samples
## - work with Indeterminate, instable, and stable types
##----------------+ 
clean()
gc()
.rs.restartR()
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
clinical_glm = cohort[,c('SAMPLE_ID', 'CANCER_TYPE', 'classification', 'MSI_SCORE', 'MSI_TYPE')]
clinical_glm$classification[which(clinical_glm$classification %in% c('complete_loss', 'partial_loss', 'relative_loss'))] = 'LOY'
clinical_glm$classification[which(clinical_glm$classification %in% c('wt', 'partial_gain', 'gain_loss', 'gain'))] = 'wt'
clinical_glm$MSI_TYPE = as.character(as.factor(clinical_glm$MSI_TYPE))
clinical_glm$CANCER_TYPE = as.character(as.factor(clinical_glm$CANCER_TYPE))
colnames(clinical_glm)[3] = 'Y_call'
clinical_glm = clinical_glm[which(clinical_glm$MSI_TYPE %in% c('Instable', 'Stable')), ]
clinical_glm = clinical_glm[!is.na(clinical_glm$MSI_SCORE), ]
clinical_glm$MSI_TYPE = factor(clinical_glm$MSI_TYPE, levels = c('Instable', 'Stable'))

##-------
## MSI-MSS SCNA comparison
##-------
IMPACT_segs = read.csv('Data/04_Loss/010523/IGV_out.seg', sep = '\t')
IMPACT_instable = IMPACT_segs[which(IMPACT_segs$ID %in% clinical_glm$SAMPLE_ID[which(clinical_glm$MSI_TYPE == 'Instable')]), ]
IMPACT_stable = IMPACT_segs[which(IMPACT_segs$ID %in% clinical_glm$SAMPLE_ID[which(clinical_glm$MSI_TYPE == 'Stable')]), ]


##-------
## Can we use MSISensor
## for categorization
##-------
MSIsensor_plot = ggplot(clinical_glm, aes(x = MSI_TYPE, y = MSI_SCORE, color = MSI_TYPE)) +
  geom_quasirandom(shape = 16, size = 0.75, alpha = 0.75) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) quantile(x, .25),
               fun.ymax = function(x) quantile(x, .75),
               geom = 'errorbar', width = 0.1) +
  stat_summary_bin(geom = 'point', fun.y = 'mean', shape = 21, color = 'black', fill = 'white') + 
  scale_color_manual(values = c('Instable' = 'springgreen4',
                                'Stable' = 'grey75')) +
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey35', linewidth = 0.35) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_discrete(labels = c(paste0('MSI\n(n=', length(unique(clinical_glm$SAMPLE_ID[which(clinical_glm$MSI_TYPE == 'Instable')])),')'),
                              paste0('MSS\n(n=', length(unique(clinical_glm$SAMPLE_ID[which(clinical_glm$MSI_TYPE == 'Stable')])),')'))) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 2,
        legend.position = 'none') +
  labs(x = '', y = 'MSIsensor score', color = '')

ggsave_golden(filename = 'Figures_original/MSISensor.pdf', plot = MSIsensor_plot, width = 7)



##----------------+
## Tumor Types with both
## MSI and MSS;
## general overview
##----------------+
MSI_all = data.frame()
for(i in unique(clinical_glm$CANCER_TYPE)){
  print(i)
  try({
    cancer = i
    n = length(clinical_glm$SAMPLE_ID[which(clinical_glm$CANCER_TYPE == i)])
    n_instable = length(clinical_glm$SAMPLE_ID[which(clinical_glm$CANCER_TYPE == i & clinical_glm$MSI_TYPE == 'Instable')])
    stable_loy = table(clinical_glm$Y_call[which(clinical_glm$CANCER_TYPE == i & clinical_glm$MSI_TYPE == 'Stable')])[['LOY']]
    s_loy_fraction = stable_loy / sum(table(clinical_glm$Y_call[which(clinical_glm$CANCER_TYPE == i & clinical_glm$MSI_TYPE == 'Stable')]))
    instable_loy = table(clinical_glm$Y_call[which(clinical_glm$CANCER_TYPE == i & clinical_glm$MSI_TYPE == 'Instable')])[['LOY']]
    i_loy_fraction = instable_loy / sum(table(clinical_glm$Y_call[which(clinical_glm$CANCER_TYPE == i & clinical_glm$MSI_TYPE == 'Instable')]))
    out = data.frame(cancer = cancer,
                     category = c('MSS', 'MSI'),
                     value = c(s_loy_fraction, i_loy_fraction),
                     n = n,
                     n_MSI = n_instable)
    MSI_all = rbind(MSI_all, out)
  })
}

MSI_tumors = ggplot(MSI_all[which(MSI_all$n_MSI > 10), ], aes(x = category, y = value)) +
  geom_jitter(width = 0.2, size = 2) +
  stat_compare_means(label.y = 0.75) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) quantile(x, .25),
               fun.ymax = function(x) quantile(x, .75),
               geom = 'errorbar', width = 0.01, color = 'springgreen4') +
  stat_summary_bin(geom = 'point', fun.y = 'mean', shape = 21, color = 'springgreen4', fill = 'white', size = 2) +
  scale_y_continuous(limits = c(0, 0.8),
                     expand = c(0.01, 0),
                     labels = c(0, 20, 40, 60, 80)) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1.8,
        axis.line = element_blank()) +
  labs(x = '', y = 'LOY [%]') +
  panel_border(size = 2, color = 'black')


ggsave_golden(filename = 'Figures_original/MSI_LOY.pdf', plot = MSI_tumors, width = 6)


clinical_glm_msi = clinical_glm[which(clinical_glm$CANCER_TYPE %in% ctypes_keep), ] 
clinical_glm_msi = clinical_glm_msi[!is.na(clinical_glm_msi$Y_call), ]
clinical_glm_msi = clinical_glm_msi[which(clinical_glm_msi$CANCER_TYPE %in% c('Bladder Cancer', 'Colorectal Cancer', 'Esophagogastric Cancer',
                                                                              'Glioma', 'Non-Small Cell Lung Cancer', 'Prostate Cancer')), ]

MSI_LOY_plot = ggplot(clinical_glm_msi, aes(x = Y_call, y = MSI_SCORE, color = MSI_TYPE)) +
  geom_quasirandom(width = 0.25) +
  facet_wrap(~CANCER_TYPE, nrow = 1, scales = 'fixed') +
  stat_compare_means(comparisons = list(c('LOY', 'nonLOY')),
                     method = 'wilcox.test', label = 'p.signif', hide.ns = T, label.y = 47) +
  scale_color_manual(values = c('Instable' = 'springgreen4',
                                'Stable' = 'grey75')) +
  geom_vline(xintercept = 2.5, linetype = 'dashed', linewidth = 0.25) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_std(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(0.8), margin = margin(), angle = 0, hjust = 0.5),
        panel.spacing = unit(0, "pt"),
        legend.position = 'top',
        aspect.ratio = 1) +
  labs(x = '', y = 'MSIsensor score', color = '')
  
ggsave_golden(filename = 'Figures_original/MSI_LOY_association.pdf', plot = MSI_LOY_plot, width = 16)  






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
  try({
    da = clinical_glm[dplyr::between(x = clinical_glm$purity, left = i, right = i+0.1), ]
    if(all(da$classification %in% c('0', '1'))){
      frac = (table(da$classification)[[2]] / sum(table(da$classification))) *100
      print(frac)
      out = data.frame(group = i,
                       n = nrow(da),
                       frac = frac,
                       average_age = mean(da$Age_Sequencing, na.rm = T))
    } else {
      frac = NA
      out = data.frame(group = i,
                       n = nrow(da),
                       frac = frac,
                       average_age = mean(da$Age_Sequencing, na.rm = T))
    }
    purity_out = rbind(purity_out, out)
  })
}

purity_out$group = as.factor(as.numeric(purity_out$group))

#' Visualization
n.purity = ggplot(purity_out, aes(x = group, y = n)) +
  geom_bar(stat = 'identity', width = 0.8) +
  labs(x = '', y = 'cases (n)') + 
  theme_std(base_size = 14) +
  theme(axis.text.x = element_blank())
frac.purity = ggplot(purity_out, aes(x = group, y = frac)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_x_discrete(labels = c('[0-0.1]', '[0.1-0.2]', '[0.2-0.3]', '[0.3-0.4]', '[0.4-0.5]',
                              '[0.5-0.6]', '[0.6-0.7]', '[0.7-0.8]', '[0.8-0.9]', '[0.9-1]')) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(x = 'purity-group', y = 'Y-chromosome loss [%]') +
  theme_std(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

purity_clinical_glm = n.purity / frac.purity
ggsave_golden('Figures_original/Purity_cases_LOY.pdf', plot = purity_clinical_glm, width = 7)







#' out