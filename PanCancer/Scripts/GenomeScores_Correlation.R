##----------------+
## Genome-scores; and
## correlations with 
## LOY
## - Purity (bin;cases)
## - Purity: sample coverage
## - LOY correlation with FGA/#arm losses/Aneuploidy score
## - LOY associated with WGD
##----------------+
##
## start: 11/01/2022
## revision: 01/10/2023
## revision: 01/12/2023
## revision: 01/13/2023
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


##----------------+
## purity distribution
## - LOY enriched in low tumor purity
## - LOY enriched in high tumor purity?
##----------------+
Y_calls = cohort
Y_calls$purity[which(Y_calls$purity == 0)] = 0.0001
Y_calls = Y_calls[!is.na(Y_calls$purity), ]

purity_bin = data.frame()
for(i in seq(0, 0.9, 0.1)){
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
for(i in unique(cohort$CANCER_TYPE_ritika)){
  cancer = gsub("\\s*\\([^\\)]+\\)", "", i)
  data.sub = cohort[which(cohort$CANCER_TYPE_ritika == i), ]
  data.sub = data.sub[!is.na(data.sub$purity), ]
  n = nrow(data.sub)
  out = data.frame(cancer = cancer,
                   n = n,
                   value = data.sub$purity,
                   median_value = median(data.sub$purity, na.rm = T))
  purity_out = rbind(purity_out, out)
  rm(data.sub)
}

purity_summary = purity_out[,c('cancer', 'n', 'median_value')]
purity_summary = unique(purity_summary)
purity_summary$cancer = gsub("\\s*\\([^\\)]+\\)", "", purity_summary$cancer)
purity_summary = purity_summary[order(purity_summary$median_value, decreasing = F), ]
fn = factor(unique(x$cancer), levels = purity_summary$cancer)

##-------
## Visualization
##-------
Purity.plot = ggplot(purity_out, aes(x = reorder(cancer, median_value), y = value)) +
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


##-------
## Add info for LOY
##-------
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
fLOY = data.frame()
for(i in unique(cohort$CANCER_TYPE_ritika)){
  cancer = gsub("\\s*\\([^\\)]+\\)", "", i)
  data.sub = cohort[which(cohort$CANCER_TYPE_ritika == i), ]
  n = nrow(data.sub)
  n_loss = length(data.sub$SAMPLE_ID[which(data.sub$classification %in% c('complete_loss', 'partial_loss', 'relative_loss'))])
  out = data.frame(cancer = cancer,
                   n = n,
                   value = n_loss / n * 100)
  fLOY = rbind(fLOY, out)
  rm(data.sub)
}

fLOY$cancer = factor(fLOY$cancer, levels = levels(fn))

fLOY_plot = ggplot(fLOY, aes(x = cancer, y = value)) +
  geom_bar(stat = 'identity', fill = 'black', color = 'black') +
  geom_hline(yintercept = seq(0, 100, 25), linewidth = 0.25, color = 'grey55', linetype = 'dashed') +
  scale_y_continuous(expand = c(0, 0)) +
  theme_std(base_size = 14) +
  theme(axis.text.x = element_blank()) +
  labs(x = '', y = 'LOY [%]')


purity_fLOY = fLOY_plot / Purity.plot
ggsave_golden(filename = 'Figures_original/Purity_fractionLOY_cancer.pdf', plot = purity_fLOY, width = 12)


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

genome_scores = data.frame()
for(i in unique(cohort$CANCER_TYPE_ritika)){
  type = i
  Aneuploidy_score = mean(x = cohort$AS_score[which(cohort$CANCER_TYPE_ritika == i)], na.rm = T)
  Loss_score = mean(x = cohort$losses_n[which(cohort$CANCER_TYPE_ritika == i)], na.rm = T)
  FGA_score = mean(x = cohort$fraction_cna[which(cohort$CANCER_TYPE_ritika == i)], na.rm = T)
  n_loy = length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE_ritika == i & cohort$classification %in% c('complete_loss', 'partial_loss', 'relative_loss'))])
  fraction_LOY = n_loy / length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE_ritika == i)])
  out = data.frame(CancerType = type,
                   Aneuploidy_score = Aneuploidy_score,
                   Loss_score = Loss_score,
                   FGA_score = FGA_score,
                   n = length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE_ritika == i)]),
                   n_loy = n_loy,
                   fraction_LOY = fraction_LOY*100)
  genome_scores = rbind(genome_scores, out)
}
genome_scores = genome_scores[!genome_scores$CancerType %in% c('Ovary/Fallopian Tube (OVARY)', 'Vulva/Vagina (VULVA)'), ]
write.table(genome_scores, file = 'Data/05_Association/gene_level/GenomeScores.txt', sep = '\t', row.names = F)


##-------
## Visualization;
## - Aneuploidy score
##-------
Aneuploidy_correlation = ggplot(genome_scores, 
                                aes(x = fraction_LOY, 
                                    y = Aneuploidy_score)) +
  geom_point(data = genome_scores, aes(size = n), shape = 20) +
  scale_size_continuous(breaks = c(50, 100, 500, 1000, 1500, 3000),
                        labels = c(50, 100, 500, 1000, 1500, '>1500'),
                        name = 'n/CancerType') +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_text_repel(aes(label = CancerType), size = 3) +
  stat_smooth(method = lm) +
  #stat_smooth(method = loess, fullrange = FALSE, alpha = 0.1, span = 10) +
  theme_std(base_size = 14, base_line_size = 1) + 
  theme(aspect.ratio = 1) +
  labs(y = 'Average # of chr. arms gained & lost', x = 'Fraction LOY')

Aneuploidy_correlation


##-------
## Visualization;
## - Loss score
##-------
Loss_correlation = ggplot(genome_scores, 
                          aes(x = fraction_LOY, 
                              y = Loss_score)) +
  geom_point(data = genome_scores, aes(size = n), shape = 20) +
  scale_size_continuous(breaks = c(50, 100, 500, 1000, 1500, 3000),
                        labels = c(50, 100, 500, 1000, 1500, '>1500'),
                        name = 'n/CancerType') +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_text_repel(aes(label = CancerType), size = 3) +
  stat_smooth(method = lm) +
  #stat_smooth(method = loess, fullrange = FALSE, alpha = 0.1, span = 10) +
  theme_std(base_size = 14, base_line_size = 1) + 
  theme(aspect.ratio = 1) +
  labs(y = 'Average # of chr. arms lost', x = 'Fraction of male tumors with LOY')

Loss_correlation


##-------
## Visualization;
## - FGA score
##-------
FGA_correlation = ggplot(genome_scores, 
                         aes(x = fraction_LOY, 
                             y = FGA_score)) +
  geom_point(data = genome_scores, aes(size = n), shape = 20) +
  scale_size_continuous(breaks = c(50, 100, 500, 1000, 1500, 3000),
                        labels = c(50, 100, 500, 1000, 1500, '>1500'),
                        name = 'n/CancerType') +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_text_repel(aes(label = CancerType), size = 3) +
  stat_smooth(method = lm) +
  #stat_smooth(method = loess, fullrange = FALSE, alpha = 0.1, span = 10) +
  theme_std(base_size = 14, base_line_size = 1) + 
  theme(aspect.ratio = 1) +
  labs(y = 'Average Fraction Genome Altered', x = 'Fraction of male tumors with LOY')

FGA_correlation

Aneuploidy_correlation / Loss_correlation / FGA_correlation


##----------------+
## Exclude low purity 
## samples; is purity
## driving the curves?
##----------------+


##----------------+
## WGD and LOY association;
##----------------+

loy = data$IMPACT_Y_classification_final
clinics = data$IMPACT_clinicalAnnotation

eso = clinics[which(clinics$CANCER_TYPE == 'Esophagogastric Cancer'), 'SAMPLE_ID']
pan = clinics[which(clinics$CANCER_TYPE == 'Pancreatic Cancer'), 'SAMPLE_ID']
pro = clinics[which(clinics$CANCER_TYPE == 'Prostate Cancer'), 'SAMPLE_ID']
gli = clinics[which(clinics$CANCER_TYPE == 'Glioma'), 'SAMPLE_ID']
table(loy$classification[which(loy$sample %in% gli)])
wgd = read.csv('Data/04_Loss/QC_metrics.txt', sep = '\t')
wgd = wgd[which(wgd$ID %in% loy$sample), ]
head(wgd)


wgd_loy = merge(loy, wgd, by.x = 'sample', by.y = 'ID', all.x = T)


head(wgd_loy)
sum(table(wgd_loy$classification[which(wgd_loy$wgd == T)]))

wgd_plot = data.frame(category = rep(c('gain', 'gain_loss', 'loss', 'relative_loss', 'wt'), 2),
                      wgd = c(rep('NONE', 5), rep(1, 5)),
                      value = c(1693, 43, 2848, 17, 5055, 1281, 75, 1967, 178, 1161))
wgd_plot$wgd = factor(wgd_plot$wgd, levels = c('NONE', '1'))
wgd_plot$new_value[1:5] = (wgd_plot$value[1:5] / 9656) * 100
wgd_plot$new_value[6:10] = (wgd_plot$value[6:10] / 4662) * 100

source('Scripts/plot_theme.R')
WGD_LOY = ggplot(wgd_plot, aes(x = wgd, y = new_value, fill = category)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'loss' = '#0E3F7C',
                               'relative_loss' = '#00AEC8',
                               'gain' = '#D53833',
                               'gain_loss' = '#E3CC98'), name = '') +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_x_discrete(labels = c('none', '1x')) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(aspect.ratio = 2) +
  labs(x = '# genome doublings', y = 'Fraction of tumors')

ggsave(filename = 'Figures_original/Fraction_Y_loss_WGD.pdf', plot = WGD_LOY, device = 'pdf', width = 4)




