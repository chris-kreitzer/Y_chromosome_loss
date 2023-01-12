##----------------+
## Genome-scores; and
## correlations with 
## LOY
##----------------+
##
## start: 11/01/2022
## revision: 01/10/2023
## revision: 01/12/2023
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
## insights into plots
## above
##----------------+
arm_subset = arm_changes[which(arm_changes$CANCER_TYPE %in% c('Bladder Cancer', 'Bone Cancer', 'Melanoma', 'Pancreatic Cancer')), ]

ggplot(arm_subset, aes(x = CANCER_TYPE, y = losses_n)) +
  geom_jitter(width = 0.2) +
  geom_boxplot()
  

t.test(arm_subset$fraction_cna[which(arm_subset$CANCER_TYPE == 'Bone Cancer')], 
       arm_subset$fraction_cna[which(arm_subset$CANCER_TYPE == 'Bladder Cancer')])

##' Question is:
##' Why is the LOY rate in Pancreatic and Melanoma tumors much higher than in 
##' Bone or Bladder cancers??
##' look at the mutations and maybe the CNA alterations?



genome_scores_balanced = genome_scores[which(genome_scores$CancerType %in% c('Mesothelioma', 'Bladder Cancer', 'Glioma',
                                                                        'Skin Cancer, Non-Melanoma', 'Prostate Cancer',
                                                                        'Gastrointestinal Stromal Tumor', 'Soft Tissue Sarcoma',
                                                                        'Small Cell Lung Cancer')), ]

FGA_phase1 = ggplot(genome_scores_balanced, aes(x = fraction_LOY, y = mean_fga)) +
  geom_point(shape = 20) +
  stat_smooth(method = lm, fullrange = FALSE, alpha = 0.1, span = 10, se = F) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0.25 , 1)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_text_repel(aes(label = CancerType), size = 4) +
  theme_std(base_size = 14, base_line_size = 1) + 
  theme(aspect.ratio = 1) +
  geom_text(x = 0.12, y = 0.94, label = 'cor: 0.74', hjust = 0, size = 5) +
  labs(y = 'Average FGA')
  
FGA_phase1

##----------------+
genome_scores_negative = genome_scores[which(genome_scores$CancerType %in% c('Head and Neck Cancer', 'Hepatobiliary Cancer',
                                                                             'Melanoma', 'Renal Cell Carcinoma', 'Cancer of Unknown Primary')), ]

FGA_phase2 = ggplot(genome_scores_negative, aes(x = fraction_LOY, y = mean_fga)) +
  geom_point(shape = 20) +
  stat_smooth(method = lm, fullrange = FALSE, alpha = 0.1, span = 10, se = F) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0.25 , 1)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_text_repel(aes(label = CancerType), size = 4) +
  theme_std(base_size = 14, base_line_size = 1) + 
  theme(aspect.ratio = 1) +
  geom_text(x = 0.4, y = 0.94, label = 'cor: -0.76', hjust = 0, size = 5) +
  labs(y = 'Average FGA')

FGA_phase1 / FGA_phase2
