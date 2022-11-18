
clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')

library(ggrepel)
library(cowplot)
source('Scripts/plot_theme.R')
data = readRDS('Data/signedOut/Cohort_07132022.rds')
Y_calls = data$IMPACT_Y_classification_final
Y_calls = merge(Y_calls, data$IMPACT_binaryY_call[,c('sample.id', 'purity')],
                by.x = 'sample', by.y = 'sample.id', all.x = T)
Y_calls$purity[which(Y_calls$purity == 0)] = 0.0001

purity_bin = data.frame()
for(i in seq(0, 0.8, 0.2)){
  #print(c(i, i+0.2))
  data_sub = Y_calls[which(Y_calls$purity > i & Y_calls$purity <= i + 0.2), ]
  n_loss = length(data_sub$sample[which(data_sub$classification %in% c('loss', 'relative_loss', 'gain_loss'))])
  n_loss_rel = n_loss / nrow(data_sub)
  out = data.frame(purity_bin = i,
                   total = nrow(data_sub),
                   n_loss_rel = n_loss_rel)
  purity_bin = rbind(purity_bin, out)
}

barplot(height = c(43.4, 39.0, 37.3, 32.5, 26.0), 
        names.arg = c('[0-0.2]', '[0.2-0.4]', '[0.4-0.6]', '[0.6-0.8]', '[0.8-1]'), 
        space = 0.1, col = 'white', ylim = c(0, 100), las = 1, xlab = 'purity-bin', 
        ylab = 'Samples [%]', main = 'Chromosome Y-loss in varying purity groups')

purity_bin


##----------------+
## Getting FGA and 
## number of chromosome 
## arm alterations (AS score)
##----------------+
clean()
gc()
.rs.restartR()
setwd(dir = '~/Documents/MSKCC/10_MasterThesis/')
full_cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
cohort_pass = full_cohort$IMPACT_Y_classification_final
#arms = read.csv('Data/04_Loss/IMPACT_arm_change_out.txt', sep = '\t')
segs = read.csv('Data/04_Loss/IMPACT_copynumber_out.txt', sep = '\t')
segs_pass = segs[which(segs$ID %in% cohort_pass$sample), ]
copyNumber = facetsSuite:::copy_number_states
cn_state_loss = copyNumber$call[which(copyNumber$numeric_call %in% c(-2, -1))]

arm_changes_full_cohort = data.frame()
for(i in unique(segs_pass$ID)){
  print(i)
  try({
    id = i
    purity = unique(segs_pass$purity[which(segs_pass$ID == i)])
    ploidy = unique(segs_pass$ploidy[which(segs_pass$ID == i)])
    arm_change = facetsSuite::arm_level_changes(segs = segs_pass[which(segs_pass$ID == i), ],
                                                ploidy = ploidy,
                                                genome = 'hg19',
                                                algorithm = 'em')
    arm_change_df = arm_change$full_output
    Aneuploidy = filter(arm_change_df, cn_state != 'DIPLOID')
    loss = filter(arm_change_df, cn_state %in% cn_state_loss)
    out = data.frame(id = id,
                     purity = purity,
                     ploidy = ploidy,
                     genome_doubled = arm_change$genome_doubled,
                     fraction_cna = arm_change$fraction_cna,
                     weighted_fraction_cna = arm_change$weighted_fraction_cna,
                     AS_score = nrow(Aneuploidy),
                     losses_n = nrow(loss))
    arm_changes_full_cohort = rbind(arm_changes_full_cohort, out)
  })
}



##----------------+
## Correlation analysis
##----------------+
clean()
gc()
.rs.restartR()
cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
top20 = read.csv('Data/04_Loss/IMPACT_Y_loss_incidences_top20.txt', sep = '\t')
Y_calls = cohort$IMPACT_Y_classification_final
arm_changes = cohort$IMPACT_ARM_level_changes
arm_changes = merge(arm_changes, cohort$IMPACT_clinicalAnnotation[,c('SAMPLE_ID', 'CANCER_TYPE', 'MUTATION_COUNT')],
                    by.x = 'id', by.y = 'SAMPLE_ID', all.x = T)
arm_changes = merge(arm_changes, Y_calls[,c('sample', 'classification')],
                    by.x = 'id', by.y = 'sample', all.x = T)

genome_scores = data.frame()
for(i in unique(arm_changes$CANCER_TYPE)){
  type = i
  median_AS = median(x = arm_changes$AS_score[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  mean_AS = mean(x = arm_changes$AS_score[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  median_loss = median(x = arm_changes$losses_n[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  mean_loss = mean(x = arm_changes$losses_n[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  median_mutation = median(x = arm_changes$MUTATION_COUNT[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  mean_mutation = mean(x = arm_changes$MUTATION_COUNT[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  median_AS = median(x = arm_changes$AS_score[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  mean_AS = mean(x = arm_changes$AS_score[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  median_fga = median(x = arm_changes$fraction_cna[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  mean_fga = mean(x = arm_changes$fraction_cna[which(arm_changes$CANCER_TYPE == i)], na.rm = T)
  n_loy = length(arm_changes$id[which(arm_changes$CANCER_TYPE == i & arm_changes$classification %in% c('loss', 'relative_loss'))])
  fraction_LOY = n_loy / length(arm_changes$id[which(arm_changes$CANCER_TYPE == i)])
  out = data.frame(CancerType = type,
                   median_AS = median_AS,
                   mean_AS = mean_AS,
                   median_loss = median_loss,
                   mean_loss = mean_loss,
                   median_mutation = median_mutation,
                   mean_mutation = mean_mutation,
                   median_fga = median_fga,
                   mean_fga = mean_fga,
                   fraction_LOY = fraction_LOY)
  genome_scores = rbind(genome_scores, out)
}

genome_scores = genome_scores[!is.na(genome_scores$CancerType), ]
genome_scores = genome_scores[which(genome_scores$CancerType %in% unique(top20$CancerType)), ]

##----------------+
## Visualization
##----------------+
anouploidy_general = ggplot(genome_scores, aes(x = fraction_LOY, y = mean_AS)) +
  geom_point(shape = 20) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(10, 35)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_text_repel(aes(label = CancerType), size = 2) +
  stat_smooth(method = loess, fullrange = FALSE, alpha = 0.1, span = 10) +
  theme_std(base_size = 14, base_line_size = 1) + 
  theme(aspect.ratio = 1) +
  geom_text(x = 0.12, y = 33, label = 'spearman: 0.6', hjust = 0, size = 6) +
  labs(y = 'Average # of chrom. arms altered')
  
loss_general = ggplot(genome_scores, aes(x = fraction_LOY, y = mean_loss)) +
  geom_point(shape = 20) +
  stat_smooth(method = loess, fullrange = FALSE, alpha = 0.1, span = 10) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(7, 35)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_text_repel(aes(label = CancerType), size = 2) +
  theme_std(base_size = 14, base_line_size = 1) + 
  theme(aspect.ratio = 1) +
  geom_text(x = 0.12, y = 33, label = 'spearman: 0.59', hjust = 0, size = 6) +
  labs(y = 'Average # of chrom. arms lost')

FGA_general = ggplot(genome_scores, aes(x = fraction_LOY, y = mean_fga)) +
  geom_point(shape = 20) +
  stat_smooth(method = loess, fullrange = FALSE, alpha = 0.1, span = 10) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0.25 , 1)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_text_repel(aes(label = CancerType), size = 2) +
  theme_std(base_size = 14, base_line_size = 1) + 
  theme(aspect.ratio = 1) +
  geom_text(x = 0.12, y = 0.94, label = 'spearman: 0.53', hjust = 0, size = 6) +
  labs(y = 'Average FGA')

# Mutation_general = ggplot(genome_scores, aes(x = fraction_LOY, y = mean_mutation)) +
#   geom_point(shape = 20) +
#   stat_smooth(method = loess, fullrange = FALSE, alpha = 0.1, span = 10) +
#   scale_y_continuous(expand = c(0, 0),
#                      limits = c(0, 35)) +
#   scale_x_continuous(expand = c(0.01, 0)) +
#   geom_text_repel(aes(label = CancerType), size = 2) +
#   theme_std(base_size = 14, base_line_size = 1) + 
#   theme(aspect.ratio = 1) +
#   geom_text(x = 0.12, y = 33, label = 'spearman: 0.35', hjust = 0, size = 6) +
#   labs(y = 'Average # of mutations')


plot_grid(anouploidy_general, loss_general, FGA_general, ncol = 3)


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
