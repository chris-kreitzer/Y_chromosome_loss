##----------------+
## MADSEQ; script to assess read-depth at Y 
## Detection of mosaic loss of the Y-chromosome
## downstream analysis; 
## Obvious signs of bi-modality in the data
## distribution; investigate what's happening there
##----------------+

## start: 12/09/2022
## revision: 12/12/2022
## revision: 01/27/2023
## revision: 01/30/2023
## revision: 03/01/2023
## 
## chris-kreitzer

clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')
library(cowplot)
library(jmuOutlier)


##----------------+
## which samples where
## not considered for
## mLOY calculation
##----------------+
IMPACT = read.csv('Data/00_CohortData/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')
IMPACT_mLOY = IMPACT[!is.na(IMPACT$Normal_bam), ]
mLOY_files = list.files('Data/03_Mosaicism/NormalizedDepth/', full.names = T)
mLOY_files_short = substr(x = basename(mLOY_files), start = 1, stop = 17)
length(intersect(mLOY_files_short, IMPACT_mLOY$SAMPLE_ID))
length(intersect(mLOY_files_short, IMPACT_mLOY$SAMPLE_ID)) == length(mLOY_files)


## TODO: number may be different
## because lymphoid and myeloid cancers
## were excluded; old masterfile used 
## for mLOY calculation



##----------------+
## read the data and 
## calculate read-depth 
## ratio
##----------------+
mLOY_files = list.files('Data/03_Mosaicism/NormalizedDepth/', full.names = T)
read_mLOY = function(data){
  print(data)
  name = substr(basename(data), start = 1, stop = 17)
  file = read.csv(data, sep = '\t')
  file$width = NULL
  file$strand = NULL
  file$id = name
  file
}

x = lapply(unique(mLOY_files), function(x) read_mLOY(data = x))
x = data.table::rbindlist(x)
xx = as.data.frame(x)


##----------------+
mloy = merge(xx, cohort[,c('SAMPLE_ID', 'GENE_PANEL')], by.x = 'id', by.y = 'SAMPLE_ID', all.x = T)


##----------------+
## Downstream analysis
##----------------+
loy_ratio = function(id, data){
  name = id
  data = data[which(data$id == id), ]
  print(name)
  Y = median(data$normed_depth[which(data$seqnames == 'chrY')])
  X = median(data$normed_depth[which(data$seqnames == 'chrX')])
  A = median(data$normed_depth[!data$seqnames %in% c('chrX', 'chrY')])
  c22 = median(data$normed_depth[which(data$seqnames == 'chr22')])
  c22_ratio = c22 / A
  ratio = Y / A
  ratio_norm = Y / unique(data$ref_depth)
  X_ratio = X / A
  out = data.frame(id = name,
                   ratio = ratio,
                   ratio_norm = ratio_norm,
                   c22 = c22_ratio,
                   X_ratio = X_ratio,
                   ChromY = Y,
                   Chrom22 = c22)
  out
}

y = lapply(unique(mloy$id), function(x) loy_ratio(id = x, data = mloy))
y = data.table::rbindlist(y)
colnames(y) = c('id', 'Y_Autosome_ratio', 'Y_RefDepth_ratio', 'chromosome22_ratio', 'X_Autosome_ratio', 'ChromY', 'Chrom22')
mloy_out = na.omit(y)
write.table(x = mloy_out, file = 'Data/03_Mosaicism/SeqRatios_NEW_IMPACT.txt', sep = '\t', row.names = F, quote = F)


##----------------+
## Y DNA concentration
## based on gene panel;
## absolute and relative concentration
##----------------+
mloy_out$panel = substr(mloy_out$id, start = 15, stop = 17)

##-- absolute DNA concentration
dev.off()
par(mfrow = c(1,2))
boxplot(mloy_out$ChromY ~ mloy_out$panel,
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5)
box(lwd = 2)
mtext(text = 'Absolute DNA concentration', side = 2, line = 1)
mtext(text = 'GenePanel', side = 1, line = 2)
mtext(text = 'Chromosome Y', side = 3, line = 0.5, adj = 0)

##-- relative DNA concentration
boxplot(mloy_out$Y_Autosome_ratio ~ mloy_out$panel,
        yaxt = 'n',
        ylab = '',
        xlab = '',
        main = '',
        lwd = 1.5)
box(lwd = 2)
mtext(text = 'relative DNA concentration [target chr./autosomes]', side = 2, line = 1)
mtext(text = 'GenePanel', side = 1, line = 2)
mtext(text = 'Chromosome Y', side = 3, line = 0.5, adj = 0)



##----------------+
## Determine correction factor
## E(Y) = 1, based on experimental 
## variation and determine 95th percentile
## cutoff for every gene panel individually
## 
## Cut-off through percentiles; 
## and the surrounding 95% CI
## (lower 5th percentile)
##----------------+

mosaic_out = data.frame()
for(i in unique(mloy_out$panel)){
  data = mloy_out[which(mloy_out$panel == i), ]
  data$expected_Y = data$Y_Autosome_ratio * 2
  density.chromoY = density(data$expected_Y, na.rm = T)
  density.maxY = density.chromoY$x[which.max(density.chromoY$y)]
  correction_factorY = density.maxY - 1
  data$observed_Y = abs(data$expected_Y - correction_factorY)
  lower_cutoff = quantileCI(x = data$observed_Y, probs = 0.05, conf.level = .95)
  lower = lower_cutoff[1]
  data$cutoff = lower
  data$mLOY = ifelse(data$observed_Y < lower, 'mloy', 'wt')
  mosaic_out = rbind(mosaic_out, data)
}


##-------
## Plot; Visualization
##-------
plot_list = list()
for(i in unique(mosaic_out$panel)){
  y = mosaic_out[which(mosaic_out$panel == i), ]
  
  y$seq = seq(1, nrow(y), 1)
  jitter = ggplot(y, aes(x = observed_Y, y = seq)) +
    geom_jitter(shape = 17, size = 1) +
    scale_x_continuous(limits = c(0.5, 1.5),
                       breaks = c(0.5, 1, 1.5)) +
    geom_vline(xintercept = lower, col = '#a22231', linewidth = 0.75, linetype = 'dashed') +
    theme_std(base_size = 14, base_line_size = 0.2) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, linewidth = 1.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(x = c(0,0,0,0), units = 'mm')) +
    labs(y = 'Individuals', x = 'Ploidy')
  
  histo = ggplot(y, aes(x = observed_Y, y = ..density..)) +
    geom_histogram(bins = 400, col = 'black', fill = 'black') +
    geom_vline(xintercept = lower, col = '#a22231', linewidth = 0.75, linetype = 'dashed') +
    scale_x_continuous(limits = c(0.5, 1.5),
                       breaks = c(0.5, 1, 1.5)) +
    theme_std(base_size = 14, base_line_size = 0.2) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, linewidth = 1.5),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(x = c(0, 0, 0, 0), units = 'mm')) +
    labs(title = paste0('GenePanel: ', i, ' (n=', dim(y)[[1]], ')'))
  
  plot_combined = plot_grid(histo, jitter, nrow = 2, rel_heights = c(1,3), align = 'hv')
  plot_list[[i]] = plot_combined
}


mosaic_out_plot = plot_list[[1]] + plot_list[[2]] + plot_list[[4]] + plot_list[[3]] + plot_layout(nrow = 1, ncol = 4)
ggsave_golden(filename = 'Figures_original/Mosaic_out_genePanel.pdf', plot = mosaic_out_plot, width = 15)

write.table(x = mosaic_out, file = 'Data/03_Mosaicism/IMPACT_mLOY_summary.txt', sep = '\t', row.names = F, quote = F)


##----------------+
## Association studies;
## mosaic loss of Y and Age
## 03/02/2023
##----------------+ 
IMPACT = readRDS('Data/00_CohortData/Cohort_071322.rds')
IMPACT$Age_Sequencing = as.integer(as.character(IMPACT$Age_Sequencing))
IMPACT = IMPACT[!is.na(IMPACT$Age_Sequencing), ]
IMPACT = IMPACT[!is.na(IMPACT$observed_Y), ]


##-- simple comparison
mLOY_Age_plot = ggplot(IMPACT, 
                       aes(x = mLOY, y = Age_Sequencing)) +
  geom_violin(width = .75) + 
  geom_quasirandom(alpha = 0.1, 
                   width = 0.2) +
  stat_summary(fun.y = "mean",
               geom = "crossbar", 
               mapping=aes(ymin = ..y.., ymax = ..y..), 
               width = 0.5, 
               position = position_dodge(),
               show.legend = FALSE,
               color = 'darkgrey') +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(18, 92)) +
  scale_x_discrete(labels = c(paste0('wildtype\n(n= ', length(IMPACT$SAMPLE_ID[which(IMPACT$mLOY == 'no')]), ')'),
                              paste0('mosaic loss chromosome Y\n(n= ', length(IMPACT$SAMPLE_ID[which(IMPACT$mLOY == 'yes')]), ')'))) +
  theme_std(base_size = 14) +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 2),
    plot.margin = unit(x = c(0, 0, 0, 0), units = 'mm'),
    aspect.ratio = 1,
    axis.text = element_text(size = 16)) +
  stat_compare_means(label = 'p.format', label.y = 88, label.x = 1.3) +
  labs(x = '', y = 'Age [reported at sequencing]')

mLOY_Age_plot
ggsave_golden(filename = 'Figures_original/mLRR-Y_Age_correlation.pdf', plot = mLOY_Age_plot, width = 9.9)


## Visualization
intercept = summary(lm(IMPACT$observed_Y ~ IMPACT$Age_Sequencing))[[4]][1]
slope = summary(lm(IMPACT$observed_Y ~ IMPACT$Age_Sequencing))[[4]][2]
p.value = summary(lm(IMPACT$observed_Y ~ IMPACT$Age_Sequencing))[[4]][8]

Y_chromosome = ggplot(IMPACT, 
                      aes(x = Age_Sequencing, 
                          y = observed_Y)) +
  geom_jitter(shape = 17, size = 0.5) +
  scale_y_continuous(limits = c(0.5, 1.5),
                     breaks = c(0.5, 1, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(18, 90), 
                     breaks = seq(20, 90, 10)) +
  geom_smooth(method = 'lm') +
  theme_std(base_size = 14) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, 
                                    linewidth = 2)) +
  annotate(geom = 'text',
           x = 25,
           y = 1.40,
           label = paste0('y = ', round(intercept), round(slope, 4), '\np: ', round(p.value, 4)),
           size = 5.0) +
  labs(x = 'Age [reported at sequencing]', y = 'Chromosome Y ploidy')

Y_chromosome
ggsave_golden(filename = 'Figures_original/mLOY_Age_correlation.pdf', plot = Y_chromosome, width = 6)



##-- depreciated; misc
##-----------------
## correction factor
## - for autosomes (chr22)
## - for the Y-chromosome
## - for the X-chromosome
##----------------+
# y = read.csv('Data/03_Mosaicism/SeqRatios_IMPACT.txt', sep = '\t')
# y$E22 = y$chromosome22_ratio * 2
# density.chromo22 = density(y$E22, na.rm = T)
# density.max22 = density.chromo22$x[which.max(density.chromo22$y)]
# correction_factor22 = density.max22 - 2
# 
# y$EY = y$Y_Autosome_ratio * 2
# density.chromoY = density(y$EY, na.rm = T)
# density.maxY = density.chromoY$x[which.max(density.chromoY$y)]
# correction_factorY = density.maxY - 1
# 
# y$EX = y$X_Autosome_ratio * 2
# density.chromoX = density(y$EX, na.rm = T)
# density.maxX = density.chromoX$x[which.max(density.chromoX$y)]
# correction_factorX = density.maxX - 1
# 
# y$O22 = abs(y$E22 - correction_factor22)
# y$OY = abs(y$EY - correction_factorY)
# y$OX = abs(y$EX - correction_factorX)
# 
# 
# ## plot: observed vs corrected
# par(mfrow = c(2,2))
# plot(density(y$chromosome22_ratio),
#      yaxt = 'n',
#      ylab = '',
#      xlab = '',
#      main = '',
#      lwd = 1.5,
#      xlim = c(0, 2))
# box(lwd = 2)
# mtext(text = 'Density', side = 2, line = 1)
# mtext(text = 'DNA concentration [target chr./autosomes]', side = 1, line = 2)
# mtext(text = 'Chromosome 22: observed', side = 3, line = 0.5, adj = 0)
# 
# plot(density(y$Y_Autosome_ratio),
#      yaxt = 'n',
#      ylab = '',
#      xlab = '',
#      main = '',
#      lwd = 1.5,
#      xlim = c(0, 1))
# box(lwd = 2)
# mtext(text = 'Density', side = 2, line = 1)
# mtext(text = 'DNA concentration [target chr./autosomes]', side = 1, line = 2)
# mtext(text = 'Chromosome Y: observed', side = 3, line = 0.5, adj = 0)
# 
# plot(density(y$O22),
#      yaxt = 'n',
#      ylab = '',
#      xlab = '',
#      main = '',
#      lwd = 1.5,
#      xlim = c(0, 4))
# box(lwd = 2)
# mtext(text = 'Density', side = 2, line = 1)
# mtext(text = 'Ploidy', side = 1, line = 2)
# mtext(text = 'Chromosome 22: corrected', side = 3, line = 0.5, adj = 0)
# 
# plot(density(y$OY),
#      yaxt = 'n',
#      ylab = '',
#      xlab = '',
#      main = '',
#      lwd = 1.5,
#      xlim = c(0, 2))
# box(lwd = 2)
# mtext(text = 'Density', side = 2, line = 1)
# mtext(text = 'Ploidy', side = 1, line = 2)
# mtext(text = 'Chromosome Y: corrected', side = 3, line = 0.5, adj = 0)


##----------------+
## X-chromosome as
## an test variable; should not decay with age
##----------------+
# cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
# mloy = read.csv('Data/03_Mosaicism/IMPACT_mLOY_summary.txt', sep = '\t')
# x = merge(cohort, mloy[,c('id', 'OX')], by.x = 'SAMPLE_ID', by.y = 'id', all.x = T)
# x = x[!is.na(x$Age_Sequencing), ]
# x$Age_Sequencing = as.integer(as.character(x$Age_Sequencing))
# 
# intercept = summary(lm(x$OX ~ x$Age_Sequencing))[[4]][1]
# slope = summary(lm(x$OX ~ x$Age_Sequencing))[[4]][2]
# p.value = summary(lm(x$OX ~ x$Age_Sequencing))[[4]][8]
# 
# X_chromosome = ggplot(x, aes(x = Age_Sequencing, y = OX)) +
#   geom_jitter(shape = 17, size = 0.5) +
#   scale_y_continuous(limits = c(0.5, 1.5),
#                      breaks = c(0.5, 1, 1.5)) +
#   scale_x_continuous(expand = c(0.01, 0.05),
#                      breaks = seq(10, 90, 10)) +
#   geom_smooth(method = 'lm') +
#   theme(panel.background = element_blank(),
#         panel.border = element_rect(fill = NA, 
#                                     linewidth = 2),
#         axis.text = element_text(size = 12, 
#                                  color = 'black')) +
#   annotate(geom = 'text',
#            x = 20,
#            y = 1.40,
#            label = paste0('y = ', round(intercept, 3), round(slope, 3), '\np: ', round(p.value, 3)),
#            size = 5.0) +
#   labs(x = 'Age [reported at sequencing]', y = 'Chromosome Y ploidy', title = 'MSK-IMPACT: Ploidy decrease with age')
# 
# X_chromosome


#' out