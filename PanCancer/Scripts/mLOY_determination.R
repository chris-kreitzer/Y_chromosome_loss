##----------------+
## MADSEQ; script to assess read-depth at Y 
## Detection of mosaic loss of the Y-chromosome
## downstream analysis; 
##----------------+

## start: 12/09/2022
## revision: 12/12/2022
## revision: 01/27/2023
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
## Downstream anlysis
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
                   X_ratio = X_ratio)
  out
}

y = lapply(unique(xx$id), function(x) loy_ratio(id = x, data = xx))
y = data.table::rbindlist(y)
colnames(y) = c('id', 'Y_Autosome_ratio', 'Y_RefDepth_ratio', 'chromosome22_ratio', 'X_Autosome_ratio')
y = na.omit(y)
write.table(x = y, file = 'Data/03_Mosaicism/SeqRatios_IMPACT.txt', sep = '\t', row.names = F, quote = F)


##----------------+
## expected ploidy:
## correction factor
##----------------+
par(mfrow = c(1,2))
plot(density(y$chromosome22_ratio),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5,
     xlim = c(0, 2))
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'DNA concentration [target chr./autosomes]', side = 1, line = 2)
mtext(text = 'Chromosome 22', side = 3, line = 0.5, adj = 0)

plot(density(y$Y_Autosome_ratio),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5,
     xlim = c(0, 1))
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'DNA concentration [target chr./autosomes]', side = 1, line = 2)
mtext(text = 'Chromosome Y', side = 3, line = 0.5, adj = 0)
# save: (device size): 8.61 x 3.82


##----------------+
## correction factor
## - for autosomes (chr22)
## - for the Y-chromosome
## - for the X-chromosome
##----------------+
y$E22 = y$chromosome22_ratio * 2
density.chromo22 = density(y$E22, na.rm = T)
density.max22 = density.chromo22$x[which.max(density.chromo22$y)]
correction_factor22 = density.max22 - 2

y$EY = y$Y_Autosome_ratio * 2
density.chromoY = density(y$EY, na.rm = T)
density.maxY = density.chromoY$x[which.max(density.chromoY$y)]
correction_factorY = density.maxY - 1

y$EX = y$X_Autosome_ratio * 2
density.chromoX = density(y$EX, na.rm = T)
density.maxX = density.chromoX$x[which.max(density.chromoX$y)]
correction_factorX = density.maxX - 1

y$O22 = abs(y$E22 - correction_factor22)
y$OY = abs(y$EY - correction_factorY)
y$OX = abs(y$EX - correction_factorX)


## plot: observed vs corrected
par(mfrow = c(2,2))
plot(density(y$chromosome22_ratio),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5,
     xlim = c(0, 2))
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'DNA concentration [target chr./autosomes]', side = 1, line = 2)
mtext(text = 'Chromosome 22: observed', side = 3, line = 0.5, adj = 0)

plot(density(y$Y_Autosome_ratio),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5,
     xlim = c(0, 1))
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'DNA concentration [target chr./autosomes]', side = 1, line = 2)
mtext(text = 'Chromosome Y: observed', side = 3, line = 0.5, adj = 0)

plot(density(y$O22),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5,
     xlim = c(0, 4))
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'Ploidy', side = 1, line = 2)
mtext(text = 'Chromosome 22: corrected', side = 3, line = 0.5, adj = 0)

plot(density(y$OY),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5,
     xlim = c(0, 2))
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'Ploidy', side = 1, line = 2)
mtext(text = 'Chromosome Y: corrected', side = 3, line = 0.5, adj = 0)

# save: (device size): 8.61 x 5



##----------------+
## determine the cut-off
## through percentiles; 
## and the surrounding 95% CI
## (lower 2.5% percentile)
##----------------+
lower_cutoff = quantileCI(x = y$OY, probs = 0.025, conf.level = .95)
lower = lower_cutoff[1]

##-------
## Plot
y$seq = seq(1, nrow(y), 1)
jitter = ggplot(y, aes(x = OY, y = seq)) +
  geom_jitter(shape = 17, size = 0.5) +
  geom_vline(xintercept = lower, col = '#a22231', linewidth = 0.75, linetype = 'dashed') +
  #geom_vline(xintercept = quantile(y$OY, probs = 0.99)[[1]], col = 'red', linewidth = 0.75) +
  theme_std(base_size = 14, base_line_size = 0.2) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 1.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(x = c(0,0,0,0), units = 'mm')) +
  labs(y = 'Individuals', x = 'Ploidy')

histo = ggplot(y, aes(x = OY, y = ..density..)) +
  geom_histogram(bins = 400, col = 'black', fill = 'black') +
  geom_vline(xintercept = lower, col = '#a22231', linewidth = 0.75, linetype = 'dashed') +
  #geom_vline(xintercept = quantile(y$OY, probs = 0.99)[[1]], col = 'red', linewidth = 0.75) +
  theme_std(base_size = 14, base_line_size = 0.2) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 1.5),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(x = c(0, 0, 0, 0), units = 'mm'))

plot_grid(histo, jitter, nrow = 2, rel_heights = c(1,3), align = 'hv')

# save: (custom device)



##----------------+
## cut-off threshold
##----------------+
lower = quantileCI(x = y$OY, probs = 0.025, conf.level = .95)[1]
y$mLOY = ifelse(y$OY < lower, 'mLOY', 'no_mLOY')
write.table(x = y, file = 'Data/03_Mosaicism/IMPACT_mLOY_summary.txt', sep = '\t', row.names = F, quote = F)




##----------------+
## Association studies;
## Age
##----------------+ 
IMPACT = read.csv('Data/00_CohortData/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')
IMPACT = IMPACT[!IMPACT$Age_Sequencing %in% c('0', '1', '2', '3', '4', '5', '6', '7'), ]
IMPACT = IMPACT[!is.na(IMPACT$Age_Sequencing), ]
IMPACT$Age_Sequencing = as.integer(as.character(IMPACT$Age_Sequencing))

LOY_Age = merge(y, IMPACT[, c('SAMPLE_ID', 'Age_Sequencing')], by.x = 'id', by.y = 'SAMPLE_ID', all.x = T)


## Visualization
intercept = summary(lm(LOY_Age$OY ~ LOY_Age$Age_Sequencing))[[4]][1]
slope = summary(lm(LOY_Age$OY ~ LOY_Age$Age_Sequencing))[[4]][2]
p.value = summary(lm(LOY_Age$OY ~ LOY_Age$Age_Sequencing))[[4]][8]

Y_chromosome = ggplot(LOY_Age, aes(x = Age_Sequencing, y = OY)) +
  geom_jitter(size = 0.2) +
  scale_y_continuous(limits = c(0.5, 1.5),
                     breaks = c(0.5, 1, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.05),
                     breaks = seq(10, 90, 10)) +
  geom_smooth(method = 'lm') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, 
                                    linewidth = 2),
        axis.text = element_text(size = 12, 
                                 color = 'black')) +
  annotate(geom = 'text',
           x = 20,
           y = 1.40,
           label = paste0('y = ', round(intercept, 3), round(slope, 3), '\np: ', round(p.value, 3)),
           size = 5.0) +
  labs(x = 'Age [reported at sequencing]', y = 'Chromosome Y ploidy', title = 'MSK-IMPACT: Ploidy decrease with age')

Y_chromosome


## X-chromosome
intercept = summary(lm(LOY_Age$OX ~ LOY_Age$Age_Sequencing))[[4]][1]
slope = summary(lm(LOY_Age$OX ~ LOY_Age$Age_Sequencing))[[4]][2]
p.value = summary(lm(LOY_Age$OX ~ LOY_Age$Age_Sequencing))[[4]][8]

X_chromosome = ggplot(LOY_Age, aes(x = Age_Sequencing, y = OX)) +
  geom_jitter(size = 0.2) +
  scale_y_continuous(limits = c(0.5, 1.5),
                     breaks = c(0.5, 1, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.05),
                     breaks = seq(10, 90, 10)) +
  geom_smooth(method = 'lm') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, 
                                    linewidth = 2),
        axis.text = element_text(size = 12, 
                                 color = 'black')) +
  annotate(geom = 'text',
           x = 20,
           y = 1.40,
           label = paste0('y = ', round(intercept, 3), round(slope, 3), '\np: ', round(p.value, 3)),
           size = 5.0) +
  labs(x = 'Age [reported at sequencing]', y = 'Chromosome Y ploidy', title = 'MSK-IMPACT: Ploidy decrease with age')

X_chromosome


#' out
