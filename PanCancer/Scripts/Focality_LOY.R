##----------------+
## Focality of chrom. Y
## alterations;
##----------------+
##
## start: 11/07/2022
## revision: 01/10/2023
## 
## chris-kreitzer


clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/')
source('Scripts/UtilityFunctions.R')
source('Scripts/plot_theme.R')
library(patchwork)


##----------------+
## working with both
## QC TRUE and FALSE
##----------------+
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
CopyNumberData = read.csv('Data/04_Loss/010523/CopyNumberStates.txt', sep = '\t')
Categorical_Classification = read.csv('Data/04_Loss/010523/Categorial_Classification_Y.txt', sep = '\t')
CopyLOY = CopyNumberData[which(CopyNumberData$chrom == 24), ]
segments_table = table(CopyLOY$id)
mono_segments = names(segments_table[which(segments_table == 1)])
multi_segments = names(segments_table[which(segments_table > 1)])


##-------
## All;
##-------
Barchart_all = data.frame(x = 1,
                          value = c(length(unique(mono_segments)),
                                    length(unique(multi_segments))),
                          category = c('1-segment',
                                       'multi-segments'))
Barchart_all$category = factor(Barchart_all$category, levels = c('multi-segments', '1-segment'))
n = sum(Barchart_all$value)

n_segment_plot = ggplot(Barchart_all, aes(x = x, y = value/n, fill = category, label = paste0((round(value/n* 100, 1)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  geom_text(size = 3, position = position_stack(vjust = 0.5), fontface = 'bold') +
  scale_fill_manual(values = c('1-segment' = '#D7D8DA',
                               'multi-segments' = '#d0664b'),
                               name = '') +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 1, 0.2)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 2,
        legend.position = 'top') +
    labs(x = '', y = paste0('Fraction of samples\n(n=', n, ')'))


##-------
## multi;
##-------
multi_cohort = Categorical_Classification[which(Categorical_Classification$sample %in% multi_segments), ]
Barchart_multi = as.data.frame(table(multi_cohort$classification))
colnames(Barchart_multi) = c('category', 'value')
Barchart_multi$x = 1
n_new = sum(Barchart_multi$value)

multi_segment_plot = ggplot(Barchart_multi, aes(x = x, y = value/n_new, 
                           fill = category, 
                           label = paste0((round(value/n_new* 100, 1)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  geom_text(size = 3, position = position_stack(vjust = 0.5), fontface = 'bold') +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'complete_loss' = '#0E3F7C',
                               'partial_loss' = '#00AEC8',
                               'relative_loss' = '#474089',
                               'gain' = '#D53833',
                               'partial_gain' = '#E3CC98',
                               'gain_loss' = '#f9f8d5'),
                    name = '') +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 1, 0.2)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 2,
        legend.position = 'top') +
  labs(x = '', y = paste0('Fraction of samples\n(n=', n_new, ')'))



##-------
## multi and
## gain_loss;
## partial_gain;
## partial_loss
##-------
multi_gain_loss = Categorical_Classification$sample[which(Categorical_Classification$sample %in% multi_segments & Categorical_Classification$classification == 'gain_loss')]
multi_gain_loss = CopyNumberData[which(CopyNumberData$id %in% multi_gain_loss & CopyNumberData$chrom == 24), c('start', 'end', 'tcn.em', 'id', 'QC')]
Yp_start = 1
Centromere = 14000000
length_p = Centromere - Yp_start
multi_gain_loss$arm = NA

for(i in 1:nrow(multi_gain_loss)){
  print(i)
  vec = multi_gain_loss$start[i]:multi_gain_loss$end[i]
  length_vec = length(vec)
  coverage = sum(vec > Yp_start & vec < Centromere) / length_vec
  multi_gain_loss$arm[i] = ifelse(coverage > 0.75, 'Yp', 'Yq')
}

gain_loss = data.frame(x = 'gain_loss',
                       value = c(length(multi_gain_loss$id[which(multi_gain_loss$tcn.em > 0 & multi_gain_loss$arm == 'Yp')]),
                                 length(multi_gain_loss$id[which(multi_gain_loss$tcn.em > 0 & multi_gain_loss$arm == 'Yq')])),
                       category = c('Yp', 'Yq'))

##-------
## partial_gain
##-------
partial_gain = Categorical_Classification$sample[which(Categorical_Classification$sample %in% multi_segments & Categorical_Classification$classification == 'partial_gain')]
partial_gain = CopyNumberData[which(CopyNumberData$id %in% partial_gain & CopyNumberData$chrom == 24), c('start', 'end', 'tcn.em', 'id', 'QC')]
partial_gain = merge(partial_gain, Categorical_Classification[, c('sample', 'Y_expected')], by.x = 'id', by.y = 'sample', all.x = T)

partial_gain$arm = NA
for(i in 1:nrow(partial_gain)){
  print(i)
  vec = partial_gain$start[i]:partial_gain$end[i]
  length_vec = length(vec)
  coverage = sum(vec > Yp_start & vec < Centromere) / length_vec
  partial_gain$arm[i] = ifelse(coverage > 0.75, 'Yp', 'Yq')
}

partial_gain_df = data.frame(x = 'partial_gain',
                          value = c(length(partial_gain$id[which(partial_gain$tcn.em > partial_gain$Y_expected & partial_gain$arm == 'Yp')]),
                                 length(partial_gain$id[which(partial_gain$tcn.em > partial_gain$Y_expected & partial_gain$arm == 'Yq')])),
                       category = c('Yp', 'Yq'))


##-------
## partial loss
##-------
partial_loss = Categorical_Classification$sample[which(Categorical_Classification$sample %in% multi_segments & Categorical_Classification$classification == 'partial_loss')]
partial_loss = CopyNumberData[which(CopyNumberData$id %in% partial_loss & CopyNumberData$chrom == 24), c('start', 'end', 'tcn.em', 'id', 'QC')]
partial_loss = merge(partial_loss, Categorical_Classification[, c('sample', 'Y_expected')], by.x = 'id', by.y = 'sample', all.x = T)

partial_loss$arm = NA
for(i in 1:nrow(partial_loss)){
  print(i)
  vec = partial_loss$start[i]:partial_loss$end[i]
  length_vec = length(vec)
  coverage = sum(vec > Yp_start & vec < Centromere) / length_vec
  partial_loss$arm[i] = ifelse(coverage > 0.75, 'Yp', 'Yq')
}

partial_loss_df = data.frame(x = 'partial_loss',
                             value = c(length(partial_loss$id[which(partial_loss$tcn.em < partial_loss$Y_expected & partial_loss$arm == 'Yp')]),
                                       length(partial_loss$id[which(partial_loss$tcn.em < partial_loss$Y_expected & partial_loss$arm == 'Yq')])),
                             category = c('Yp', 'Yq'))


Focality_out = rbind(gain_loss, partial_gain_df, partial_loss_df)


##----------------+
## Visualization
##----------------+
n_new = length(unique(multi_segments))
Focality_plot = ggplot(Focality_out, aes(x = x, 
                         y = value/n_new, 
                         fill = category)) +
                         #label = paste0((round(value/n_new* 100, 1)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  #geom_hline(yintercept = seq(0, 1, 0.2), linewidth = 0.25, color = 'white', linetype = 'dashed') +
  scale_fill_manual(values = c('Yp' = '#D7D8DA',
                               'Yq' = '#d0664b'),
                    name = '') +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 1, 0.2)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(aspect.ratio = 2,
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = paste0('Fraction of samples\n(n=', n_new, ')'))


##-------
## combine + out
##-------
plot_grid = cowplot::plot_grid(plotlist = list(n_segment_plot, multi_segment_plot, Focality_plot), nrow = 1, ncol = 3,
                   align = 'h')
ggsave_golden(filename = 'Figures_original/Focality_LOY.pdf', plot = plot_grid, width = 14)

#' out