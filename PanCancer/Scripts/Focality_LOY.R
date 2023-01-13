##----------------+
## Focality of chrom. Y
## alterations;
##----------------+
##
## AND
## 
## ##----------------+
## run GISTIC on n=546 samples
## where we have multiple segments
## from the Y-chromosome;
## determine whether there is any
## focal alteration above the background?
##----------------+
##
## start: 11/07/2022
## revision: 01/13/2023
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



##----------------+
## GISTIC on double segments
##----------------+
##-------
## Modify GISTIC
## output; to be plotted
##-------
gene.lookup = read.csv('Data/01_Coverage_Depth/BAM_quality_out.txt', sep = '\t')
gene.lookup = unique(gene.lookup$gene)

all_output_files = list.files(path = 'Data/04_Loss/GISTIC_results/', full.names = T)
all_data_genes = read.csv(file = grep('.thresholded.by_gene.*', all_output_files, value = T), sep = '\t')

## fetch genes from chromosome Y
sample_gene.lookup = all_data_genes[all_data_genes$Gene.Symbol %in% as.character(gene.lookup),, drop = F]
sample_gene.lookup = as.data.frame(t(sample_gene.lookup))
colnames(sample_gene.lookup) = sample_gene.lookup[1, ]
sample_gene.lookup = droplevels(sample_gene.lookup[-c(1,2,3),, drop = F])
sample_gene.lookup$ID = row.names(sample_gene.lookup)
rownames(sample_gene.lookup) = NULL
sample_gene.lookup[1:12] = lapply(sample_gene.lookup[1:12], str_trim) 

all_genes = data.frame()
for(i in (1:(length(sample_gene.lookup) - 1))){
  print(i)
  n = dim(sample_gene.lookup)[1]
  gene = colnames(sample_gene.lookup)[i]
  crosstab = as.data.frame(table(sample_gene.lookup[,i]))
  crosstab$gene = gene
  crosstab$freq = crosstab$Freq / n
  colnames(crosstab) = c('Category', 'Frequency', 'gene', 'rel_freq')
  all_genes = rbind(all_genes, crosstab)
}

xx = all_genes[which(all_genes$Category == '-1'),, drop = F ]
xx = xx[order(xx$rel_freq), ]
fn = factor(xx$gene, levels = xx$gene)
all_genes$gene = factor(all_genes$gene, levels = levels(fn))

all_genes$Category = as.character(all_genes$Category)
all_genes$Category = ifelse(all_genes$Category == '2', '1',
                            ifelse(all_genes$Category == '-2', '-1', all_genes$Category))
all_genes$Category = factor(all_genes$Category, levels = rev(c('-1', '0', '1')))


Genes_GISTIC = ggplot(all_genes,
                      aes(x = gene,
                          y = rel_freq,
                          fill = Category, 
                          label = paste0((round(rel_freq* 100, 2)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  coord_flip() +
  scale_fill_manual(values = c('0' = '#D7D8DA',
                               '-1' = '#00AEC8',
                               '1' = '#E3CC98'),
                    name = '') +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 1, 0.2),
                     sec.axis = dup_axis()) +
  theme_std(base_size = 14, base_line_size = 1) +
  labs(x = '', y = 'Fraction')


ggsave_golden(filename = 'Figures_original/Genes_GISTIC_alterations.pdf', plot = Genes_GISTIC, width = 6)

## all_lesions.conf.XX
all_lesions = read.csv(file = grep('.lesions.conf.*', all_output_files, value = T), sep = '\t')
all_lesions = all_lesions[grep('0.*', all_lesions$Amplitude.Threshold),, drop = F]
all_lesions_subset1 = all_lesions[, c('Unique.Name', 'Descriptor', 'Peak.Limits', 'q.values')]
all_lesions_subset2 = all_lesions[, c(10:ncol(all_lesions))]
all_lesions_subset2$X = NULL

for(i in 1:ncol(all_lesions_subset2)){
  all_lesions_subset2[, i] = ifelse(all_lesions_subset2[, i] == 2, 1,
                                    ifelse(all_lesions_subset2[, i] == 1, 1, 0))
}
sample.summary = data.frame(Sample.summary = rowSums(all_lesions_subset2))
all_lesions_out = cbind(all_lesions_subset1, sample.summary)
all_lesions_out$Peak.Limits = gsub(pattern = '.p.*', replacement = '', all_lesions_out$Peak.Limits)

## work on Gistic output
gistic = read.csv(file = grep('.cores.gistic.*', all_output_files, value = T), sep = '\t')


##-------
## Visualization;
##-------
plot_Gscores_Gistic = ggplot() +
  geom_path(data = subset(gistic, gistic$Type == 'Amp'), aes(x = Start, y = G.score), lineend = 'butt', linejoin = 'round', 
            linemitre = 50, linetype = 'solid', linewidth = 1.5, col = 'red') +
  geom_path(data = subset(gistic, gistic$Type == 'Del'), aes(x = Start, y = G.score), lineend = 'butt', linejoin = 'round', 
            linemitre = 50, linetype = 'solid', linewidth = 1.5, col = 'blue') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.8)) +
  theme(aspect.ratio = 3, axis.text.x = element_blank()) +
  labs(x = '', y = 'G score') +
  
  facet_wrap(~Chromosome, scales = 'free_x', nrow = 1) +
  
  theme(panel.spacing = unit(0.05, "lines"),
        axis.ticks = element_blank(),
        strip.placement = 'inside', 
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  panel_border(color = 'black')

plot_Gscores_Gistic
ggsave_golden(filename = 'Figures_original/GISTIC_scores_Y.pdf', plot = plot_Gscores_Gistic, width = 16)

#' out