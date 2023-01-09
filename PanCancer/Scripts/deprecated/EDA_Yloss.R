clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/')
source('Scripts/plot_theme.R')

data = readRDS('Data/signedOut/Cohort_07132022.rds')


CNA = data$IMPACT_Y_classification_final

n = length(unique(CNA$sample))
prop_out = data.frame()
for(i in unique(CNA$classification)){
  fraction = length(CNA$sample[which(CNA$classification == i)]) / n
  out = data.frame(category = i,
                   fraction = fraction)
  prop_out = rbind(prop_out, out)
}

prop_out$arrange = 1
prop_out$category = factor(prop_out$category, levels = rev(c('wt', 'loss', 'relative_loss', 'gain', 'gain_loss')))
fraction_Y_loss = ggplot(prop_out, aes(x = arrange, y = fraction, fill = category, label = paste0((round(fraction* 100, 2)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  geom_text(size = 3, position = position_stack(vjust = 0.5), fontface = 'bold') +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'loss' = '#0E3F7C',
                               'relative_loss' = '#00AEC8',
                               'gain' = '#D53833',
                               'gain_loss' = '#E3CC98'),
                    name = '') +
  scale_y_continuous(expand = c(0.01,0),
                     breaks = seq(0, 1, 0.2)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 2) +
  labs(x = '', y = 'Fraction of samples\n(n=14,322)')

ggsave(filename = 'Figures_original/Fraction_Y_loss.pdf', plot = fraction_Y_loss, device = 'pdf', width = 4)



##-----------------
## CNA: IGV across the
## whole cohort:
##-----------------
IGV_all = read.csv('Data/04_Loss/IMPACT_IGV_out.seg', sep = '\t')
IGV = IGV_all[which(IGV_all$ID %in% unique(CNA$sample)), ]

















