##-----------------
## Association Studies
## Starting with a general overview
## dig deeper continously;
## 
## start: 08/20/2022
## revision: 10/28/2022
## chris-kreitzer
##----------------


clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/')
source('Scripts/plot_theme.R')
library(patchwork)


##-----------------
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
## Fraction of Y-chromosome
## loss across cancer types
##-----------------
data = readRDS('Data/signedOut/Cohort_07132022.rds')
CNA = data$IMPACT_Y_classification_final

Race = read.csv('Data/signedOut/RaceCategory_cBIO.tsv', sep = '\t')
Race = Race[,c('Sample.ID', 'Ethnicity.Category', 'Race.Category')]
Clinical = data$IMPACT_clinicalAnnotation

Clinical = merge(Clinical, Race, by.x = 'SAMPLE_ID', by.y = 'Sample.ID', all.x = T)

Y_clinical = merge(CNA, Clinical[, c('SAMPLE_ID', 'CANCER_TYPE', 'Race.Category')],
                   by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)

##-----------------
all_out = data.frame()
for(i in unique(Y_clinical$CANCER_TYPE)){
  data.sub = Y_clinical[which(Y_clinical$CANCER_TYPE == i), ]
  n = length(Y_clinical$sample[which(Y_clinical$CANCER_TYPE == i)])
  for(j in unique(data.sub$classification)){
    fraction = length(data.sub$sample[which(data.sub$classification == j)]) / n
    out = data.frame(cancer = i,
                     category = j,
                     fraction = fraction)
    all_out = rbind(all_out, out)
  }
}

##-----------------
all_out = all_out[which(all_out$cancer %in% names(IMPACT_top_20)), ]
loss = all_out[which(all_out$category == 'loss'), ]
loss = loss[order(loss$fraction), ]
all_out$cancer = factor(all_out$cancer, levels = rev(loss$cancer))

#' Vis:
Fraction_across_cancerTypes = ggplot(all_out, aes(x = cancer, y = fraction, fill = category, label = paste0((round(fraction* 100, 2)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  #geom_text(size = 3, position = position_stack(vjust = 0.5), fontface = 'bold') +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'loss' = '#0E3F7C',
                               'relative_loss' = '#00AEC8',
                               'gain' = '#D53833',
                               'gain_loss' = '#E3CC98'),
                    name = '') +
  scale_y_continuous(expand = c(0.01,0),
                     breaks = seq(0, 1, 0.2)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = 'Fraction of samples\n(n=14,322)')
  

## combine
plot.out = fraction_Y_loss + theme(legend.position = 'none') + Fraction_across_cancerTypes + labs(y = '')
ggsave(filename = 'Figures_original/Fraction_Y_loss_acrossCancer.pdf', plot = plot.out, device = 'pdf', width = 14)



##-----------------
## Race category
##-----------------
Ancestry = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/admixture_results.50k.2021-09-14.txt', sep = '\t')
CNA = data$IMPACT_Y_classification_final
CNA$PatientID = substr(CNA$sample, start = 1, stop = 9)
Y_ancestry = merge(CNA, Ancestry[,c('Patient', 'ancestry_label')], by.x = 'PatientID', by.y = 'Patient', all.x = T) 

# African (AFR), 
# European (EUR), 
# East Asian (EAS), 
# Native American (NAM) 
# and South Asian (SAS), 
# Ashkenazi Jewish (ASJ) 

all_out = data.frame()
for(i in unique(Y_ancestry$ancestry_label)){
  data.sub = Y_ancestry[which(Y_ancestry$ancestry_label == i), ]
  n = length(Y_ancestry$sample[which(Y_ancestry$ancestry_label == i)])
  for(j in unique(data.sub$classification)){
    fraction = length(data.sub$sample[which(data.sub$classification == j)]) / n
    out = data.frame(ancestry = i,
                     category = j,
                     fraction = fraction)
    all_out = rbind(all_out, out)
  }
}

##-----------------
## Ancestry Visualization:
##-----------------
AncestryYloss = ggplot(all_out, aes(x = ancestry, 
                                    y = fraction, 
                                    fill = category, 
                                    label = paste0((round(fraction* 100, 2)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  #geom_text(size = 3, position = position_stack(vjust = 0.5), fontface = 'bold') +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'loss' = '#0E3F7C',
                               'relative_loss' = '#00AEC8',
                               'gain' = '#D53833',
                               'gain_loss' = '#E3CC98'),
                    name = '') +
  scale_y_continuous(expand = c(0.01,0),
                     breaks = seq(0, 1, 0.2)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = 'Fraction of samples\n(n=14,322)')

ggsave(filename = 'Figures_original/Fraction_Y_loss_Race.pdf', plot = AncestryYloss, device = 'pdf', width = 8)


## proportion test
# fractions = c(0.33, 0.32, 0.33)
# res = chisq.test(fractions, p = c(1/3, 1/3, 1/3))
# 
# test = read.csv('Data/04_Loss/IMPACT_copynumber_out.txt', sep = '\t')

