##----------------+
## Chromosome Y incidences 
## among different cancer types;
##----------------+
## Association Studies:
## Starting with
## - categories
## - cancer types ritika (oncotree-code)
## - cancer types detailed
## - sample site
## - Ancestry
## - Age
##----------------+
## start: 09/01/2021
## revision: 08/25/2022
## revision: 01/08/2023
## revision: 01/10/2023
## chris-kreitzer


clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
source('Scripts/UtilityFunctions.R')
source('Scripts/plot_theme.R')
library(RColorBrewer)
library(patchwork)



##----------------+
## General categories overview;
##----------------+
##
## CAUTION:
## 
## ALL SAMPLES INCLUDES (QC = T AND F)
##----------------+

cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')

n = length(unique(cohort$SAMPLE_ID))
prop_out = data.frame()
for(i in unique(cohort$classification)){
  fraction = length(cohort$SAMPLE_ID[which(cohort$classification == i)]) / n
  out = data.frame(category = i,
                   fraction = fraction)
  prop_out = rbind(prop_out, out)
}

prop_out$arrange = 1
prop_out$n = n
prop_out$category = factor(prop_out$category, levels = rev(c('wt', 'complete_loss', 'relative_loss', 'partial_loss', 'gain', 'partial_gain', 'gain_loss')))
fraction_Y_loss = ggplot(prop_out, aes(x = arrange, y = fraction, fill = category, label = paste0((round(fraction* 100, 2)), '%'))) +
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
  scale_y_continuous(expand = c(0.01,0),
                     breaks = seq(0, 1, 0.2)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 2) +
  labs(x = '', y = paste0('Fraction of samples\n(n=', prop_out$n, ')'))

ggsave(filename = 'Figures_original/Fraction_Y_loss.pdf', plot = fraction_Y_loss, device = 'pdf', width = 4)


##----------------+
## Cancer Types
##----------------+
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
CancerTypes = data.frame()
for(i in unique(cohort$CANCER_TYPE_ritika)){
  data.sub = cohort[which(cohort$CANCER_TYPE_ritika == i), ]
  n = length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE_ritika == i)])
  for(j in unique(data.sub$classification)){
    fraction = length(data.sub$SAMPLE_ID[which(data.sub$classification == j)]) / n
    out = data.frame(cancer = i,
                     category = j,
                     fraction = fraction,
                     n = n)
    CancerTypes = rbind(CancerTypes, out)
  }
}

##----------------+
CancerTypes = CancerTypes[!is.na(CancerTypes$category), ]
CancerTypes = CancerTypes[!CancerTypes$cancer %in% c('Vulva/Vagina (VULVA)', 'Ovary/Fallopian Tube (OVARY)'), ]
CancerTypes$cancer = gsub("\\s*\\([^\\)]+\\)", "", CancerTypes$cancer)
CancerTypes$cancer = paste0(CancerTypes$cancer, ' (n=', CancerTypes$n, ')')
loss = CancerTypes[which(CancerTypes$category %in% c('complete_loss')), ]
loss = loss[order(loss$fraction), ]
CancerTypes$cancer = factor(CancerTypes$cancer, levels = rev(loss$cancer))
CancerTypes$category = factor(CancerTypes$category, levels = c('complete_loss',
                                                               'relative_loss',
                                                               'partial_loss',
                                                               'gain',
                                                               'partial_gain',
                                                               'gain_loss',
                                                               'wt'))

##----------------+ 
## Visualization;
##----------------+
Fraction_across_cancerTypes = ggplot(CancerTypes, 
                                     aes(x = cancer, 
                                         y = fraction, 
                                         fill = category, 
                                         label = paste0((round(fraction* 100, 2)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                              'complete_loss' = '#0E3F7C',
                              'partial_loss' = '#00AEC8',
                              'relative_loss' = '#474089',
                              'gain' = '#D53833',
                              'partial_gain' = '#E3CC98',
                              'gain_loss' = '#f9f8d5'),
                    name = '') +
  scale_y_continuous(expand = c(0.01,0),
                     breaks = seq(0, 1, 0.2)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = paste0('Fraction of samples\n(n=', sum(CancerTypes[!duplicated(CancerTypes$cancer), 'n']), ')'))


## combine
plot.out = fraction_Y_loss + theme(legend.position = 'none') + Fraction_across_cancerTypes + labs(y = '')
ggsave(filename = 'Figures_original/LOY_CancerType_Ritika.pdf', plot = plot.out, device = 'pdf', width = 14)



##----------------+
## CancerType detailed
##----------------+
CancerTypesDetailed = data.frame()
for(i in unique(cohort$CANCER_TYPE_DETAILED_ritika)){
  data.sub = cohort[which(cohort$CANCER_TYPE_DETAILED_ritika == i), ]
  n = length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE_DETAILED_ritika == i)])
  for(j in unique(data.sub$classification)){
    fraction = length(data.sub$SAMPLE_ID[which(data.sub$classification == j)]) / n
    out = data.frame(cancer = i,
                     category = j,
                     fraction = fraction,
                     n = n,
                     n_all = n/length(unique(cohort$SAMPLE_ID)) * 100)
    CancerTypesDetailed = rbind(CancerTypesDetailed, out)
  }
}


##----------------+
CancerTypesDetailed = CancerTypesDetailed[which(CancerTypesDetailed$n_all >= 0.5), ]
CancerTypesDetailed = CancerTypesDetailed[!is.na(CancerTypesDetailed$category), ]
CancerTypesDetailed$cancer[which(CancerTypesDetailed$cancer == 'Undifferentiated Pleomorphic Sarcoma/Malignant Fibrous Histiocytoma/High-Grade Spindle Cell Sarcoma (MFH)')] = 'High-Grade Spindle Cell Sarcoma (MFH)'
#CancerTypes$cancer = gsub("\\s*\\([^\\)]+\\)", "", CancerTypes$cancer)
#CancerTypesDetailed$cancer = paste(CancerTypesDetailed$cancer, paste0('(n=', CancerTypes$n, ')'), sep = '\n')
loss = CancerTypesDetailed[which(CancerTypesDetailed$category %in% c('complete_loss')), ]
loss = loss[order(loss$fraction), ]
CancerTypesDetailed$cancer = factor(CancerTypesDetailed$cancer, levels = rev(loss$cancer))
CancerTypesDetailed$category = factor(CancerTypesDetailed$category, levels = rev(c('complete_loss',
                                                               'relative_loss',
                                                               'partial_loss',
                                                               'gain',
                                                               'partial_gain',
                                                               'gain_loss',
                                                               'wt')))


##----------------+
## Visualization
##----------------+
Fraction_CancerTypesDetailed = ggplot(CancerTypesDetailed, 
                                     aes(x = cancer, 
                                         y = fraction, 
                                         fill = category, 
                                         label = paste0((round(fraction* 100, 2)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  coord_flip() +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'complete_loss' = '#0E3F7C',
                               'partial_loss' = '#00AEC8',
                               'relative_loss' = '#474089',
                               'gain' = '#D53833',
                               'partial_gain' = '#E3CC98',
                               'gain_loss' = '#f9f8d5'),
                    name = '') +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 1, 0.2),
                     sec.axis = dup_axis()) +
  theme_std(base_size = 14, base_line_size = 1) +
  labs(x = '', y = paste0('Fraction of samples\n(n=', sum(CancerTypes[!duplicated(CancerTypes$cancer), 'n']), ')'))

Fraction_CancerTypesDetailed

## combine
plot.out = fraction_Y_loss + theme(legend.position = 'none') + Fraction_CancerTypesDetailed + labs(y = '')
ggsave(filename = 'Figures_original/LOY_CancerType_Detailed_Ritika.pdf', plot = plot.out, device = 'pdf', width = 14)



##-----------------
## Sample site;
## Metastatic or Primary samples
##-----------------
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')

prop_out_site = data.frame()
for(i in unique(cohort$SAMPLE_TYPE)){
  data.sub = cohort[which(cohort$SAMPLE_TYPE == i), ]
  n = length(cohort$SAMPLE_ID[which(cohort$SAMPLE_TYPE == i)])
  for(j in unique(data.sub$classification)){
    fraction = length(data.sub$SAMPLE_ID[which(data.sub$classification == j)]) / n
    out = data.frame(cancer = i,
                     category = j,
                     fraction = fraction,
                     n = n,
                     n_all = n/length(unique(cohort$SAMPLE_ID)) * 100)
    prop_out_site = rbind(prop_out_site, out)
  }
}

prop_out_site = prop_out_site[!is.na(prop_out_site$category), ]
prop_out_site$cancer = factor(prop_out_site$cancer, levels = c('Primary', 'Metastasis', 'Local Recurrence', 'Unknown'))
prop_out_site$category = factor(prop_out_site$category, levels = rev(c('complete_loss',
                                                                       'relative_loss',
                                                                       'partial_loss',
                                                                       'gain',
                                                                       'partial_gain',
                                                                       'gain_loss',
                                                                       'wt')))

##----------------+
## Visualization
##----------------+
Fraction_CancerSite = ggplot(prop_out_site,
                             aes(x = cancer, 
                                 y = fraction, 
                                 fill = category, 
                                 label = paste0((round(fraction* 100, 2)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
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
  labs(x = '', y = paste0('Fraction of samples\n(n=', sum(prop_out_site[!duplicated(prop_out_site$cancer), 'n']), ')'))

ggsave(filename = 'Figures_original/LOY_CancerSite.pdf', plot = Fraction_CancerSite, device = 'pdf', width = 14)



##-----------------
## Race category
##-----------------

# African (AFR), 
# European (EUR), 
# East Asian (EAS), 
# Native American (NAM) 
# and South Asian (SAS), 
# Ashkenazi Jewish (ASJ) 

Ancestry = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/admixture_results.50k.2021-09-14.txt', sep = '\t')
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
Ancestry_short = merge(cohort[,c('Normal_file', 'classification')], Ancestry[,c('Sample', 'ancestry_label')],
                       by.x = 'Normal_file', by.y = 'Sample', all.x = T)
Ancestry_short = Ancestry_short[!is.na(Ancestry_short$ancestry_label), ]


Ancestry_out = data.frame()
for(i in unique(Ancestry_short$ancestry_label)){
  data.sub = Ancestry_short[which(Ancestry_short$ancestry_label == i), ]
  n = length(Ancestry_short$Normal_file[which(Ancestry_short$ancestry_label == i)])
  for(j in unique(data.sub$classification)){
    fraction = length(data.sub$Normal_file[which(data.sub$classification == j)]) / n
    out = data.frame(ancestry = i,
                     category = j,
                     fraction = fraction)
    Ancestry_out = rbind(Ancestry_out, out)
  }
}

Ancestry_out = Ancestry_out[!is.na(Ancestry_out$category), ]
loss = Ancestry_out[which(Ancestry_out$category %in% c('complete_loss')), ]
loss = loss[order(loss$fraction), ]
Ancestry_out$ancestry = factor(Ancestry_out$ancestry, levels = rev(loss$ancestry))
Ancestry_out$category = factor(Ancestry_out$category, levels = rev(c('complete_loss',
                                                                                   'relative_loss',
                                                                                   'partial_loss',
                                                                                   'gain',
                                                                                   'partial_gain',
                                                                                   'gain_loss',
                                                                                   'wt')))

##-----------------
## Ancestry Visualization:
##-----------------
Ancestry_LOY = ggplot(Ancestry_out, 
                      aes(x = ancestry, 
                          y = fraction, 
                          fill = category, 
                          label = paste0((round(fraction* 100, 2)), '%'))) +
  geom_bar(stat = 'identity', position = 'fill') +
  geom_hline(yintercept = seq(0, 1, 0.2), color = 'white', linewidth = 0.25, linetype = 'dashed') +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'complete_loss' = '#0E3F7C',
                               'partial_loss' = '#00AEC8',
                               'relative_loss' = '#474089',
                               'gain' = '#D53833',
                               'partial_gain' = '#E3CC98',
                               'gain_loss' = '#f9f8d5'),
                    name = '') +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 1, 0.2)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = 'Fraction of samples\n(n=14,322)')

Ancestry_LOY

ggsave(filename = 'Figures_original/Fraction_LOY_Ancestry.pdf', plot = AncestryYloss, device = 'pdf', width = 8)


##----------------+
## proportion test;
## any association with ancestry
##----------------+
# fractions = c(0.33, 0.32, 0.33)
# res = chisq.test(fractions, p = c(1/3, 1/3, 1/3))
# test = read.csv('Data/04_Loss/IMPACT_copynumber_out.txt', sep = '\t')


##-----------------
## Association with age;
##-----------------
cohort = readRDS('Data/signedOut/Cohort_07132022.rds')

impact = merge(cohort$IMPACT_cohort, cohort$IMPACT_clinicalAnnotation[, c('SAMPLE_ID', 'AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)')],
               by = 'SAMPLE_ID', all.x = T)

colnames(impact)[ncol(impact)] = 'Age_Sequencing'
impact = impact[which(impact$Y_CNA == TRUE & !is.na(impact$Y_call)), ]


Loss_age = data.frame()
for(i in seq(15, 90, 5)){
  da = impact[between(impact$Age_Sequencing, i-4, i), ]
  da = da[!is.na(da$Age_Sequencing), ]
  
  if(length(table(da$Y_call)) == 2){
    ratio = (table(da$Y_call)[2] / sum(table(da$Y_call))) * 100
    n = sum(table(da$Y_call))
    out = data.frame(group = i,
                     ratio = ratio,
                     n = n)
    
  } else if (length(table(da$Y_call)) == 1) {
    table.ratio = table(da$Y_call)
    ratio = ifelse(names(table.ratio) == 'intact_Y_chrom', 0, 100)
    n = table.ratio[[1]]
    out = data.frame(group = i,
                     ratio = ratio,
                     n = n)
  } else {
    ratio = NA
    n = dim(da)[[1]]
    out = data.frame(group = i,
                     ratio = ratio, 
                     n = n)
  }
  Loss_age = rbind(Loss_age, out)
}

Loss_age$group_plot = paste0(Loss_age$group - 4, '-', Loss_age$group)
Loss_age$group = factor(Loss_age$group, levels = Loss_age$group)

## Visualization

Age_distribution = ggplot(Loss_age, aes(x = group, y = ratio)) +
  geom_bar(stat = 'identity', color = 'white', fill = 'grey35') +
  scale_x_discrete(labels = Loss_age$group_plot, expand = c(0.05, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 40)) +
  geom_hline(yintercept = seq(10, 30, 10), linetype = 'solid', color = 'white') +
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, color = 'black'),
        panel.background = element_rect(fill = NA, size = 1.9),
        line = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.y = element_text(angle = 90)) +
  labs(x = 'Age Group', y = 'Chromosome Y Loss incidence [%]')

Age_distribution
ggsave_golden(filename = 'Figures_original/Y_loss_AgeDistribution.pdf', plot = Age_distribution, width = 6)







IMPACT_incidence = data.frame()
for(i in unique(IMPACT_data$CANCER_TYPE_ritika)){
  n.all = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE %in% c('Primary', 'Metastasis') & IMPACT_data$CANCER_TYPE_ritika == i)])
  n.primary = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Primary' & IMPACT_data$CANCER_TYPE_ritika == i)])
  n.metastasis = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Metastasis' & IMPACT_data$CANCER_TYPE_ritika == i)])
  
  primary.frequency = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$SAMPLE_TYPE == 'Primary' & IMPACT_data$CANCER_TYPE == i)]) / n.primary
  metastatic.frequency = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$SAMPLE_TYPE == 'Metastasis' & IMPACT_data$CANCER_TYPE == i)]) / n.metastasis
  cancer_median = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$CANCER_TYPE == i)]) / n.all
  cancer_summary = data.frame(CancerType = i,
                              Type = c('Primary', 'Metastasis'),
                              value = c(primary.frequency, metastatic.frequency),
                              median_cancerType = cancer_median)
  
  IMPACT_incidence = rbind(IMPACT_incidence, cancer_summary)
  rm(n.all, cancer_summary)
  rm(primary.frequency)
  rm(metastatic.frequency)
  rm(cancer_median)
}

#' add sample amount
IMPACT_incidence$mergedCancerType = NA
for(i in unique(names(IMPACT_top_20))){
  IMPACT_incidence$mergedCancerType[which(IMPACT_incidence$CancerType == i)] = IMPACT_top_20[which(names(IMPACT_top_20) == i)][[1]]
}

IMPACT_incidence$mergedCancerType = paste0(IMPACT_incidence$CancerType, ' (n=', IMPACT_incidence$mergedCancerType, ')')
write.table(IMPACT_incidence, file = 'Data/04_Loss/IMPACT_Y_loss_incidences_top20.txt', sep = '\t', row.names = F, quote = F)


##-----------------
## Visualization
##-----------------
color_selected = brewer.pal(n = 6, name = 'Paired')
color_selected = color_selected[c(5,6,1,2)]

IMPACT_incidence_plot = ggplot(IMPACT_incidence, 
                               aes(x = reorder(mergedCancerType, median_cancerType), 
                                   y = (value * 100), 
                                   fill = Type)) + 
  
  geom_bar(stat = 'identity', position = position_dodge(0.82), width = 0.8) +
  scale_fill_manual(values = c('Primary' = color_selected[4],
                               'Metastasis' = color_selected[2]),
                    guide = guide_legend(direction = 'horizontal',
                                         title = 'Site',
                                         label.theme = element_text(size = 14))) +
  scale_y_continuous(expand = c(0.005, 0)) +
  geom_hline(yintercept = seq(0, 60, 20), color = 'grey35', linetype = 'dashed', size = 0.2) +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.line.x = element_line(size = 0.4),
        axis.text = element_text(size = 10, color = 'black')) +
  coord_flip() +
  labs(y = '% Y chromosome loss', x = '')

IMPACT_incidence_plot  

ggsave_golden(IMPACT_incidence_plot, filename = 'Figures_original/IMPACT_Y_Loss_top20.pdf', width = 9)


##-----------------
## Y-loss across SampleType
##-----------------
IMPACT_incidence$Type = factor(IMPACT_incidence$Type, levels = c('Primary', 'Metastasis'))
Incidence_site = ggplot(IMPACT_incidence, aes(x = Type, y = (value*100), color = Type)) +
  geom_boxplot(width = 0.4, size = 0.85) +
  geom_jitter(width = 0.15) +
  scale_color_manual(values = c('Primary' = color_selected[4],
                                'Metastasis' = color_selected[2])) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 80)) +
  theme(legend.position = 'none',
        #panel.background = element_rect(fill = NA),
        axis.text = element_text(size = 12, color = 'black'),
        aspect.ratio = 2.2,
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  labs(x = '', y = 'Y-chromosome loss [%]')

Incidence_site

ggsave_golden(filename = 'Figures_original/Y_loss_SITE_top20.pdf', plot = Incidence_site, width = 6)






