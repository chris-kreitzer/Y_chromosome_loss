##----------------+
## Chromosome Y incidences 
## among different cancer types;
##----------------+
## Association Studies:
## Starting with
## - Categories
## - CANCER_TYPE
## - CANCER_TYPE_detailed
## - sample site
## - Ancestry
## - Age
##----------------+
## start: 09/01/2021
## revision: 08/25/2022
## revision: 01/08/2023
## revision: 01/10/2023
## revision: 01/16/2023
## revision: 02/01/2023
## 
## chris-kreitzer


clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(cowplot)

##----------------+
## General categories overview;
##----------------+
##
## TODO: 
## 
## CAUTION:
## - ALL SAMPLES INCLUDES (QC = T AND F)
## - Site Sequenced: include local recurrence or not (n~238) and UNKNOWN
##----------------+

cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]

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
prop_out$category = factor(prop_out$category, levels = rev(c('complete_loss',
                                                             'relative_loss',
                                                             'partial_loss',
                                                             'gain',
                                                             'partial_gain',
                                                             'gain_loss',
                                                             'wt')))
fraction_Y_loss = ggplot(prop_out, 
                         aes(x = arrange, 
                             y = fraction, 
                             fill = category, 
                             label = paste0((round(fraction* 100, 1)), '%'))) +
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
  theme_std(base_size = 16, base_line_size = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        aspect.ratio = 2) +
  labs(x = '', y = paste0('Fraction of samples\n(n=', prop_out$n, ')'))

ggsave(filename = 'Figures_original/Fraction_Y_loss.pdf', plot = fraction_Y_loss, device = 'pdf', width = 4)



##----------------+
## Fraction LOY by 
## Cancer type
##----------------+
CancerTypes = data.frame()
for(i in unique(cohort$CANCER_TYPE)){
  total = nrow(cohort)
  data.sub = cohort[which(cohort$CANCER_TYPE == i), ]
  n = length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE == i)])
  for(j in unique(data.sub$classification)){
    fraction = length(data.sub$SAMPLE_ID[which(data.sub$classification == j)]) / n
    out = data.frame(cancer = i,
                     category = j,
                     fraction = fraction,
                     n = n,
                     rel_contribution_cohort = (n/total)*100)
    CancerTypes = rbind(CancerTypes, out)
  }
}


##----------------+
CancerTypes = CancerTypes[!is.na(CancerTypes$category), ]
CancerTypes = CancerTypes[which(CancerTypes$cancer %in% ctypes_keep), ]
CancerTypes$cancer_n = paste0(CancerTypes$cancer, ' (n=', CancerTypes$n, ')')
loss = CancerTypes[which(CancerTypes$category %in% c('complete_loss')), ]
loss = loss[order(loss$fraction), ]
CancerTypes$cancer_n = factor(CancerTypes$cancer_n, levels = rev(loss$cancer_n))
CancerTypes$cancer = factor(CancerTypes$cancer, levels = rev(loss$cancer))
CancerTypes$category = factor(CancerTypes$category, levels = rev(c('complete_loss',
                                                               'relative_loss',
                                                               'partial_loss',
                                                               'gain',
                                                               'partial_gain',
                                                               'gain_loss',
                                                               'wt')))

##----------------+ 
## Visualization;
##----------------+
Fraction_across_cancerTypes = ggplot(CancerTypes, 
                                     aes(x = cancer_n, 
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
  geom_hline(yintercept = seq(0, 1, 0.2), linetype = 'dashed', color = 'white', linewidth = 0.25) +
  theme_std(base_size = 16, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = paste0('Fraction of samples\n(n=', sum(CancerTypes[!duplicated(CancerTypes$cancer), 'n']), ')'))


## combine
plot.out = fraction_Y_loss + theme(legend.position = 'none') + Fraction_across_cancerTypes + labs(y = '')
ggsave(filename = 'Figures_original/LOY_CancerType.pdf', plot = plot.out, device = 'pdf', width = 14)



##-----------------
## Site sequenced;
## Metastatic or Primary samples
##-----------------
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort_ss = cohort[which(cohort$Study_include == 'yes'), ]
cohort_ss$SAMPLE_TYPE[which(cohort_ss$SAMPLE_TYPE == 'Local Recurrence')] = 'Primary'
cohort_ss$SAMPLE_TYPE[which(cohort_ss$SAMPLE_TYPE == 'Unknown')] = 'Metastasis'

prop_out_site = data.frame()
for(i in unique(cohort_ss$CANCER_TYPE)){
  data.sub = cohort_ss[which(cohort_ss$CANCER_TYPE == i), ]
  for(j in unique(data.sub$SAMPLE_TYPE)){
    data.sub2 = data.sub[which(data.sub$SAMPLE_TYPE == j), ]
    n = length(data.sub$SAMPLE_ID[which(data.sub$SAMPLE_TYPE == j)])
    for(k in unique(data.sub2$classification)){
      fraction = length(data.sub2$SAMPLE_ID[which(data.sub2$classification == k)]) / n
      n2 = length(data.sub2$SAMPLE_ID[which(data.sub2$classification == k)])
      out = data.frame(cancer = i,
                       site_sequenced = j,
                       category = k,
                       fraction = fraction,
                       total_cancertype = n,
                       n_ss = n2)
      prop_out_site = rbind(prop_out_site, out)
    }
  }
  rm(data.sub, n, n2, fraction)
}

prop_out_site = prop_out_site[!is.na(prop_out_site$category), ]
prop_out_site = prop_out_site[which(prop_out_site$category == 'complete_loss'), ]
prop_out_site = prop_out_site[which(prop_out_site$cancer %in% ctypes_keep), ]
prop_out_site = rbind(prop_out_site, data.frame(cancer = 'CNS Cancer',
                                                site_sequenced = 'Metastasis',
                                                category = 'complete_loss',
                                                fraction = 0,
                                                total_cancertype = 0,
                                                n_ss = 0))


##-------
## assign binomial confidence interval;
## for plotting;
## Error bars represent binomial CIs
##-------
prop_out_site$lower = ci_lower(n = prop_out_site$n_ss, N = prop_out_site$total_cancertype)
prop_out_site$upper = ci_upper(n = prop_out_site$n_ss, N = prop_out_site$total_cancertype)
prop_out_site$lower[which(prop_out_site$cancer == 'CNS Cancer' & prop_out_site$site_sequenced == "Metastasis")] = 0
prop_out_site$upper[which(prop_out_site$cancer == 'CNS Cancer' & prop_out_site$site_sequenced == "Metastasis")] = 0

prop_out_site$site_sequenced = factor(prop_out_site$site_sequenced, levels = rev(c('Metastasis', 'Primary')))
for(i in unique(prop_out_site$cancer)){
  prop_out_site$name[which(prop_out_site$cancer == i)] = paste0(prop_out_site$cancer[which(prop_out_site$cancer == i)], 
                                                                ' (', 
                                                                prop_out_site$total_cancertype[which(prop_out_site$cancer == i & prop_out_site$site_sequenced == 'Primary')],
                                                                ', ', 
                                                                prop_out_site$total_cancertype[which(prop_out_site$cancer == i & prop_out_site$site_sequenced == 'Metastasis')], ')')
}


##----------------+
## Visualization
##----------------+
prop_out_site$cancer = factor(prop_out_site$cancer, levels = levels(CancerTypes$cancer))
LOY_site_Sequenced = ggplot(prop_out_site, 
                            aes(x = cancer, 
                                y = fraction, 
                                ymin = lower, 
                                ymax = upper)) +
  geom_pointrange(aes(color = site_sequenced), position = position_dodge2(width = .5), size = as_points(1), fatten = 4) +
  stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = ..x.. - 0.2, yend = ..y..)) +
  stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = ..x.. + 0.2, yend = ..y..)) +
  scale_color_manual(values = c('Primary' = '#7c93b5',
                                'Metastasis' = '#35538c'),
                     name = 'Site sequenced') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) + #labels = function(x) 100*x) +
  theme_std(base_size = 16, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'top',
        legend.justification = c(1, 1)) +
  #aspect.ratio = 0.2) +
  labs(x = NULL, y = '', color = NULL)


##-------
## All tumors combined
##-------
Site_Sequenced_all = ggplot(prop_out_site, aes(x = site_sequenced, y = fraction, 
                                               color = site_sequenced)) +
  geom_quasirandom(width = 0.25, alpha = 0.95) +
  stat_summary(fun.y = "mean", geom = "crossbar", size = 0.5, color = 'black', width = 0.35) +
  stat_summary(fun.y = "mean", geom = "point", size = 3, color = 'red') +
  stat_compare_means(label.x = 1.2,
                     label.y = 0.9) +
  scale_color_manual(values = c('Primary' = '#7c93b5',
                                'Metastasis' = '#35538c'),
                     name = 'Site sequenced') + 
  scale_y_continuous(expand = c(0.01,0),
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1)) +
  theme_std(base_size = 16, base_line_size = 1) +
  theme(legend.position = 'none') +
  labs(x = '', y = paste0('Fraction of samples\n(n=', prop_out$n, ')'))



##----------------+
## combine the two plots
## LOY rates and Site_sequenced
##----------------+

aa = Site_Sequenced_all + LOY_site_Sequenced + theme(axis.text.x = element_blank()) + plot_layout(widths = c(1, 5))
bb = fraction_Y_loss + Fraction_across_cancerTypes + theme(legend.position = 'bottom', axis.title.y = element_blank())

aa
bb

dd = aa / bb + plot_layout(ncol = 2, nrow = 2)
ggsave_golden(filename = 'CohortOverview.pdf', plot = dd, width = 26)


AA = Site_Sequenced_all
BB = LOY_site_Sequenced + theme(axis.text.x = element_blank())
CC = fraction_Y_loss + theme(legend.position = 'none')
DD = Fraction_across_cancerTypes + labs(y = '') + theme(legend.position = 'bottom')
XX = (AA+BB) / (CC+DD)
ggsave_golden(filename = 'Figures_original/Cohort_Incidences_combined.pdf', plot = XX, width = 21)



##----------------+
## CancerType detailed
##----------------+
CancerTypesDetailed = data.frame()
for(i in unique(cohort$CANCER_TYPE)){
  data.sub = cohort[which(cohort$CANCER_TYPE == i), ]
  n = length(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE == i)])
  for(j in unique(data.sub$CANCER_TYPE_DETAILED_ritika)){
    data.sub2 = data.sub[which(data.sub$CANCER_TYPE_DETAILED_ritika == j), ]
    n2 = length(data.sub2$SAMPLE_ID[which(data.sub2$CANCER_TYPE_DETAILED_ritika == j)])
    cancertype = j
    for(k in unique(data.sub2$classification)){
      fraction = length(data.sub2$SAMPLE_ID[which(data.sub2$classification == k)]) / n2
      out = data.frame(cancer = i,
                       cancerdetailed = cancertype,
                       category = k,
                       fraction = fraction,
                       n = n2)
      CancerTypesDetailed = rbind(CancerTypesDetailed, out)
    }
  }
  rm(data.sub, data.sub2, n, fraction)
}


##-------
## Further summary
##-------
CancerTypesDetailed_loss = CancerTypesDetailed[which(CancerTypesDetailed$category == 'complete_loss'), ]
CancerTypesDetailed_loss$cancerdetailed2 = NA
for(i in 1:nrow(CancerTypesDetailed_loss)){
  CancerTypesDetailed_loss$cancerdetailed2[i] = gsub(pattern = CancerTypesDetailed_loss$cancer[i], 
                                                     replacement = '', 
                                                     x = CancerTypesDetailed_loss$cancerdetailed[i])
}

CancerTypesDetailed_loss$cancerdetailed2 = substr(x = CancerTypesDetailed_loss$cancerdetailed2, start = 1, stop = nchar(CancerTypesDetailed_loss$cancerdetailed2)-1)
CancerTypesDetailed_loss$name = paste0(CancerTypesDetailed_loss$cancerdetailed2, ', ', CancerTypesDetailed_loss$n, ')')

##-------
## at least:
## - 2 sub-types
## - 10 cases per sub-type 
##-------
CancerTypesDetailed_loss = CancerTypesDetailed_loss[which(CancerTypesDetailed_loss$n >= 10), ]
names_exclude = sort(table(CancerTypesDetailed_loss$cancer), decreasing = T)
names_exclude = names(names_exclude[names_exclude > 1])
CancerTypesDetailed_loss = CancerTypesDetailed_loss[which(CancerTypesDetailed_loss$cancer %in% names_exclude), ]
CancerTypesDetailed_loss$fraction = CancerTypesDetailed_loss$fraction * 100

plot_list = list()
for(i in unique(CancerTypesDetailed_loss$cancer)){
  plot_i = ggplot(CancerTypesDetailed_loss[which(CancerTypesDetailed_loss$cancer == i), ], 
                aes(x = reorder(name, fraction), y = fraction)) +
    geom_bar(stat = 'identity', color = 'black', fill = 'black') +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, round(max(CancerTypesDetailed_loss$fraction[which(CancerTypesDetailed_loss$cancer == i)]), 1) + 1),
                       sec.axis = dup_axis()) +
    geom_hline(yintercept = seq(0, 100, 20), color = 'white', linetype = 'dashed', linewidth = 0.25) +
    theme_std(base_size = 16) +
    theme(axis.line.x.bottom = element_blank(),
          axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          axis.title.x.bottom = element_blank(),
          plot.margin = unit(c(0,0,0,0), units = 'mm')) +
    labs(x = '', y = 'Percent of cases with LOY', title = i)
  
  plot_list[[i]] = plot_i
}

for(i in 1:length(plot_list)){
  ggsave_golden(filename = paste0('Figures_original/CTD_', names(plot_list[i]), '.pdf'), plot = plot_list[[i]], width = 8)
}


##-------
## save individual plot;
##-------
ggsave_golden(filename = 'Figures_original/CTD_Esophagastric.pdf', plot = plot_list$`Esophagogastric Cancer`, width = 8)
plot_list$`Esophagogastric Cancer` / plot_list$`Renal Cell Carcinoma` / plot_list$`Pancreatic Cancer` / plot_list$`Ampullary Cancer`





##----------------+
## TCGA to MSK cohort comparison
##----------------+
TCGA_MSK = read.csv('Data/04_Loss/TCGA_MSK_LOYrates.txt', sep = '\t')
TCGA_MSK$Fraction_LOY = gsub(pattern = ',', replacement = '\\.', x = TCGA_MSK$Fraction_LOY)
TCGA_MSK$Fraction_LOY.1 = gsub(pattern = ',', replacement = '\\.', x = TCGA_MSK$Fraction_LOY.1)
TCGA_MSK$Fraction_LOY = as.numeric(as.character(TCGA_MSK$Fraction_LOY))
TCGA_MSK$Fraction_LOY.1 = as.numeric(as.character(TCGA_MSK$Fraction_LOY.1))
TCGA_MSK$name = paste0(TCGA_MSK$TCGA_Abb, ' - ', TCGA_MSK$TumorType_MSKCC, ' (n=', TCGA_MSK$n.1, ')')
comparison_plot = ggplot(TCGA_MSK,
                         aes(x = Fraction_LOY.1,
                             y = Fraction_LOY,
                             color = name, 
                             size = n,
                             label = paste0(TCGA_Abb, ' (n=', n, ')'))) +
  geom_point() +
  geom_text_repel() +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed',
              linewidth = 0.35,
              color = 'grey35') +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 1)) +
  scale_size_continuous(breaks = c(50, 100, 300, 500), range = c(1,8)) +
  theme_std(base_size = 16, base_line_size = 0.5) +
  annotate('text',
           x = 0.10,
           y = 0.9,
           label = paste0('r = ', round(cor.test(TCGA_MSK$Fraction_LOY, TCGA_MSK$Fraction_LOY.1)[[4]][[1]], 2), 
                          '\np = ', round(cor.test(TCGA_MSK$Fraction_LOY, TCGA_MSK$Fraction_LOY.1)$p.value, 10)),
           hjust = 0, 
           vjust = 0,
           family = 'ArialMT',
           size = 6) +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, linewidth = 2),
        panel.background = element_rect(fill = 'white'),
        legend.position = 'top') +
  guides(color = 'none') +
  labs(x = 'Fraction LOY MSKCC',
       y = 'Fraction LOY TCGA',
       color = 'TCGA Abbrevation - MSKCC Cancer type',
       size = 'Sample size: TCGA study')

comparison_plot

ggsave_golden(filename = 'Figures_original/TCGA_MSK_LOY_comparison.pdf', plot = comparison_plot, width = 12)


##----------------+
## Race category
##----------------+

# African (AFR), 
# European (EUR), 
# East Asian (EAS), 
# Native American (NAM) 
# and South Asian (SAS), 
# Ashkenazi Jewish (ASJ) 
Ancestry = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/admixture_results.50k.2021-09-14.txt', sep = '\t')
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
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
                     fraction = fraction,
                     n = n)
    Ancestry_out = rbind(Ancestry_out, out)
  }
}

Ancestry_out = Ancestry_out[!is.na(Ancestry_out$category), ]
Ancestry_out = Ancestry_out[!Ancestry_out$ancestry %in% 'ADMIX_OTHER', ]
Ancestry_out$ancestry = paste0(Ancestry_out$ancestry, ' (n=', Ancestry_out$n, ')')
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
  theme_std(base_size = 16, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = paste0('Fraction of samples\n(n=', sum(unique(Ancestry_out$n)),')'))


ggsave(filename = 'Figures_original/Fraction_LOY_Ancestry.pdf', plot = Ancestry_LOY, device = 'pdf', width = 8)


## lm model if there is an difference among the groups



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
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
cohort$Y_call = ifelse(cohort$classification %in% c('partial_loss', 'relative_loss', 'complete_loss'), 'LOY', 'nonLOY')
impact = cohort
impact$Age_Sequencing = as.integer(impact$Age_Sequencing)

Loss_age = data.frame()
for(i in seq(15, 90, 5)){
  da = impact[between(impact$Age_Sequencing, i-4, i), ]
  da = da[!is.na(da$Age_Sequencing), ]
  
  if(length(table(da$Y_call)) == 2){
    ratio = (table(da$Y_call)[[1]] / sum(table(da$Y_call))) * 100
    n = sum(table(da$Y_call))
    out = data.frame(group = i,
                     LOY = ratio,
                     n = n)
    
  } else if (length(table(da$Y_call)) == 1) {
    table.ratio = table(da$Y_call)
    ratio = ifelse(names(table.ratio) == 'nonLOY', 0, 100)
    n = table.ratio[[1]]
    out = data.frame(group = i,
                     LOY = ratio,
                     n = n)
  } else {
    ratio = NA
    n = dim(da)[[1]]
    out = data.frame(group = i,
                     LOY = ratio, 
                     n = n)
  }
  Loss_age = rbind(Loss_age, out)
}

Loss_age$group_plot = paste0(Loss_age$group - 4, '-', Loss_age$group)
Loss_age$group = factor(Loss_age$group, levels = Loss_age$group)
Loss_age = Loss_age[-1, ]
Loss_age$group_plot[1] = '18-20'

##----------------+ 
## Visualization;
##----------------+
Age_distribution = ggplot(Loss_age, aes(x = group, y = LOY)) +
  geom_bar(stat = 'identity', color = 'white', fill = 'grey35') +
  scale_x_discrete(labels = Loss_age$group_plot, expand = c(0.05, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 40)) +
  geom_hline(yintercept = seq(10, 30, 10), linetype = 'solid', color = 'white') +
  theme_std(base_size = 16, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = 'Age-group', y = 'LOY [%]')

Age_distribution

ggsave_golden(filename = 'Figures_original/LOY_AgeDistribution.pdf', plot = Age_distribution, width = 6)




#' out