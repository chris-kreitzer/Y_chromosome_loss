##----------------+
## Look into Y-chromosome incidences among different cancer types
## Make an overview graph; as Bastien initially did (back in 2019)
## 
## start: 09/01/2021
## revision: 08/25/2022
## 
## chris-kreitzer

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
source('Scripts/UtilityFunctions.R')

## Libraries and Input
library(RColorBrewer)
cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
samples_loss = unique(as.character(cohort$IMPACT_cohort$SAMPLE_ID[which(cohort$IMPACT_cohort$LOY == 'no')]))

##-----------------
## Sample consideration
##-----------------
Binary = read.csv('Data/04_Loss/IMPACT.binaryLoss_out.txt', sep = '\t')
Binary = Binary[which(Binary$sample.id %in% samples_loss), ]

#' Just use the FACETS QC TRUE samples
Facets_QC = read.csv('Data/04_Loss/QC_metrics.txt', sep = '\t')
samples_pass = Facets_QC[which(Facets_QC$QC == 'TRUE'), ]
Y_CNA = Binary[which(Binary$sample.id %in% samples_pass$ID), ]

IMPACT = cohort$IMPACT_cohort
IMPACT$Y_CNA = NA

for(i in 1:nrow(IMPACT)){
  if(IMPACT$SAMPLE_ID[i] %in% Y_CNA$sample.id){
    IMPACT$Y_CNA[i] = TRUE
  } else {
    IMPACT$Y_CNA[i] = FALSE
  }
}

cohort = list(IMPACT_cohort = IMPACT,
              IMPACT_clinicalAnnotation = cohort$IMPACT_clinicalAnnotation,
              IMPACT_LOY = cohort$IMPACT_LOY,
              IMPACT_binaryY_call = Y_CNA)

saveRDS(cohort, file = '~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')


##-----------------
## Incidences
##-----------------
cohortData = readRDS('Data/signedOut/Cohort_07132022.rds')
IMPACT_data = cohortData$IMPACT_cohort

#' top 20 cancer types
IMPACT_data = IMPACT_data[which(IMPACT_data$Y_CNA == TRUE), ]
IMPACT_top_20 = table(IMPACT_data$CANCER_TYPE)
IMPACT_top_20 = IMPACT_top_20[order(IMPACT_top_20, decreasing = T)]
IMPACT_top_20 = IMPACT_top_20[1:20]

IMPACT_incidence = data.frame()
for(i in unique(names(IMPACT_top_20))){
  n.all = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE %in% c('Primary', 'Metastasis') & IMPACT_data$CANCER_TYPE == i)])
  n.primary = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Primary' & IMPACT_data$CANCER_TYPE == i)])
  n.metastasis = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Metastasis' & IMPACT_data$CANCER_TYPE == i)])
  
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


##-----------------
## Minority CancerTypes
##-----------------
IMPACT_data = IMPACT_data[which(IMPACT_data$Y_CNA == TRUE), ]
IMPACT_minority = table(IMPACT_data$CANCER_TYPE)
IMPACT_minority = IMPACT_minority[order(IMPACT_minority, decreasing = T)]
IMPACT_minority = IMPACT_minority[21:36]

IMPACT_minority_df = data.frame()
for(i in unique(names(IMPACT_minority))){
  n.all = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE %in% c('Primary', 'Metastasis') & IMPACT_data$CANCER_TYPE == i)])
  n.primary = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Primary' & IMPACT_data$CANCER_TYPE == i)])
  n.metastasis = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Metastasis' & IMPACT_data$CANCER_TYPE == i)])
  
  primary.frequency = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$SAMPLE_TYPE == 'Primary' & IMPACT_data$CANCER_TYPE == i)]) / n.primary
  metastatic.frequency = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$SAMPLE_TYPE == 'Metastasis' & IMPACT_data$CANCER_TYPE == i)]) / n.metastasis
  cancer_median = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$CANCER_TYPE == i)]) / n.all
  cancer_summary = data.frame(CancerType = i,
                              Type = c('Primary', 'Metastasis'),
                              value = c(primary.frequency, metastatic.frequency),
                              median_cancerType = cancer_median)
  
  IMPACT_minority_df = rbind(IMPACT_minority_df, cancer_summary)
  rm(n.all, cancer_summary)
  rm(primary.frequency)
  rm(metastatic.frequency)
  rm(cancer_median)
}

#' add sample amount
IMPACT_minority_df$mergedCancerType = NA
for(i in unique(names(IMPACT_minority))){
  IMPACT_minority_df$mergedCancerType[which(IMPACT_minority_df$CancerType == i)] = IMPACT_minority[which(names(IMPACT_minority) == i)][[1]]
}

IMPACT_minority_df$mergedCancerType = paste0(IMPACT_minority_df$CancerType, ' (n=', IMPACT_minority_df$mergedCancerType, ')')


##-----------------
## Minority Plot
##-----------------
IMPACT_minority_plot = ggplot(IMPACT_minority_df, 
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

IMPACT_minority_plot

ggsave_golden(filename = 'Figures_original/Y_loss_SITE_minorities.pdf', plot = IMPACT_minority_plot, width = 6)



##-----------------
## Association with AGE
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





##-----------------
## Old misc scripts
##-----------------
#' cohortData = readRDS('Data_out/cohort_data.rds')
#' IMPACT_data = cohortData$IMPACT.cohort
#' WES_data = cohortData$WES.cohort
#' TCGA_data = cohortData$male_TCGA_cohort
#' 
#' #' make overview of available data resources (IMPACT, WES and TCGA)
#' #' IMPACT
#' top_20_cancer_type = table(IMPACT_data$CANCER_TYPE)
#' top_20_cancer_type = top_20_cancer_type[order(top_20_cancer_type, decreasing = T)]
#' top_20_cancer_type = top_20_cancer_type[1:20]
#' 
#' pdf(paste0('Figures/IMPACT_Top20_cancer_type.pdf'), width = 9, height = 6)
#' par(oma = c(0,10,0,0))
#' barplot(top_20_cancer_type[20:1], 
#'         horiz = T, 
#'         border = NA, 
#'         axes = F, 
#'         las = 2, 
#'         xlim = c(0, 1800), 
#'         cex.axis = 0.75,
#'         main = paste0('MSK-IMPACT (n = ', dim(IMPACT_data)[1], ')'))
#' axis(side = 1, at = c(0,500,1000,1500,2000), labels = c(0,500,1000,1500,2000))
#' dev.off()
#' 
#' #' WES
#' top_20_cancer_type = table(WES_data$Parental_Tumor_Type)
#' top_20_cancer_type = top_20_cancer_type[order(top_20_cancer_type, decreasing = T)]
#' top_20_cancer_type = top_20_cancer_type[1:20]
#' 
#' pdf(paste0('Figures/WES_Top20_cancer_type.pdf'), width = 9, height = 6)
#' par(oma = c(0,10,0,0))
#' barplot(top_20_cancer_type[20:1], 
#'         horiz = T, 
#'         border = NA, 
#'         axes = F, 
#'         las = 2, 
#'         xlim = c(0, 150), 
#'         cex.axis = 0.75,
#'         main = paste0('MSK-WES (n = ', dim(WES_data)[1], ')'))
#' axis(side = 1, at = c(0, 50, 100, 150), labels = c(0, 50, 100, 150))
#' dev.off()
#' 
#' #' TCGA
#' top_20_cancer_type = table(TCGA_data$type)
#' top_20_cancer_type = top_20_cancer_type[order(top_20_cancer_type, decreasing = T)]
#' top_20_cancer_type = top_20_cancer_type[1:20]
#' 
#' pdf(paste0('Figures/TCGA_Top20_cancer_type.pdf'), width = 9, height = 6)
#' par(oma = c(0,10,0,0))
#' barplot(top_20_cancer_type[20:1], 
#'         horiz = T, 
#'         border = NA, 
#'         axes = F, 
#'         las = 2, 
#'         xlim = c(0, 500), 
#'         cex.axis = 0.75,
#'         main = paste0('TCGA-WES (n = ', dim(TCGA_data)[1], ')'))
#' axis(side = 1, at = seq(0, 500, 100), labels = seq(0, 500, 100))
#' dev.off()
#' 
#' 
#' 
#' #' Y-chromosome loss among different cancer types; stratified according to primary and metastatic samples;
#' #' start with IMPACT
#' IMPACT_data = cohortData$IMPACT.cohort
#' IMPACT_top_20 = table(IMPACT_data$CANCER_TYPE)
#' IMPACT_top_20 = IMPACT_top_20[order(IMPACT_top_20, decreasing = T)]
#' IMPACT_top_20 = IMPACT_top_20[1:20]
#' 
#' IMPACT_incidence = data.frame()
#' for(i in unique(names(IMPACT_top_20))){
#'   n.all = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE %in% c('Primary', 'Metastasis') & IMPACT_data$CANCER_TYPE == i)])
#'   n.primary = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Primary' & IMPACT_data$CANCER_TYPE == i)])
#'   n.metastasis = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Metastasis' & IMPACT_data$CANCER_TYPE == i)])
#'   
#'   primary.frequency = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$SAMPLE_TYPE == 'Primary' & IMPACT_data$CANCER_TYPE == i)]) / n.primary
#'   metastatic.frequency = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$SAMPLE_TYPE == 'Metastasis' & IMPACT_data$CANCER_TYPE == i)]) / n.metastasis
#'   cancer_median = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$CANCER_TYPE == i)]) / n.all
#'   cancer_summary = data.frame(CancerType = i,
#'                               Type = c('Primary', 'Metastasis'),
#'                               value = c(primary.frequency, metastatic.frequency),
#'                               median_cancerType = cancer_median)
#'   
#'   IMPACT_incidence = rbind(IMPACT_incidence, cancer_summary)
#'   rm(n.all)
#'   rm(primary.frequency)
#'   rm(metastatic.frequency)
#'   rm(cancer_median)
#' }
#' 
#' #' add sample amount
#' IMPACT_incidence$mergedCancerType = NA
#' for(i in unique(names(IMPACT_top_20))){
#'   IMPACT_incidence$mergedCancerType[which(IMPACT_incidence$CancerType == i)] = IMPACT_top_20[which(names(IMPACT_top_20) == i)][[1]]
#' }
#' 
#' IMPACT_incidence$mergedCancerType = paste0(IMPACT_incidence$CancerType, ' (n=', IMPACT_incidence$mergedCancerType, ')')
#' write.table(IMPACT_incidence, file = 'Data_out/IMPACT/IMPACT_Y.loss_incidence.txt', sep = '\t', row.names = F, quote = F)
#' 
#' #' make the visualization
#' #' color code from Subhi:
#' color_selected = brewer.pal(n = 6, name = 'Paired')
#' color_selected = color_selected[c(5,6,1,2)]
#' 
#' IMPACT_incidence_plot = ggplot(IMPACT_incidence, 
#'                                aes(x = reorder(mergedCancerType, median_cancerType), 
#'                                    y = value, 
#'                                    fill = Type)) + 
#'                                  
#'   geom_bar(stat = 'identity', position = position_dodge(0.82), width = 0.8) +
#'   scale_fill_manual(values = c('Primary' = color_selected[4],
#'                                'Metastasis' = color_selected[2]),
#'                     guide = guide_legend(direction = 'horizontal',
#'                                          title = 'Site',
#'                                          label.theme = element_text(size = 14))) +
#'   scale_y_continuous(expand = c(0.005, 0)) +
#'   geom_hline(yintercept = seq(0, 0.6, 0.2), color = 'grey35', linetype = 'dashed', size = 0.2) +
#'   theme(legend.position = 'top',
#'         panel.background = element_rect(fill = NA),
#'         axis.ticks.y = element_blank(),
#'         axis.ticks.x = element_line(size = 0.2),
#'         axis.line.x = element_line(size = 0.4),
#'         axis.text = element_text(size = 10, color = 'black')) +
#'   coord_flip() +
#'   labs(y = '% Y chromosome loss', x = '')
#'   
#' IMPACT_incidence_plot  
#' 
#' ggsave_golden(IMPACT_incidence_plot, filename = 'Figures/IMPACT_Y_Loss.pdf', width = 9)
#' 
#' 
#' #' look at the overall distribution of Y-chromosome loss across tissue site
#' IMPACT_incidence$Type = factor(IMPACT_incidence$Type, levels = c('Primary', 'Metastasis'))
#' Incidence_site = ggplot(IMPACT_incidence, aes(x = Type, y = value, color = Type)) +
#'   geom_boxplot(width = 0.4, size = 0.85) +
#'   geom_jitter(width = 0.15) +
#'   scale_color_manual(values = c('Primary' = color_selected[4],
#'                                 'Metastasis' = color_selected[2])) +
#'   scale_y_continuous(expand = c(0, 0),
#'                      limits = c(0, 0.8)) +
#'   theme_Y + 
#'   theme(legend.position = 'none',
#'         aspect.ratio = 2.2,
#'         axis.text.x = element_text(angle = 45, hjust = 1)) +
#'   labs(x = '', y = 'Y-chromosome loss')
#' 
#' ggsave_golden(filename = 'Figures/Y_loss_SITE.pdf', plot = Incidence_site, width = 6)
#' 
#' #' statistics
#' t.test(IMPACT_incidence$value ~ IMPACT_incidence$Type)
#' 
#' 
#' ###############################################################################
#' ###############################################################################
#' #' look at the age distribution, and whether age may be associated with higher Y chromosome loss
#' IMPACT_data # raw
#' IMPACT_data$AGE_AT_SEQ_REPORTED_YEARS = as.numeric(as.character(IMPACT_data$AGE_AT_SEQ_REPORTED_YEARS))
#' 
#' plot_list = list()
#' for(i in unique(names(IMPACT_top_20))){
#'   data.plot = IMPACT_data[which(IMPACT_data$CANCER_TYPE == i), ]
#'   type.density = density(data.plot$AGE_AT_SEQ_REPORTED_YEARS, na.rm = T)
#'   median.age = type.density$x[which.max(type.density$y)]
#'   print(median.age)
#'   peak.max = type.density$y[which.max(type.density$y)]
#'   #median.age = median(data.plot$AGE_AT_SEQ_REPORTED_YEARS, na.rm = T)
#'   upper.limit.plot = peak.max / 15
#'   
#'   plot = ggplot(data.plot, aes(x = AGE_AT_SEQ_REPORTED_YEARS)) + 
#'     geom_density(size = 1, bw = 'SJ') +
#'     theme_minimal() +
#'     theme(panel.grid = element_blank(),
#'           axis.title = element_blank(),
#'           axis.text.y = element_blank(),
#'           aspect.ratio = 1) +
#'     scale_y_continuous(expand = c(0.01, 0)) +
#'     scale_x_continuous(limits = c(30, 90)) +
#'     geom_segment(x = median.age, y = 0, xend = median.age, yend = upper.limit.plot, color = 'red', size = 0.8) +
#'     labs(title = i)
#'   
#'   plot_list[[i]] = plot
#'   
#'   rm(data.plot,
#'      median.age,
#'      peak.max,
#'      upper.limit.plot)
#' }
#' 
#' Age_distribution_IMPACT = plot_list$`Esophagogastric Cancer` / plot_list$`Pancreatic Cancer` / plot_list$`Colorectal Cancer` / plot_list$`Soft Tissue Sarcoma` / plot_list$`Prostate Cancer` / plot_list$Glioma
#' ggsave_golden(plot = Age_distribution_IMPACT, filename = 'Figures/Age_Distribution_IMPACT.pdf', width = 12)
#' 
#' 
#' #' look at the general trend between Y-chromosome loss and intact among age groups;
#' i_without = IMPACT_data[!is.na(IMPACT_data$AGE_AT_SEQ_REPORTED_YEARS) & !is.na(IMPACT_data$Y_call), ]
#' 
#' Age_cohort = ggplot(i_without) + 
#'   geom_histogram(aes(x = AGE_AT_SEQ_REPORTED_YEARS, 
#'                      color = Y_call, 
#'                      fill = Y_call), bins = 35, binwidth = 0.999, na.rm = T) +
#'   scale_y_continuous(expand = c(0.01, 0.01)) +
#'   scale_x_continuous(expand = c(0.01, 0.01),
#'                      limits = c(0, 90),
#'                      breaks = seq(0, 90, 30)) +
#'   scale_fill_manual(values = c('Y_chrom_loss' = color_selected[2],
#'                              'intact_Y_chrom' = color_selected[4]),
#'                     labels = c('Y loss', 'Intact Y chromosome'),
#'                     guide = guide_legend(direction = 'horizontal',
#'                                          title = '',
#'                                          label.theme = element_text(size = 12))) +
#'   scale_color_manual(values = c('Y_chrom_loss' = color_selected[2],
#'                                 'intact_Y_chrom' = color_selected[4]), guide = 'none') +
#'   theme_Y +
#'   theme(axis.text.y = element_blank(),
#'         legend.position = 'top',
#'         aspect.ratio = 1) +
#'   labs(x = 'Age', y = 'Counts')
#'         
#' ggsave_golden(Age_cohort, filename = 'Figures/Age_Distribution_Y_loss.pdf', width = 7)
#' 
#' 
#' #' out