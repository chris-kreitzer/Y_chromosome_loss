## Look into Y-chromosome incidences among different cancer types
## Make an overview graph; as Bastien initially did (back in 2019)
## 
## 09/01/2021
## chris-kreitzer
## 

set.seed(112)
rm(list = ls())
.rs.restartR()

## Libraries and Input
cohortData = readRDS('Data_out/cohort_data.rds')
IMPACT_data = cohortData$IMPACT.cohort
WES_data = cohortData$WES.cohort
TCGA_data = cohortData$male_TCGA_cohort

#' make overview of available data resources (IMPACT, WES and TCGA)
#' IMPACT
top_20_cancer_type = table(IMPACT_data$CANCER_TYPE)
top_20_cancer_type = top_20_cancer_type[order(top_20_cancer_type, decreasing = T)]
top_20_cancer_type = top_20_cancer_type[1:20]

pdf(paste0('Figures/IMPACT_Top20_cancer_type.pdf'), width = 9, height = 6)
par(oma = c(0,10,0,0))
barplot(top_20_cancer_type[20:1], 
        horiz = T, 
        border = NA, 
        axes = F, 
        las = 2, 
        xlim = c(0, 1800), 
        cex.axis = 0.75,
        main = paste0('MSK-IMPACT (n = ', dim(IMPACT_data)[1], ')'))
axis(side = 1, at = c(0,500,1000,1500,2000), labels = c(0,500,1000,1500,2000))
dev.off()

#' WES
top_20_cancer_type = table(WES_data$Parental_Tumor_Type)
top_20_cancer_type = top_20_cancer_type[order(top_20_cancer_type, decreasing = T)]
top_20_cancer_type = top_20_cancer_type[1:20]

pdf(paste0('Figures/WES_Top20_cancer_type.pdf'), width = 9, height = 6)
par(oma = c(0,10,0,0))
barplot(top_20_cancer_type[20:1], 
        horiz = T, 
        border = NA, 
        axes = F, 
        las = 2, 
        xlim = c(0, 150), 
        cex.axis = 0.75,
        main = paste0('MSK-WES (n = ', dim(WES_data)[1], ')'))
axis(side = 1, at = c(0, 50, 100, 150), labels = c(0, 50, 100, 150))
dev.off()

#' TCGA
top_20_cancer_type = table(TCGA_data$type)
top_20_cancer_type = top_20_cancer_type[order(top_20_cancer_type, decreasing = T)]
top_20_cancer_type = top_20_cancer_type[1:20]

pdf(paste0('Figures/TCGA_Top20_cancer_type.pdf'), width = 9, height = 6)
par(oma = c(0,10,0,0))
barplot(top_20_cancer_type[20:1], 
        horiz = T, 
        border = NA, 
        axes = F, 
        las = 2, 
        xlim = c(0, 500), 
        cex.axis = 0.75,
        main = paste0('TCGA-WES (n = ', dim(TCGA_data)[1], ')'))
axis(side = 1, at = seq(0, 500, 100), labels = seq(0, 500, 100))
dev.off()



#' Y-chromosome loss among different cancer types; stratified according to primary and metastatic samples;
#' start with IMPACT
IMPACT_data = cohortData$IMPACT.cohort
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
  rm(n.all)
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
write.table(IMPACT_incidence, file = 'Data_out/IMPACT/IMPACT_Y.loss_incidence.txt', sep = '\t', row.names = F, quote = F)

#' make the visualization
IMPACT_incidence_plot = ggplot(IMPACT_incidence, 
                               aes(x = reorder(mergedCancerType, median_cancerType), 
                                   y = value, 
                                   fill = Type)) + 
                                 
  geom_bar(stat = 'identity', position = position_dodge(0.82), width = 0.8) +
  scale_fill_manual(values = c('Primary' = '#c31f21',
                               'Metastasis' = '#4977b0'),
                    guide = guide_legend(direction = 'horizontal',
                                         title = 'Site',
                                         label.theme = element_text(size = 12))) +
  scale_y_continuous(expand = c(0.005, 0)) +
  geom_hline(yintercept = seq(0, 0.6, 0.2), color = 'grey85', linetype = 'dashed', size = 0.2) +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.line.x = element_line(size = 0.4),
        axis.text = element_text(size = 10, color = 'black')) +
  coord_flip() +
  labs(y = '% Y chromosome loss', x = '')
  
IMPACT_incidence_plot  

ggsave_golden(IMPACT_incidence_plot, filename = 'Figures/IMPACT_Y_Loss.pdf', width = 9)


###############################################################################
###############################################################################
#' look at the age distribution, and whether age may be associated with higher Y chromosome loss
IMPACT_data # raw
IMPACT_data$AGE_AT_SEQ_REPORTED_YEARS = as.numeric(as.character(IMPACT_data$AGE_AT_SEQ_REPORTED_YEARS))

plot_list = list()
for(i in unique(names(IMPACT_top_20))){
  data.plot = IMPACT_data[which(IMPACT_data$CANCER_TYPE == i), ]
  type.density = density(data.plot$AGE_AT_SEQ_REPORTED_YEARS, na.rm = T)
  median.age = type.density$x[which.max(type.density$y)]
  print(median.age)
  peak.max = type.density$y[which.max(type.density$y)]
  #median.age = median(data.plot$AGE_AT_SEQ_REPORTED_YEARS, na.rm = T)
  upper.limit.plot = peak.max / 15
  
  plot = ggplot(data.plot, aes(x = AGE_AT_SEQ_REPORTED_YEARS)) + 
    geom_density(size = 1, bw = 'SJ') +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          aspect.ratio = 1) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(30, 90)) +
    geom_segment(x = median.age, y = 0, xend = median.age, yend = upper.limit.plot, color = 'red', size = 0.8) +
    labs(title = i)
  
  plot_list[[i]] = plot
  
  rm(data.plot,
     median.age,
     peak.max,
     upper.limit.plot)
}

Age_distribution_IMPACT = plot_list$`Esophagogastric Cancer` / plot_list$`Pancreatic Cancer` / plot_list$`Colorectal Cancer` / plot_list$`Soft Tissue Sarcoma` / plot_list$`Prostate Cancer` / plot_list$Glioma
ggsave_golden(plot = Age_distribution_IMPACT, filename = 'Figures/Age_Distribution_IMPACT.pdf', width = 12)


#' look at the general trend between Y-chromosome loss and intact among age groups;
i_without = IMPACT_data[!is.na(IMPACT_data$AGE_AT_SEQ_REPORTED_YEARS), ]
Age_cohort = ggplot(i_without) + 
  geom_histogram(aes(x = AGE_AT_SEQ_REPORTED_YEARS, color = Y_call, fill = Y_call), bins = 35, binwidth = 0.999, na.rm = T) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 90),
                     breaks = seq(0, 90, 30)) +
  scale_fill_manual(values = c('Y_chrom_loss' = '#c31f21',
                               'intact_Y_chrom' = '#4977b0'),
                    labels = c('Y loss', 'Intact Y chromosome'),
                    guide = guide_legend(direction = 'horizontal',
                                         title = '',
                                         label.theme = element_text(size = 12))) +
  scale_color_manual(values = c('Y_chrom_loss' = '#c31f21',
                               'intact_Y_chrom' = '#4977b0')) +
  guides(color = FALSE) +
                     
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_line(size = 0.4, color = 'black'),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, color = 'black'),
        panel.border = element_rect(fill = NA, size = 0.6),
        legend.position = c(0.1, 0.9)) +
  labs(x = 'Age')
  
ggsave_golden(Age_cohort, filename = 'Figures/Age_Distribution_Y_loss.pdf', width = 9)


#' out