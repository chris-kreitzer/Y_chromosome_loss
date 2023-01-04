##-----------------
## IGV plots on Y-chromosome
## cancers; modified
##-----------------

clean()
gc()
.rs.restartR()
setwd(dir = '~/Documents/MSKCC/10_MasterThesis/')
source('Scripts/IGV-like_plot.R')
ggsave_golden = function(filename, plot, width, ...){
  ggsave(filename = filename, plot = plot, device = cairo_pdf, width = width, height = width / 1.61803398875)
}

library(tidyverse)
library(tidyr)
library(dplyr)
library(plyr)
library(facetsSuite)
library(cowplot)

##----------------- 
## Data
IGV = read.csv('Data/04_Loss/IMPACT_IGV_out.seg', sep = '\t')
data = readRDS('Data/signedOut/Cohort_07132022.rds')
CNA = data$IMPACT_Y_classification_final
clinical = data$IMPACT_clinicalAnnotation

cancers_keep = names(sort(table(clinical$CANCER_TYPE), decreasing = T))[1:20]

## make a for loop over the top20 cancers
plot_list = list()
for(i in unique(cancers_keep)){
  samples = clinical$SAMPLE_ID[which(clinical$CANCER_TYPE == i)]
  CNA_samples = CNA[which(CNA$sample %in% samples), ]
  IGV_samples = IGV[which(IGV$ID %in% CNA_samples$sample), ]
  data_merged = merge(IGV_samples, CNA_samples, by.x = 'ID', by.y = 'sample', all.x = T)
  print(paste0(i, ' ', length(unique(data_merged$ID))))
  
  #' define categories:
  loss = unique(data_merged$ID[which(data_merged$classification == 'loss')])
  relative_loss = unique(data_merged$ID[which(data_merged$classification == 'relative_loss')])
  wt = unique(data_merged$ID[which(data_merged$classification == 'wt')])
  gain_loss = unique(data_merged$ID[which(data_merged$classification == 'gain_loss')])
  gain = unique(data_merged$ID[which(data_merged$classification == 'gain')])
  
  data_plot = data_merged %>%
    arrange(factor(ID, levels = c(loss, relative_loss, wt, gain_loss, gain)))
  
  plot_IGV = plot_segmentation(data_plot, cap_log_ratios = 2)
  plot_line = plot_IGV[[2]] +
    geom_hline(yintercept = length(loss) * plot_IGV[[1]], size = 0.4, linetype = 'dashed', color = 'grey15') +
    geom_hline(yintercept = sum(length(loss) * plot_IGV[[1]], length(relative_loss) * plot_IGV[[1]]), size = 0.4, linetype = 'dashed', color = 'grey15') +
    geom_hline(yintercept = sum(length(loss) * plot_IGV[[1]], length(relative_loss) * plot_IGV[[1]], length(wt) * plot_IGV[[1]]), size = 0.4, linetype = 'dashed', color = 'grey15') +
    geom_hline(yintercept = sum(length(loss) * plot_IGV[[1]], length(relative_loss) * plot_IGV[[1]], 
                                length(wt) * plot_IGV[[1]], length(gain_loss) * plot_IGV[[1]]), size = 0.4, linetype = 'dashed', color = 'grey15') +
    labs(title = i)
  
  ggsave_golden(filename = paste0('Data/05_Association/IGV_plots/', i, '.pdf'), plot = plot_line, width = 9)
  plot_list[[i]] = plot_line
  
}


##-----------------
## Further Associations:
##-----------------