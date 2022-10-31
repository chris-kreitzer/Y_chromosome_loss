clean()
.rs.restartR()

library(tidyverse)
library(tidyr)
library(dplyr)
library(plyr)
library(facetsSuite)

plot_segmentation = function(segs,
                             all_chromosomes = TRUE,
                             sample_order = NULL,
                             cap_log_ratios = TRUE,
                             colors = c('darkblue', 'white', 'red'),
                             return_object = FALSE){
  
  if (all_chromosomes) { 
    segs = filter(segs, chrom %in% seq(1,24,1))
    chrom_order = seq(1:24)
  } else {
    chrom_order = seq(1:23)
  }
  segs = mutate(segs, chrom = factor(chrom, levels = chrom_order, ordered = T))
  
  if (!is.null(sample_order)) {
    if (!all(segs$ID %in% sample_order)) stop('Samples missing from provided sample order', call. = FALSE)
    #segs = mutate(segs, ID = factor(ID, levels = sample_order, ordered = T))
    segs = segs[order(factor(segs$ID, levels=unique(sample_order))),]
  }
  
  # Determine lengths of chromosomes and adjust X-axis accordingly
  chrom_lengths = map_dfr(unique(segs$chrom), function(x) {
    chrom_max = max(segs$loc.end[segs$chrom == x], na.rm = T)
    chrom_min = min(segs$loc.end[segs$chrom == x], na.rm = T)
    list(chrom = x,
         chrom_length = as.numeric(chrom_max - chrom_min))
  }) %>%
    mutate(rel_length = chrom_length/sum(chrom_length)) # Edit, get rid of .data
  segs = left_join(segs, chrom_lengths, by = 'chrom')
  
  # Cap log-ratios and set colorlegend
  if (cap_log_ratios != FALSE) {
    if (is.numeric(cap_log_ratios)) {
      segs$seg.mean[which(segs$seg.mean > cap_log_ratios)] = cap_log_ratios
      segs$seg.mean[which(segs$seg.mean < -cap_log_ratios)] = -cap_log_ratios
      legend_breaks = c(-cap_log_ratios, -cap_log_ratios/2, 0, cap_log_ratios/2, cap_log_ratios)
    } else {
      segs$seg.mean[which(segs$seg.mean > 1)] = 1
      segs$seg.mean[which(segs$seg.mean < -1)] = -1
      legend_breaks = seq(-1, 1, 0.4)
    }
  } else {
    legend_breaks = c(min(segs$seg.mean), min(segs$seg.mean)/2, 0, max(segs$seg.mean)/2, max(segs$seg.mean))
  }
  
  # Set Y-axis
  sample_number = length(unique(segs$ID))
  increments = 100 / sample_number
  print(increments)
  max_vec = cumsum(rep(increments, sample_number))
  min_vec = c(0, max_vec[-length(max_vec)])
  segs = mutate(segs,
                y_min = min_vec[match(ID, factor(unique(segs$ID)))],
                y_max = max_vec[match(ID, factor(unique(segs$ID)))])
  
  y_labs = distinct(segs, ID, .keep_all = T) %>%
    mutate(pos = (y_max-y_min)/2 + y_min)
  
  # Plot seg files
  seg_plot = ggplot(segs, aes(xmin = loc.start, xmax = loc.end, 
                              ymin = y_min, ymax = y_max, 
                              fill = seg.mean)) +
    geom_rect() +
    scale_fill_gradient2(low = colors[1], mid = colors[2], high = colors[3], guide = 'colourbar',
                         'Log ratio', breaks = legend_breaks, labels = round(legend_breaks)) +
    scale_y_continuous(expand = c(0,0), breaks = y_labs$pos, labels = y_labs$ID) +
    facet_grid(.~chrom, space = 'free_x', scales = 'free_x', switch = 'x') +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.spacing.x = unit(-.5, 'lines'),
      panel.spacing.y = unit(0, 'lines'),
      strip.background = element_rect(fill = 'white'),
      panel.background = element_rect(fill = 'white'),
      legend.position = "bottom"
    ) + 
    geom_vline(aes(xintercept = 0), colour = "black", size = 0.4, data = subset(segs, chrom == 1)) +
    geom_hline(aes(yintercept = 0), colour = "black", size = 0.4) +
    geom_hline(aes(yintercept = 100), colour = "black", size = 0.4)
  
  if (return_object == TRUE) {
    seg_plot
  } else {
    print(seg_plot)
  }
}

#' IGV-like plot:
plot_segmentation(segs = IGV, plotX = F, cap_log_ratios = 3)





##------ tests
setwd('~/Documents/MSKCC/10_MasterThesis/')
IGV = read.csv('Data/04_Loss/IMPACT_IGV_out.seg', sep = '\t')
data = readRDS('Data/signedOut/Cohort_07132022.rds')
CNA = data$IMPACT_Y_classification_final
clinical = data$IMPACT_clinicalAnnotation

##-----
## Esophagogastric cancer
eso = clinical$SAMPLE_ID[which(clinical$CANCER_TYPE == 'Esophagogastric Cancer')]
CNA_eso = CNA[which(CNA$sample %in% eso), ]
IGV_eso = IGV[which(IGV$ID %in% CNA_eso$sample), ]
IGV_esophagastric = merge(IGV_eso, CNA_eso, by.x = 'ID', by.y = 'sample', all.x = T)





## 
sa = base::sample(size = 10, x = IGV_esophagastric$ID, replace = F)
test = IGV_esophagastric[which(IGV_esophagastric$ID %in% sa), ]

unique(IGV_esophagastric$classification)
loss = unique(IGV_esophagastric$ID[which(IGV_esophagastric$classification == 'loss')])
gain = unique(IGV_esophagastric$ID[which(IGV_esophagastric$classification == 'gain')])
wt = unique(IGV_esophagastric$ID[which(IGV_esophagastric$classification == 'wt')])
relative_loss = unique(IGV_esophagastric$ID[which(IGV_esophagastric$classification == 'relative_loss')])
gain_loss = unique(IGV_esophagastric$ID[which(IGV_esophagastric$classification == 'gain_loss')])

df = IGV_esophagastric %>%
  arrange(factor(ID, levels = c(loss, relative_loss, wt, gain_loss, gain)))




plot_segmentation(df, cap_log_ratios = 2) +
  geom_hline(yintercept = length(unique(loss)) * 0.1360544) +
  geom_hline(yintercept = sum(length(loss) * 0.1360544, length(relative_loss) * 0.1360544)) +
  geom_hline(yintercept = sum(length(loss) * 0.1360544, length(relative_loss) * 0.1360544, length(wt) * 0.1360544)) +
  geom_hline(yintercept = sum(length(loss) * 0.1360544, length(relative_loss) * 0.1360544, length(wt) * 0.1360544, 
                              length(relative_loss) * 0.1360544))




##---
## stratify according to Y_chrom loss and order the patients
## within



IGV = IGV[which(IGV$ID %in% Clinical$SAMPLE_ID[which(Clinical$CANCER_TYPE == 'Bladder Cancer')]), ]
length(unique(IGV$ID))
n = sample(x = unique(IGV$ID), size = 10, replace = )

x = IGV[which(IGV$ID %in% n), ]
plot_segmentation(segs = x, cap_log_ratios = 2, plotX = F) +
  geom_hline(yintercept = 10) +
  geom_hline(yintercept = 20, linetype = 'dashed', size = 2)







