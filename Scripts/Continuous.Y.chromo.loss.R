## PanCancer analysis of Y-chromosome loss:
## qualitatively call Y-chromosome loss from WES samples (binary state; 0/1)
## Transform qualitative measure (50% rule) into continious scale via median Copy number log Ratio
## starting with Prostate samples
## 
## start analysis: 16/03/2021
## start refinement: 08/04/2021 (based on discussion with Niki and Subhi)
## output data fetched from Juno (run via bsub)


rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')
source('Y_chromosome_loss/Plotting_theme.R')

# Libraries
library(stringi)
library(ggplot2)
library(cowplot)
library(egg)


# Input:
# Processed on Juno: 08. April
# Whole exame data analysis: Prostate
WES.prostate.processed_df = read.csv('Data_out/WES.Prostate.processed.out.15.4.txt', sep = '\t')
WES.prostate.binary_df = read.csv('Data_out/WES.Prostate.binary.15.4.txt', sep = '\t')
WES.prostate.QC_df = read.csv('Data_out/WES.prostate.QC.15.4.txt', sep = '\t')
WES.prostate.cnlr_df = read.csv('Data_out/Cnlr.WES.prostate.out.15.4.txt', sep = '\t')

#' FacetsY for IMPACT:
IMPACT.prostate.processed_df = read.csv('Data_out/IMPACT.Prostate.processed.out.16.4.txt', sep = '\t')
IMPACT.prostate.binary_df = read.csv('Data_out/IMPACT.Prostate.binary.16.4.txt', sep = '\t')
IMPACT.prostate.QC_df = read.csv('Data_out/IMPACT.prostate.QC.16.4.txt', sep = '\t')
IMPACT.prostate.cnlr_df = read.csv('Data_out/Cnlr.IMPACT.prostate.out.16.4.txt', sep = '\t')

#' WES to IMPACT ID matches
WES.IMPACT.ids_df = read.csv('Data/WES_cBio_ID.match.txt', sep = '\t')
WES.prostate.samples.out = substr(unique(WES.prostate.processed_df$ID), start = 3, 
                                  stop = 19)


## Functions:



## Processing
WES.IMPACT.ids_df$WES.ID = NA
for(i in 1:nrow(WES.IMPACT.ids_df)){
  WES.IMPACT.ids_df$WES.ID[i] = stringr::str_split(string = WES.IMPACT.ids_df$Facet_Path, pattern = '/')[[i]][1]
}

WES.IMPACT.ids_df = WES.IMPACT.ids_df[, c('DMP_Sample_ID', 'Sex', 'Tumor_Type', 
                                      'Parental_Tumor_Type', 'Sample_Class', 'WES.ID')]
WES.IMPACT.ids_df$WES.ID = substr(WES.IMPACT.ids_df$WES.ID, start = 1, stop = 17)

WES.prostate.samples.retain = WES.IMPACT.ids_df$DMP_Sample_ID[which(WES.IMPACT.ids_df$WES.ID %in% WES.prostate.samples.out)]

#' modify sample IDs in IMPACT cohort
IMPACT.prostate.binary_df$sample.id = paste0(substr(IMPACT.prostate.binary_df$sample.id, start = 8, stop = 16),
                                     '-T', 
                                     substr(IMPACT.prostate.binary_df$sample.id, start = 1, stop = 6))
IMPACT.prostate.binary_df = IMPACT.prostate.binary_df[!duplicated(IMPACT.prostate.binary_df$sample.id), ]

IMPACT.prostate.cnlr_df$ID = paste0(substr(IMPACT.prostate.cnlr_df$ID, start = 8, stop = 16),
                                    '-T', 
                                    substr(IMPACT.prostate.cnlr_df$ID, start = 1, stop = 6))
IMPACT.prostate.processed_df$ID = paste0(substr(IMPACT.prostate.processed_df$ID, start = 8, stop = 16),
                                         '-T', 
                                         substr(IMPACT.prostate.processed_df$ID, start = 1, stop = 6))
IMPACT.prostate.QC_df$ID = paste0(substr(IMPACT.prostate.QC_df$ID, start = 8, stop = 16),
                                  '-T', 
                                  substr(IMPACT.prostate.QC_df$ID, start = 1, stop = 6))
IMPACT.prostate.QC_df = IMPACT.prostate.QC_df[!duplicated(IMPACT.prostate.QC_df$ID), ]

#' modify sample IDs in WES cohort
WES.prostate.binary_df$sample.id = substr(WES.prostate.binary_df$sample.id, start = 3, stop = 19)
WES.prostate.cnlr_df$ID = substr(WES.prostate.cnlr_df$ID, start = 3, stop = 19)
WES.prostate.processed_df$ID = substr(WES.prostate.processed_df$ID, start = 3, stop = 19)
WES.prostate.QC_df$ID = substr(WES.prostate.QC_df$ID, start = 3, stop = 19)

## Analysis:



#' match WES and Panel results:
DMP.CMO.match_df = WES.IMPACT.ids_df[which(WES.IMPACT.ids_df$DMP_Sample_ID %in% WES.prostate.samples.retain), 
                                     c('DMP_Sample_ID', 'WES.ID'), drop = F]

DMP.CMO.match_df$W.50_call = NA
DMP.CMO.match_df$W.QC = NA
DMP.CMO.match_df$I.QC = NA
DMP.CMO.match_df$W.WGD = NA
DMP.CMO.match_df$I.WGD = NA
DMP.CMO.match_df$I.50_call = NA
DMP.CMO.match_df$W.Cnlr_segment = NA
DMP.CMO.match_df$I.Cnlr_segment = NA
DMP.CMO.match_df$W.mCnlr = NA
DMP.CMO.match_df$I.mCnlr = NA
DMP.CMO.match_df$W.mCnlr.sd = NA
DMP.CMO.match_df$I.mCnlr.sd = NA
DMP.CMO.match_df$W.purity = NA
DMP.CMO.match_df$I.purity = NA
DMP.CMO.match_df$W.ploidy = NA
DMP.CMO.match_df$I.ploidy = NA

#' loop through every sample and fetch information:
for(i in 1:nrow(DMP.CMO.match_df)){
  #' WES data
  DMP.CMO.match_df$W.50_call[i] = WES.prostate.binary_df$Y_call[which(WES.prostate.binary_df$sample.id == DMP.CMO.match_df$WES.ID[i])]
  DMP.CMO.match_df$W.QC[i] = WES.prostate.QC_df$QC[which(WES.prostate.QC_df$ID == DMP.CMO.match_df$WES.ID[i])]
  DMP.CMO.match_df$W.WGD[i] = WES.prostate.QC_df$wgd[which(WES.prostate.QC_df$ID == DMP.CMO.match_df$WES.ID[i])]
  DMP.CMO.match_df$W.ploidy[i] = WES.prostate.binary_df$ploidy[which(WES.prostate.binary_df$sample.id == DMP.CMO.match_df$WES.ID[i])]
  DMP.CMO.match_df$W.purity[i] = WES.prostate.binary_df$purity[which(WES.prostate.binary_df$sample.id == DMP.CMO.match_df$WES.ID[i])]
  DMP.CMO.match_df$W.Cnlr_segment[i] = median(WES.prostate.processed_df$cnlr.median[which(WES.prostate.processed_df$ID == DMP.CMO.match_df$WES.ID[i] &
                                                                                            WES.prostate.processed_df$chrom == 24)], na.rm = T)
  DMP.CMO.match_df$W.mCnlr[i] = median(WES.prostate.cnlr_df$cnlr[which(WES.prostate.cnlr_df$ID == DMP.CMO.match_df$WES.ID[i])], na.rm = T)
  DMP.CMO.match_df$W.mCnlr.sd[i] = sd(WES.prostate.cnlr_df$cnlr[which(WES.prostate.cnlr_df$ID == DMP.CMO.match_df$WES.ID[i])], na.rm = T)
  #' IMPACT data
  DMP.CMO.match_df$I.QC[i] = IMPACT.prostate.QC_df$QC[which(IMPACT.prostate.QC_df$ID == DMP.CMO.match_df$DMP_Sample_ID[i])]
  DMP.CMO.match_df$I.WGD[i] = IMPACT.prostate.QC_df$wgd[which(IMPACT.prostate.QC_df$ID == DMP.CMO.match_df$DMP_Sample_ID[i])]
  DMP.CMO.match_df$I.ploidy[i] = IMPACT.prostate.binary_df$ploidy[which(IMPACT.prostate.binary_df$sample.id == DMP.CMO.match_df$DMP_Sample_ID[i])]
  DMP.CMO.match_df$I.purity[i] = IMPACT.prostate.binary_df$purity[which(IMPACT.prostate.binary_df$sample.id == DMP.CMO.match_df$DMP_Sample_ID[i])]
  DMP.CMO.match_df$I.50_call[i] = IMPACT.prostate.binary_df$Y_call[which(IMPACT.prostate.binary_df$sample.id == DMP.CMO.match_df$DMP_Sample_ID[i])]
  DMP.CMO.match_df$I.Cnlr_segment[i] = median(IMPACT.prostate.processed_df$cnlr.median[which(IMPACT.prostate.processed_df$ID == DMP.CMO.match_df$DMP_Sample_ID[i] &
                                                                                            IMPACT.prostate.processed_df$chrom == 24)], na.rm = T)
  DMP.CMO.match_df$I.mCnlr[i] = median(IMPACT.prostate.cnlr_df$cnlr[which(IMPACT.prostate.cnlr_df$ID == DMP.CMO.match_df$DMP_Sample_ID[i])], na.rm = T)
  DMP.CMO.match_df$I.mCnlr.sd[i] = sd(IMPACT.prostate.cnlr_df$cnlr[which(IMPACT.prostate.cnlr_df$ID == DMP.CMO.match_df$DMP_Sample_ID[i])], na.rm = T)
  
}


#' number samples match:
sum(DMP.CMO.match_df$W.50_call == DMP.CMO.match_df$I.50_call)
sum(DMP.CMO.match_df$W.50_call != DMP.CMO.match_df$I.50_call)
DMP.CMO.match_df$mismatch = ifelse(DMP.CMO.match_df$W.50_call == DMP.CMO.match_df$I.50_call, 'match', 'mismatch')

# write.table(DMP.CMO.match_df, file = 'Data_out/WES_IMPACT.combinedCalls.txt', sep = '\t', row.names = F, quote = F)

##-----------------------------------------------------------------------------
## Either one (IMPACT or WES needs to be FacetsQC TRUE)
## ----------------------------------------------------------------------------
DMP.CMO.match.QC_TRUE = DMP.CMO.match_df[which(DMP.CMO.match_df$W.QC == TRUE | DMP.CMO.match_df$I.QC == TRUE), ]
# write.table(DMP.CMO.match.QC_TRUE, file = 'Data_out/WES_IMPACT.combinedCalls.QC.TRUE.txt', sep = '\t', row.names = F, quote = F)

## Visualization:
#' overall concordance between WES and IMPACT calls
purity.plot = ggplot(DMP.CMO.match.QC_TRUE, aes(x = W.purity, y = I.purity)) +
  geom_point() + 
  geom_abline(slope = coef(lm(DMP.CMO.match.QC_TRUE$W.purity ~ DMP.CMO.match.QC_TRUE$I.purity))[[2]], 
              intercept = coef(lm(DMP.CMO.match.QC_TRUE$W.purity ~ DMP.CMO.match.QC_TRUE$I.purity))[[1]]) +
  theme_Y(base_size = 14) +
  scale_x_continuous(expand = c(0.01,0.01),
                     limits = c(0, 1)) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0, 1)) +
  theme(aspect.ratio = 1) +
  labs(x = 'WES [purity]', y = 'IMPACT [purity]') +
  annotate('text',
           x = 0.05,
           y = 0.95,
           label = paste0('p[Pearson]=', round(cor.test(DMP.CMO.match.QC_TRUE$I.purity, DMP.CMO.match.QC_TRUE$W.purity)[[4]][[1]], 3)),
           hjust = 0, 
           vjust = 0,
           family = 'RobotoCondensed-Regular')
  
#' ploidy plot
ploidy.plot = ggplot(DMP.CMO.match.QC_TRUE, aes(x = W.ploidy, y = I.ploidy)) +
  geom_point() + 
  geom_abline(slope = coef(lm(DMP.CMO.match.QC_TRUE$W.ploidy ~ DMP.CMO.match.QC_TRUE$I.ploidy))[[2]], 
              intercept = coef(lm(DMP.CMO.match.QC_TRUE$W.ploidy ~ DMP.CMO.match.QC_TRUE$I.ploidy))[[1]]) +
  theme_Y(base_size = 14) +
  scale_x_continuous(expand = c(0.01,0.01),
                     limits = c(0, 6)) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0, 6)) +
  theme(aspect.ratio = 1) +
  labs(x = 'WES [ploidy]', y = 'IMPACT [ploidy]') +
  annotate('text',
           x = 0.3,
           y = 5.7,
           label = paste0('p[Pearson]=', round(cor.test(DMP.CMO.match.QC_TRUE$I.ploidy, DMP.CMO.match.QC_TRUE$W.ploidy)[[4]][[1]], 3)),
           hjust = 0, 
           vjust = 0,
           family = 'RobotoCondensed-Regular')


#' median CnLog Ratio along the Y-axis.
Y_continious_loss = ggplot(DMP.CMO.match.QC_TRUE, 
       aes(x = W.mCnlr, 
           y = I.mCnlr, 
           color = mismatch)) +
  geom_point() +
  geom_abline(slope = coef(lm(DMP.CMO.match.QC_TRUE$W.mCnlr ~ DMP.CMO.match.QC_TRUE$I.mCnlr))[[2]], 
              intercept = coef(lm(DMP.CMO.match.QC_TRUE$W.mCnlr ~ DMP.CMO.match.QC_TRUE$I.mCnlr))[[1]],
              linetype = 'dashed',
              size = 0.1,
              color = 'grey55') +
  geom_vline(xintercept = -.2,
             color = 'grey55',
             size = 0.1,
             linetype = 'dashed') +
  geom_hline(yintercept = -.2, 
             size = 0.1, 
             linetype = 'dashed',
             color = 'grey55') +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(-4.5, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(-4.5, 1.5)) +
  scale_color_manual(values = c('mismatch' = '#F04E33',
                                'match' = '#006CA2'),
                     name = 'Method',
                     guide = guide_legend(direction = 'horizontal',
                                          title.theme = element_text(size = 12, face = 'bold', family = 'RobotoCondensed-Regular'),
                                          label.theme = element_text(size = 12, family = 'RobotoCondensed-Regular'))) +
  
  annotate('text',
           x = -4.2,
           y = 1.2,
           label = paste0('p[Pearson]=', round(cor.test(DMP.CMO.match.QC_TRUE$I.mCnlr, DMP.CMO.match.QC_TRUE$W.mCnlr)[[4]][[1]], 3)),
           hjust = 0, 
           vjust = 0,
           family = 'RobotoCondensed-Regular',
           size = 4) +

  
  theme_Y(base_size = 14) +
  theme(aspect.ratio = 1,
        legend.position = 'top') +
  
  labs(x = 'WES [median Copy Number Log Ratio]',
       y = 'IMPACT [median Copy Number Log Ratio]')
  

panelA = ploidy.plot + purity.plot + Y_continious_loss
ggsave_golden(filename = 'Figures/Y.loss.continious.pdf', plot = panelA, width = 16)


#' Visualize the mis-matched samples --
Y.call_mismatch = ggplot(DMP.CMO.match.QC_TRUE, 
                           aes(x = W.mCnlr, 
                               y = I.mCnlr, 
                               color = mismatch)) +
  geom_point(aes(alpha = mismatch)) +
  geom_vline(xintercept = -.2,
             color = 'grey55',
             size = 0.1,
             linetype = 'dashed') +
  geom_hline(yintercept = -.2, 
             size = 0.1, 
             linetype = 'dashed',
             color = 'grey55') +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(-4.5, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(-4.5, 1.5)) +
  
  scale_alpha_manual(values = c('match' = 0.2,
                                'mismatch' = 1),
                     guide = FALSE) + 
  
  scale_color_manual(values = c('mismatch' = '#F04E33',
                                'match' = 'grey65'),
                     name = 'Method',
                     guide = guide_legend(direction = 'horizontal',
                                          title.theme = element_text(size = 12, face = 'bold', family = 'RobotoCondensed-Regular'),
                                          label.theme = element_text(size = 12, family = 'RobotoCondensed-Regular'))) +
  
  theme_Y(base_size = 14) +
  theme(aspect.ratio = 1,
        legend.position = 'top') +
  
  labs(x = 'WES [median Copy Number Log Ratio]',
       y = 'IMPACT [median Copy Number Log Ratio]')

ggsave_golden(filename = 'Figures/Y.call_mismatch.pdf', plot = Y.call_mismatch, width = 14)

## out
