##----------------+
## Confirm IMPACT chromosome Y 
## calls with WES sequenced
## samples; using CnLR as surrogate
## for binary Y-chromosome calls;
##----------------+ 
## 
## Start (revision): 08/30/2021
## revision: 07/26/2022
## revision: 08/21/2022
## revision: 12/16/2022
## revision: 01/30/2023
##
## chris-kreitzer

clean()
gc()
.rs.restartR()
set.seed(100)
setwd('~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')


## Input
WES_cnLR = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/WES/WES_Cnlr_Y.txt', sep = '\t')
IMPACT_cnLR = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/04_Loss/010523/Cnlr_Y.txt', sep = '\t')
WES.IMPACT.ids_df = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/WES_cBio_ID.match.txt', sep = '\t')



##----------------+
## median LRR for WES 
## chromosome Y
##-----------------+
WES_cnLR$mLRR = NA
for(i in unique(WES_cnLR$ID)){
  print(i)
  if(length(WES_cnLR$ID[which(WES_cnLR$ID == i)]) >= 10){
    WES_cnLR$mLRR[which(WES_cnLR$ID == i)] = median(WES_cnLR$cnlr[which(WES_cnLR$ID == i)], na.rm = T)
  } else {
    WES_cnLR$mLRR[which(WES_cnLR$ID == i)] = 0
  }
}


##----------------+
## median LRR for IMPACT
##----------------+
IMPACT_cnLR$mLRR = NA
for(i in unique(IMPACT_cnLR$ID)){
  print(i)
  if(length(IMPACT_cnLR$ID[which(IMPACT_cnLR$ID == i)]) >= 10){
    median.LRR = median(IMPACT_cnLR$cnlr[which(IMPACT_cnLR$ID == i)])
    IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == i)] = median.LRR
  } else {
    median.LRR = 0
    IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == i)] = median.LRR
  }
}


##----------------+
## merge the matched
## validation cohort
## - add further info
## - purity, QC
##----------------+
WES_cnLR$ID = substr(x = WES_cnLR$ID, start = 3, stop = 17)
Seq_merged = WES.IMPACT.ids_df[which(WES.IMPACT.ids_df$CMO_Sample_ID %in% unique(WES_cnLR$ID)), c('DMP_Sample_ID', 'CMO_Sample_ID', 'Sex')]

Seq_merged$mLRR_IMPACT = NA
Seq_merged$mLRR_WES = NA

for(i in 1:nrow(Seq_merged)){
  print(i)
  if(Seq_merged$DMP_Sample_ID[i] %in% IMPACT_cnLR$ID){
    Seq_merged$mLRR_IMPACT[i] = unique(IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == Seq_merged$DMP_Sample_ID[i])])
  } else {
    Seq_merged$mLRR_IMPACT[i] = NA
  }
}

for(i in 1:nrow(Seq_merged)){
  print(i)
  if(Seq_merged$CMO_Sample_ID[i] %in% WES_cnLR$ID){
    Seq_merged$mLRR_WES[i] = unique(WES_cnLR$mLRR[which(WES_cnLR$ID == Seq_merged$CMO_Sample_ID[i])])
  } else {
    Seq_merged$mLRR_WES[i] = NA
  }
}


##-------
## add. annotation
##-------
IMPACT_cna = read.csv('Data/04_Loss/010523/CopyNumberStates.txt', sep = '\t')
WES_cna = read.csv('Data/02_Method_Validation/WES/WES_CopyNumberStates.txt', sep = '\t')
WES_cna$id = substr(x = WES_cna$id, start = 3, stop = 17)

Seq_merged = merge(Seq_merged, IMPACT_cna[, c('id', 'QC', 'purity', 'ploidy')], by.x = 'DMP_Sample_ID', by.y = 'id', all.x = T)
Seq_merged = Seq_merged[!duplicated(Seq_merged), ]
colnames(Seq_merged)[6:8] = c('IMPACT_QC', 'IMPACT_purity', 'IMPACT_ploidy')

Seq_merged = merge(Seq_merged, WES_cna[,c('id', 'QC', 'purity', 'ploidy')], by.x = 'CMO_Sample_ID', by.y = 'id', all.x = T)
Seq_merged = Seq_merged[!duplicated(Seq_merged), ]
colnames(Seq_merged)[9:11] = c('WES_QC', 'WES_purity', 'WES_ploidy')


Seq_merged = Seq_merged[!with(Seq_merged, is.na(mLRR_IMPACT) | is.na(mLRR_WES)), ]
Seq_merged = Seq_merged[which(Seq_merged$IMPACT_QC == TRUE), ]
Seq_merged = Seq_merged[which(Seq_merged$WES_QC == TRUE), ]
write.table(Seq_merged, file = '~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/ValidationCohort_out_n769.txt', sep = '\t', quote = F, row.names = F)




##----------------+
## Visualization
## - first mLRR of both cohorts
## - then purity
##----------------+
Y_continious_loss = ggplot(Seq_merged, 
                           aes(x = mLRR_WES, 
                               y = mLRR_IMPACT)) +
  geom_point() +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed',
              linewidth = 0.35,
              color = 'grey35') +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(-5, 2.5)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(-5, 2.5)) +
  theme_std(base_size = 14, base_line_size = 0.5) +
  annotate('text',
           x = -3.8,
           y = 1.5,
           label = paste0('r = ', round(cor.test(Seq_merged$mLRR_IMPACT, Seq_merged$mLRR_WES)[[4]][[1]], 2), 
                          '\np < 2.2e-16'),
           hjust = 0, 
           vjust = 0,
           family = 'ArialMT',
           size = 6) +
  
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, linewidth = 2),
        panel.background = element_rect(fill = 'white')) +
  
  labs(x = 'WES-recapture [median Cnlr]',
       y = 'MSK-IMPACT [median Cnlr]',
       title = paste0('Y-chromosome; n = ', length(unique(Seq_merged$CMO_Sample_ID))))


Y_continious_loss


##----------------+
## Purity comparison
##----------------+
Purity_plot = ggplot(Seq_merged, 
                     aes(x = WES_purity, 
                         y = IMPACT_purity)) +
  geom_point() +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed',
              linewidth = 0.35,
              color = 'grey35') +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 1)) +
  theme_std(base_size = 14, base_line_size = 0.5) +
  annotate('text',
           x = 0.15,
           y = 0.85,
           label = paste0('r = ', round(cor.test(Seq_merged$IMPACT_purity, Seq_merged$WES_purity)[[4]][[1]], 2), 
                          '\np < 2.2e-16'),
           hjust = 0, 
           vjust = 0,
           family = 'ArialMT',
           size = 6) +
  
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, linewidth = 2),
        panel.background = element_rect(fill = 'white')) +
  
  labs(x = 'purity WES-recapture',
       y = 'purity MSK-IMPACT',
       title = paste0('Y-chromosome; n = ', length(unique(Seq_merged$CMO_Sample_ID))))



ggsave_golden(filename = 'Figures_original/Method_validation_mCnLR.pdf', plot = Y_continious_loss, width = 8)
ggsave_golden(filename = 'Figures_original/Method_Purity_validation.pdf', plot = Purity_plot, width = 8)


method_validation = Y_continious_loss + Purity_plot
ggsave_golden(filename = '~/Documents/MSKCC/10_MasterThesis/Figures_original/Method_validation_MSK.pdf', plot = method_validation, width = 12)



##-----------------
## SAMPLE OUTLIER:
##-----------------
# P-0019114-T01-IM6
# P-0037069-T01-IM6
# P-0047344-T01-IM6



#' out