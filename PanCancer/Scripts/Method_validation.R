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

## chris-kreitzer

clean()
gc()
.rs.restartR()
set.seed(100)
setwd('~/Documents/GitHub/Y_chromosome_loss/PanCancer')
source('Scripts/UtilityFunctions.R')


## Input
WES_cnLR = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/WES_cnLR_out.txt', sep = '\t')
IMPACT_cnLR = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/IMPACT_Cnlr_out.txt', sep = '\t')
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
WES_cnLR$ID = substr(x = WES_cnLR$ID, start = 3, stop = 17)

##----------------+
## median LRR for IMPACT
##----------------+
IMPACT_cnLR$mLRR = NA
for(i in unique(IMPACT_cnLR$ID)){
  if(length(IMPACT_cnLR$ID[which(IMPACT_cnLR$ID == i)]) >= 10){
    median.LRR = median(IMPACT_cnLR$cnlr[which(IMPACT_cnLR$ID == i)])
    IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == i)] = median.LRR
  } else {
    median.LRR = 0
    IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == i)] = median.LRR
  }
}


##----------------+
## merge datasets
##----------------+
Seq_merged = WES.IMPACT.ids_df[, c('DMP_Sample_ID', 'CMO_Sample_ID', 'Sex')]
Seq_merged$mLRR_IMPACT = NA
Seq_merged$mLRR_WES = NA

for(i in 1:nrow(Seq_merged)){
  print(i)
  if(Seq_merged$DMP_Sample_ID[i] %in% IMPACT_cnLR$ID){
    Seq_merged$mLRR_IMPACT[i] = IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == Seq_merged$DMP_Sample_ID[i])]
  } else {
    Seq_merged$mLRR_IMPACT[i] = NA
  }
}

for(i in 1:nrow(Seq_merged)){
  print(i)
  if(Seq_merged$CMO_Sample_ID[i] %in% WES_cnLR$ID){
    Seq_merged$mLRR_WES[i] = WES_cnLR$mLRR[which(WES_cnLR$ID == Seq_merged$CMO_Sample_ID[i])]
  } else {
    Seq_merged$mLRR_WES[i] = NA
  }
}

Seq_merged = Seq_merged[!with(Seq_merged, is.na(mLRR_IMPACT) & is.na(mLRR_WES)), ]
write.table(Seq_merged, file = '~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/mLRR_IMPACT_WES_n950.txt', sep = '\t', quote = F, row.names = F)

##----------------+
## Visualization
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

  # geom_vline(xintercept = 0,
  #            color = 'grey55',
  #            linewidth = 0.1,
  #            linetype = 'dashed') +
  # geom_hline(yintercept = 0, 
  #            linewidth = 0.1, 
  #            linetype = 'dashed',
  #            color = 'grey55') +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(-5, 2.5)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(-5, 2.5)) +
  annotate('text',
           x = -3.5,
           y = 1.5,
           label = paste0('r = ', round(cor.test(Seq_merged$mLRR_IMPACT, Seq_merged$mLRR_WES)[[4]][[1]], 3), 
                          '\np < 2.2e-16'),
           hjust = 0, 
           vjust = 0,
           family = 'ArialMT',
           size = 6) +
  
  theme(aspect.ratio = 1,
        legend.position = 'top',
        panel.border = element_rect(fill = NA, linewidth = 1.5),
        axis.ticks = element_blank(),
        axis.text = element_text(colour = 'black', size = 12),
        panel.background = element_rect(fill = 'white')) +
  
  labs(x = 'Exome re-sequencing [Y-CnLR]',
       y = 'MSK-IMPACT [Y-CnLR]',
       title = 'Y-chromosome CnLR; n = 950')


Y_continious_loss


ggsave_golden(filename = 'Figures/Method_validation.pdf', plot = Y_continious_loss, width = 10)



##----------------+
## Base R solution; plot
##----------------+
dev.off()
par(mfrow = c(1, 3))
plot(Seq_merged$mLRR_IMPACT, Seq_merged$mLRR_WES,
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = '',
     pch = 19)
axis(side = 1, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3), 
     cex = 1.3,
     line = 0.3, 
     lwd = 1.5, 
     lwd.ticks = 1.5)
axis(side = 2, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3), 
     cex = 1.3, 
     las = 2, 
     line = 0.3, 
     lwd = 1.5,
     lwd.ticks = 1.5)
abline(a = 0, b = 1, lty = 'dashed', col = 'grey35')
box(lwd = 2)

mtext(text = 'Whole Exome Sequencing [mLRR]', side = 1, line = 2.3)
mtext(text = 'MSK-IMPACT Sequencing [mLRR]', side = 2, line = 2.3)
text(x = -4.7, y = 1.7, labels = 'r = 97.5\n p < 2.2e-16')


##-----------------
## SAMPLE OUTLIER:
##-----------------
# P-0019114-T01-IM6
# P-0037069-T01-IM6
# P-0047344-T01-IM6




#' Purity and Ploidy estimates of both methods (supplement for validation)
cohortData = readRDS(file = 'Data_out/cohort_data.rds')
IMPACT = cohortData$IMPACT.cohort
IMPACT$counts_file = NULL
WES = cohortData$WES.cohort
WES$Facet_Countfile = NULL
WES$Facet_Path = NULL
WES$CMO_Sample_ID = NULL

ID_match = cohortData$ID_matches
ID_match$Facet_Countfile = NULL
ID_match$Facet_Path = NULL

#' create matching data.frame for IMPACT vs WES
match_out = data.frame()
for(i in 1:nrow(WES)){
  print(WES$User_Sample_ID[i])
  WES.sample = WES$User_Sample_ID[i]
  WES.purity = WES$purity[i]
  WES.ploidy = WES$ploidy[i]
  IMPACT.sample = ID_match$DMP_Sample_ID[which(ID_match$User_Sample_ID == WES.sample)]
  print(IMPACT.sample)
  IMPACT.purity = IMPACT$purity[which(IMPACT$SAMPLE_ID == IMPACT.sample)]
  IMPACT.ploidy = IMPACT$ploidy[which(IMPACT$SAMPLE_ID == IMPACT.sample)]
  
  out = data.frame(WES = WES.sample,
                   WES.purity = WES.purity,
                   WES.ploidy = WES.ploidy,
                   IMPACT = IMPACT.sample,
                   IMPACT.purity = IMPACT.purity,
                   IMPACT.ploidy = IMPACT.ploidy)

  match_out = rbind(match_out, out)
  
}

#' Visualization
match_out = match_out[which(match_out$WES.purity != 0 & match_out$IMPACT.purity != 0), ]

#' purity agreement plot
plot(match_out$WES.purity, match_out$IMPACT.purity,
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = '',
     pch = 19,
     xlim = c(0, 1),
     ylim = c(0, 1))
axis(side = 1, at = seq(0, 1, 0.2), 
     cex = 1.3,
     line = 0.3, 
     lwd = 1.5, 
     lwd.ticks = 1.5)
axis(side = 2, at = seq(0, 1, 0.2), 
     cex = 1.3, 
     las = 2, 
     line = 0.3, 
     lwd = 1.5,
     lwd.ticks = 1.5)
abline(a = 0, b = 1, lty = 'dashed', col = 'grey35')
box(lwd = 2)

mtext(text = 'Purity Whole Exome Sequencing [%]', side = 1, line = 2.3)
mtext(text = 'Purity MSK-IMPACT Sequencing [%]', side = 2, line = 2.3)
text(x = 0.15, y = 0.85, labels = 'r = 90.0\n p < 2.2e-16')
cor.test(match_out$WES.purity, match_out$IMPACT.purity)

## ploidy
plot(match_out$WES.ploidy, match_out$IMPACT.ploidy,
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = '',
     pch = 19,
     xlim = c(0, 6),
     ylim = c(0,6))
axis(side = 1, at = seq(1, 6, 1), 
     cex = 1.3,
     line = 0.3, 
     lwd = 1.5, 
     lwd.ticks = 1.5)
axis(side = 2, 
     at = seq(1, 6, 1),
     labels = seq(1, 6, 1),
     cex = 1.3, 
     las = 2, 
     line = 0.3, 
     lwd = 1.5,
     lwd.ticks = 1.5)
abline(a = 0, b = 1, lty = 'dashed', col = 'grey35')
box(lwd = 2)

mtext(text = 'Ploidy Whole Exome Sequencing', side = 1, line = 2.3)
mtext(text = 'Ploidy MSK-IMPACT Sequencing', side = 2, line = 2.3)
text(x = 1.1, y = 5.5, labels = 'r = 62.5\n p < 2.2e-16')
cor.test(match_out$WES.ploidy, match_out$IMPACT.ploidy)


purity.plot = ggplot(match_out, aes(x = WES.purity, y = IMPACT.purity)) +
  geom_jitter(size = 0.6) +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed',
              size = 0.2,
              color = 'grey35') +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  theme_Y +
  theme(aspect.ratio = 1) +
  labs(x = 'WES recapture', y = 'IMPACT targeted panel', 
       title = paste0('Purity; n = 946\n r = ', round(cor.test(match_out$WES.purity, match_out$IMPACT.purity)[[4]][1], 3)))


ggsave_golden(filename = 'Figures/PurityValidation.pdf', plot = purity.plot, width = 8)


#' ploidy agreement plot
ploidy.plot = ggplot(match_out, aes(x = WES.ploidy, y = IMPACT.ploidy)) +
  geom_jitter() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(1, 5)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(1, 5)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  labs(x = 'WES recapture', y = 'IMPACT targeted panel', 
       title = paste0('Ploidy agreement plot\nn = 946, r = ', round(cor.test(match_out$WES.ploidy, match_out$IMPACT.ploidy)[[4]][1], 3)))

Supplement1 = ploidy.plot + purity.plot
ggsave_golden(Supplement1, filename = 'Figures/PurityPloidy_Validation.pdf', width = 12)


#' out