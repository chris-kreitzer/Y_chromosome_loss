## PanCancer Analysis of Y-chromosome loss
## Investigate mis-matched samples (WES vs IMPACT)
## 
## Gaussian mixture model; 
## 
## Start analysis: 04/14/2021
## Revision: 04/19/2021


set.seed(3491)
rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')
source('Y_chromosome_loss/Plotting_theme.R')


## Libraries
library(mixtools)
library(ggplot2)


## Input
WES.IMPACT.combined.calls = read.csv('Data_out/WES_IMPACT.combinedCalls.txt', sep = '\t')
Combined.QC.TRUE = read.csv('Data_out/WES_IMPACT.combinedCalls.QC.TRUE.txt', sep = '\t')
IMPACT.IGV = read.csv('Data_out/IMPACT.prostate.IGV.16.4.seg', sep = '\t')
IMPACT.cnlr = read.csv('Data_out/Cnlr.IMPACT.prostate.out.16.4.txt', sep = '\t')
WES.cnlr = read.csv('Data_out/Cnlr.WES.prostate.out.15.4.txt', sep = '\t')


## Functions
plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}


## Processing:
#' overall overview
Matched_QC.TRUE = Combined.QC.TRUE[which(Combined.QC.TRUE$W.50_call == Combined.QC.TRUE$I.50_call), ]
Mismatched_QC.TRUE = Combined.QC.TRUE[which(Combined.QC.TRUE$W.50_call != Combined.QC.TRUE$I.50_call), ]

#' IGV from IMPACT QC true samples
IMPACT.QC.TRUE = Combined.QC.TRUE$DMP_Sample_ID[which(Combined.QC.TRUE$I.QC == TRUE)]
IMPACT.IGV$ID = IMPACT.IGV$ID = paste0(substr(IMPACT.IGV$ID, start = 8, stop = 16), '-T', 
                                                  substr(IMPACT.IGV$ID, start = 1, stop = 6))
IMPACT.IGV.QC.TRUE = IMPACT.IGV[which(IMPACT.IGV$ID %in% IMPACT.QC.TRUE), ]
# write.table(IMPACT.IGV.QC.TRUE, file = 'Data_out/IMPACT.IGV_QC.TRUE.seg', sep = '\t', row.names = F, quote = F)

IMPACT.cnlr$ID = IMPACT.cnlr$ID = paste0(substr(IMPACT.cnlr$ID, start = 8, stop = 16), '-T', 
                                                     substr(IMPACT.cnlr$ID, start = 1, stop = 6))
WES.cnlr$ID = substr(WES.cnlr$ID, start = 3, stop = 19)


## Analysis:
##-----------------------------------------------------------------------------
## Analyse one mis-matched sample: less extreme
##-----------------------------------------------------------------------------
a = IMPACT.cnlr[which(IMPACT.cnlr$ID == 'P-0020482-T02-IM6'), ]
b = WES.cnlr[which(WES.cnlr$ID == 's_C_001131_P001_d'), ]

cnlr = a$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
IMPACT.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio',
       title = 'IMPACT Cn-LogR Distribution across the Y-chromosome',
       subtitle = '94% of markers belong to \'red\' distribution with mu = 0.10') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))

ggsave_golden(filename = 'Figures/IMPACT.GMM.example.pdf', plot = IMPACT.mix.example, width = 12)

#' corresponding WES example:
cnlr = b$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio',
       title = 'WES Cn-LogR Distribution across the Y-chromosome',
       subtitle = '88% of markers belong to \'red\' distribution with mu = -0.18') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))

ggsave_golden(filename = 'Figures/WES.GMM.example.pdf', plot = WES.mix.example, width = 12)


##-----------------------------------------------------------------------------
## Analyse a second mis-matched example: more extreme
## ----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Analyse one mis-matched sample:
##-----------------------------------------------------------------------------
a = IMPACT.cnlr[which(IMPACT.cnlr$ID == 'P-0033767-T01-IM6'), ]
b = WES.cnlr[which(WES.cnlr$ID == 's_C_UF36VT_M001_d'), ]

cnlr = a$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
IMPACT.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio',
       title = 'IMPACT Cn-LogR Distribution across the Y-chromosome',
       subtitle = '51.7% of markers belong to \'blue\' distribution with mu = -0.16') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))

ggsave_golden(filename = 'Figures/IMPACT.GMM.example2.extreme.pdf', plot = IMPACT.mix.example, width = 12)

#' corresponding WES example:
cnlr = b$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio',
       title = 'WES Cn-LogR Distribution across the Y-chromosome',
       subtitle = '51.1% of markers belong to \'red\' distribution with mu = -3.25') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))

ggsave_golden(filename = 'Figures/WES.GMM.example2.extreme.pdf', plot = WES.mix.example, width = 12)


#------------------------------------------------------------------------------
# Writing the analysis pipeline:
# define ruleset
GMM.recover = data.frame()
for(i in 1:nrow(Mismatched_QC.TRUE)){
  try({
    WES.QC = Mismatched_QC.TRUE$W.QC[i]
    IMPACT.QC = Mismatched_QC.TRUE$I.QC[i]
    
    #' select with which sample to go in GMM
    if(WES.QC & IMPACT.QC){
      selected.sample = Mismatched_QC.TRUE$WES.ID[i]
      GMM.input = WES.cnlr[which(WES.cnlr$ID == selected.sample), 'cnlr']
    } else if(WES.QC){
      selected.sample = Mismatched_QC.TRUE$WES.ID[i]
      GMM.input = WES.cnlr[which(WES.cnlr$ID == selected.sample), 'cnlr']
    } else if (!WES.QC & IMPACT.QC) {
      selected.sample = Mismatched_QC.TRUE$DMP_Sample_ID[i]
      GMM.input = IMPACT.cnlr[which(IMPACT.cnlr$ID == selected.sample), 'cnlr']
    }
    
    #' GMM on selected sample
    GMM.processed = mixtools::normalmixEM(x = GMM.input, k = 2)
    GMM.lambda = GMM.processed$lambda
    GMM.mu = GMM.processed$mu
    
    if(any(GMM.lambda >= 0.75)){
      called.mean = GMM.mu[which.max(GMM.lambda)]
      lambda.max = GMM.lambda[which.max(GMM.lambda)]
      Y_call = ifelse(called.mean <= -0.20, 'loss', 'intact')
    } else {
      called.mean = GMM.mu[which.max(GMM.lambda)]
      lambda.max = GMM.lambda[which.max(GMM.lambda)]
      Y_call = 'ambigiuous'
    }
    
    processed.sample = data.frame(ID = selected.sample,
                                  lambda.max = lambda.max,
                                  mean.max = called.mean,
                                  Y_call = Y_call)
    
    GMM.recover = rbind(GMM.recover, processed.sample)
    
  })
  
}


#' full combined output
Subset1 = Matched_QC.TRUE[,c('DMP_Sample_ID', 'I.50_call')]
colnames(Subset1)[2] = 'Y_chromosome'
Subset1$Y_chromosome = ifelse(Subset1$Y_chromosome == 'intact_Y_chrom', 'intact', 'loss')
Subset2 = GMM.recover[,c('ID', 'Y_call')]

Subset2.wes = Subset2[!grepl(pattern = '^P-.*', x = Subset2$ID), ]
Subset2.wes = merge(Subset2.wes, Combined.QC.TRUE[,c('DMP_Sample_ID', 'WES.ID')],
                    by.x = 'ID', by.y = 'WES.ID', all.x = T)
Subset2_a = Subset2.wes[,c('DMP_Sample_ID', 'Y_call')]
colnames(Subset2_a)[2] = 'Y_chromosome'
Subset2_b = Subset2[grepl(pattern = '^P-.*', x = Subset2$ID), c('ID', 'Y_call')]
colnames(Subset2_b)[1] = 'DMP_Sample_ID'
colnames(Subset2_b)[2] = 'Y_chromosome'

calls.possible = unique(c(Subset1$DMP_Sample_ID, Subset2_a$DMP_Sample_ID, Subset2_b$DMP_Sample_ID))
QC.FALSE = setdiff(WES.IMPACT.combined.calls$DMP_Sample_ID, calls.possible)
QC.FALSE = data.frame(DMP_Sample_ID = QC.FALSE,
                      Y_chromosome = 'N/A')

Y_calls_prostate = rbind(Subset1, Subset2_a, Subset2_b, QC.FALSE)
write.table(Y_calls_prostate, file = 'Data_out/Y_loss_calls_prostate.19.04.txt', sep = '\t', row.names = F)

## out