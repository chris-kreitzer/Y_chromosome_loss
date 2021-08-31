## Method validation of WES and IMPACT sequencing results
## Show that both methods are suitable in calling Y-chromosome loss
## We are using the median CnLR as surrogate for the binary Y-chromosome call;
## Just show whether there is a significant association or not; WES provides us with a higher resolution;
## 
## Furthermore, we can expand the correlation analysis to include ploidy and purity estimates as well;
## 
## Start (revision): 08/30/2021
## chris-kreitzer
##

set.seed(12)
setwd('~/Documents/GitHub/Y_chromosome_loss/PanCancer')

## Input
WES_cnLR = read.csv('Data_out/WES/WES_cnLR_out.txt', sep = '\t')
IMPACT_cnLR = read.csv('Data_out/IMPACT/Cnlr_out.txt', sep = '\t')
WES.IMPACT.ids_df = read.csv('Data_in/WES_cBio_ID.match.txt', sep = '\t')


## Data modifications:
WES_cnLR$sample_id = substr(x = WES_cnLR$ID, start = 3, stop = 17)
WES_cnLR$ID = NULL


## Processing; 
#' calculate median LRR for Y-chromosome; WES
WES_cnLR$mLRR = NA
for(i in unique(WES_cnLR$sample_id)){
  if(length(WES_cnLR$sample_id[which(WES_cnLR$sample_id == i)]) >= 10){
    median.LRR = median(WES_cnLR$cnlr[which(WES_cnLR$sample_id == i)])
    WES_cnLR$mLRR[which(WES_cnLR$sample_id == i)] = median.LRR
  } else {
    median.LRR = 0
    WES_cnLR$mLRR[which(WES_cnLR$sample_id == i)] = median.LRR
  }
}

#' IMPACT
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

#' merge the two data frames;
Seq_merged = WES.IMPACT.ids_df[, c('DMP_Sample_ID', 'CMO_Sample_ID', 'Sex')]
Seq_merged$mLRR_IMPACT = NA
Seq_merged$mLRR_WES = NA

for(i in 1:nrow(Seq_merged)){
  if(Seq_merged$DMP_Sample_ID[i] %in% IMPACT_cnLR$ID){
    Seq_merged$mLRR_IMPACT[i] = IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == Seq_merged$DMP_Sample_ID[i])]
  } else {
    Seq_merged$mLRR_IMPACT[i] = 0
  }
}

for(i in 1:nrow(Seq_merged)){
  if(Seq_merged$CMO_Sample_ID[i] %in% WES_cnLR$sample_id){
    Seq_merged$mLRR_WES[i] = WES_cnLR$mLRR[which(WES_cnLR$sample_id == Seq_merged$CMO_Sample_ID[i])]
  } else {
    Seq_merged$mLRR_WES[i] = 0
  }
}

#' subset to only those samples where we have both, IMPACT and WES data
Seq_merged = Seq_merged[which(Seq_merged$mLRR_IMPACT != 0 & Seq_merged$mLRR_WES != 0), ]


## Visualization:
Y_continious_loss = ggplot(Seq_merged, 
                           aes(x = mLRR_WES, 
                               y = mLRR_IMPACT)) +
  geom_point() +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed',
              size = 0.1,
              color = 'grey55') +
  geom_vline(xintercept = 0,
             color = 'grey55',
             size = 0.1,
             linetype = 'dashed') +
  geom_hline(yintercept = 0, 
             size = 0.1, 
             linetype = 'dashed',
             color = 'grey55') +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(-5, 2.5)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(-5, 2.5)) +
  annotate('text',
           x = -4.5,
           y = 2.2,
           label = paste0('r = ', round(cor.test(Seq_merged$mLRR_IMPACT, Seq_merged$mLRR_WES)[[4]][[1]], 3), 
                          '\np < 2.2e-16'),
           hjust = 0, 
           vjust = 0,
           family = 'RobotoCondensed-Regular',
           size = 4) +
  
  theme_Y(base_size = 14) +
  theme(aspect.ratio = 1,
        legend.position = 'top') +
  
  labs(x = 'WES [median Copy Number Log Ratio]',
       y = 'IMPACT [median Copy Number Log Ratio]',
       title = 'CnLR correlation of WES and IMPACT; n = 950')


Y_continious_loss
ggsave_golden(filename = 'Figures/Method_validation.pdf', plot = Y_continious_loss, width = 16)


#' out