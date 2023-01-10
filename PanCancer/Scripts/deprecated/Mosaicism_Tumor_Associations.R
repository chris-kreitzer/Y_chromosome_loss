## Association study: Is there a connection between Mosaic loss and Y-chromosome loss in 
## solid cancer tissue; association study
## 
## 11/08/2021
## chris-kreitzer

rm(list = ls())
.rs.restartR()
setwd('~/Documents/GitHub/Y_chromosome_loss/PanCancer/')
source('Scripts/UtilityFunctions.R')
set.seed(99)
library(pctGCdata)


## Data
cohort = readRDS('Data_out/cohort_data.rds')
IMPACT_tumor = cohort$IMPACT.cohort
IMPACT_mosaic = cohort$IMPACT.mosaic

IMPACT_merged = merge(IMPACT_tumor[, c('SAMPLE_ID', 'Y_call', 'ploidy', 'purity')],
                      IMPACT_mosaic[which(IMPACT_mosaic$target == 'Y'), c('sample', 'corrected.CN', 'Y_mosaicism')], 
                      by.x = 'SAMPLE_ID', by.y = 'sample', all.x = T)

#' cross-tabulation
xtabs(~IMPACT_merged$Y_call + IMPACT_merged$Y_mosaicism)
chisq.test(xtabs(~IMPACT_merged$Y_call + IMPACT_merged$Y_mosaicism))

#' there is an association between Germline mosaic and Y-chromosome loss
#' Hypothesis: There is a technical association between Mosaic and loss in solid tumor tissue
#' the reason for this is simple: If there is mosaic - there is clearly less Y-chromosome DNA
#' content in the normal; for the tumor assessment we are basically using a T/N ratio.
#' If the N-DNA content is low -- the ratios are skewed anyway.
IMPACT_tumorCN = read.csv('Data_out/IMPACT/IMPACT_copynumber_out.txt', sep = '\t')
IMPACT_tumorCN = IMPACT_tumorCN[which(IMPACT_tumorCN$chrom == 24), ]

#' with tumor CN
IMPACT_merged = merge(IMPACT_tumorCN[, c('cnlr.median', 'tcn.em', 'lcn.em', 'ID')], IMPACT_merged,
                      by.x = 'ID', by.y = 'SAMPLE_ID', all.x = T)


#' CLASS 1: samples with mosaic Y-loss and an intact Y-chromosome in the tumor

#' Mosaic loss and Y-chromosome loss
IMPACT_merged$populationEstimate = NA

for(i in unique(IMPACT_merged$ID)){
  da = IMPACT_merged[which(IMPACT_merged$ID == i), ]
  if(nrow(da) > 1){
    IMPACT_merged$populationEstimate[which(IMPACT_merged$ID == i)] = mean(da$cnlr.median[which(da$ID == i)])
  } else {
    IMPACT_merged$populationEstimate[which(IMPACT_merged$ID == i)] = da$cnlr.median
  }
}

IMPACT_merged = IMPACT_merged[!is.na(IMPACT_merged$corrected.CN), ]

#' visualization:
class1.plot = ggplot(IMPACT_merged, aes(x = corrected.CN, y = populationEstimate)) +
  geom_jitter(size = 0.3) +
  geom_smooth(method = 'lm', se = F) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(-2, 2)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(0.5, 1.5),
                     breaks = c(0.5, 1, 1.5)) +
  theme_Y +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black', size = 1.4)) +
  labs(x = 'Germline Copy-Number', y = 'Tumor Copy-Number')

class1.plot
ggsave_golden(filename = 'Figures/GermlineTumorIMPACT_correlation.pdf', plot = class1.plot, width = 8)


summary(lm(IMPACT_merged$populationEstimate ~ IMPACT_merged$corrected.CN))

IMPACT_merged$ID[which(IMPACT_merged$corrected.CN > 1.5)]


#' look into P-0010668-T01-IM5 (high Germline DNA content; low Tumor DNA content (Y))
Germline_raw = read.csv('Mosaicism/IMPACT_coverage_bins.txt', sep = '\t')
Germline_raw = Germline_raw[which(Germline_raw$sample == 'countsMerged____P-0010668-T01-IM5_P-0010668-N01-IM5.dat.gz'), ]
cbio(cohort$IMPACT.cohort$counts_file[grepl(pattern = 'P-0010668-T01-IM5_P-0010668-N01-IM5.dat.gz', x = cohort$IMPACT.cohort$counts_file)])
counts = facetsY::readSnpMatrix('Data_out/tmp/countsMerged____P-0010668-T01-IM5_P-0010668-N01-IM5.dat.gz')
counts = facetsY::preProcSample(rcmat = counts, ndepth = 30, snp.nbhd = 0)
data.out = counts$jointseg[which(counts$jointseg$gcpct >= 0.4 & counts$jointseg$gcpct <= 0.6),]

#' out IMPACT


















#' samples where we have info on Y chromosome mosaicism
mosaic = cohort$WES.cohort[!is.na(cohort$WES.cohort$Y_mosaic), ]
mosaic$AGE_AT_SEQ_REPORTED_YEARS = as.numeric(as.character(mosaic$AGE_AT_SEQ_REPORTED_YEARS))
mosaic = mosaic[!is.na(mosaic$AGE_AT_SEQ_REPORTED_YEARS), ]

Loss_age = data.frame()
for(i in seq(5, 90, 5)){
  da = mosaic[between(mosaic$AGE_AT_SEQ_REPORTED_YEARS, i-4, i), ]
  
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

age_loss = ggplot(Loss_age[!Loss_age$group %in% c(5, 10, 15, 20), ], aes(x = group, y = ratio)) + 
  geom_hline(yintercept = seq(0, 50, 25), linetype = 'dashed', size = 0.2, color = 'grey35') +
  geom_bar(stat = 'identity') +
  stat_smooth(method = 'gam', formula = y ~ s(x, k = 3), size = 2, se = F, color = 'red', alpha = 0.2) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(25, 90, 5),
                     labels = paste0(Loss_age[!Loss_age$group %in% c(5, 10, 15, 20), 'group'], '\n',
                                     ' n=',
                                     Loss_age[!Loss_age$group %in% c(5, 10, 15, 20), 'n'])) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 70),
                     breaks = c(0, 50, 25)) +
  
  theme_bw(base_size = 22) +
  theme(panel.grid = element_blank()) +
  labs(x = 'AGE groups', y = 'Percent Y chromosome loss')

ggsave_golden(filename = 'Figures/WES_Age_Y_loss.pdf', plot = age_loss, width = 12)





#' Are there samples which have a Y-chromosome mosaicism and an intact Y-chromosome in tumor
#' You need to imagine that the assessment of Y-chromosome mosaicism is about a quantitative approach
#' mosaicism describes a state where several cells loose the Y-chromosome (measured in total density)
#' but we are still getting signals from the sequencing. 
#' Sequencing reads - or actually with the sliding window approach, I try to figure out the average
#' Y-chromosome DNA concentration. 
#' However, we can still use a reduced Y-chromosome quanity to determine the Y-chromosome status in tumors
#' via the T/N log ratio (however, keep in mind that the concentration of the normal will influence the estimation in the tumor)
#' e.g. if we have a normal sample with almost no coverage (much of mosaicism, I would expect the Y-call in the tumor to be amplified)



#' Starting with those samples with lowest Y-chromosome DNA content:
#' Fetch those samples (confirm manually) - and look into respective tumor sample
#' The hypothesis is, that Y-chromosome must either be called as lost in tumor or amplified

#' mosaic calls
WES_Germline = read.csv('Data_out/WES/MSK_WES_GermlineCN.txt', sep = '\t')
WES_Germline$sample = substr(x = WES_Germline$sample, start = 1, stop = 36)
WES_Germline = WES_Germline[which(WES_Germline$target == 24), ]

#' binary classifications
WES_all = cohort$WES.cohort
WES_all$id = basename(WES_all$Facet_Countfile)
WES_all$id = substr(x = WES_all$id, start = 1, stop = 36)

#' tcn's of Y-chromosomes
WES_tcn = read.csv('Data_out/WES/WES_copyNumber_out.txt', sep = '\t')
WES_tcn = WES_tcn[which(WES_tcn$ID %in% class1$id), ]
WES_tcn = WES_tcn[which(WES_tcn$chrom == 24), ]



#' CLASS 1: samples with mosaic Y-loss and an intact Y-chromosome in the tumor
#' n = 118 cases: look at those

#' Mosaic loss and Y-chromosome loss
WES_tcn$populationEstimate = NA

for(i in unique(WES_tcn$ID)){
  da = WES_tcn[which(WES_tcn$ID == i), ]
  if(nrow(da) > 1){
    WES_tcn$populationEstimate[which(WES_tcn$ID == i)] = mean(da$cnlr.median[which(da$ID == i)])
  } else {
    WES_tcn$populationEstimate[which(WES_tcn$ID == i)] = da$cnlr.median
  }
}

WES_Y_Estimation = WES_tcn[, c('ID', 'populationEstimate')]
WES_Y_Estimation = unique(WES_Y_Estimation)

class1 = WES_all[which(WES_all$Y_mosaic == 'yes' & WES_all$Y_call == 'intact_Y_chrom'), ]
class1 = merge(class1, WES_Germline, by.x = 'id', by.y = 'sample', all.x = T)
class1 = merge(class1, WES_Y_Estimation, by.x = 'id', by.y = 'ID', all.x = T)

class1.plot = ggplot(class1, aes(x = corrected.CN, y = populationEstimate)) +
  geom_jitter() +
  geom_smooth(method = 'lm', se = F) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(-1.5, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_bw(base_size = 18) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  labs(x = 'Germline Copy-Number', y = 'Tumor Copy-Number')

ggsave_golden(filename = 'Figures/GermlineTumor_correlation.pdf', plot = class1.plot, width = 8)


#' Class 4: We have 112 samples where we see some sort of mosaic and equally a Y-chromosome loss
#' in the solid tumor tissue. The question here is, do we see the Y-chromosome loss in the solid
#' tumor site just because of the mosaicism or is this a unique event associated with the cancer
class4 = WES_all$id[which(WES_all$Y_call == 'Y_chrom_loss' & WES_all$Y_mosaic == 'yes')]
class4 = WES_all[which(WES_all$id %in% class4), ]


head(class4)
WES_CNA = read.csv('Data_out/WES/WES_copyNumber_out.txt', sep = '\t')
WES_CNA = WES_CNA[which(WES_CNA$ID %in% class1), ]



