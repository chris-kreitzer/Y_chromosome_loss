## Investigate Y-chromosome mosaicism: Here I am looking into IMPACT samples
## This should help to distinguish whether Y-chromosome loss happened
## as a physiological cause (e.g., age) or truly because of the cancer.
## 
## start data: 08/20/2021
## modified: 09/06/2021
## revision: 07/27/2022
## chris kreitzer

## Firstly fetch the information from the cluster
## install local FacetsY and pctGCdata
require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('data.table', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
library(facetsSuite)
library(dplyr)
library(data.table)
set.seed(99)

# load FACETS countsfiles for Prostate WES samples
IMPACT.samples = read.csv('/juno/home/kreitzec/Master/Data/cohort_data.rds')
IMPACT.samples = as.character(IMPACT.samples$IMPACT.cohort$counts_file)

#' mappability issues and GC content are already considered in Facets pipeline
#' I will only include regions with GC content between 40 and 60%; equal distribution across the genome.

#' start processing the normal samples
snp.nbhd = 0
normal_fetch = function(data){
  data.in = suppressMessages(data.table::fread(data))
  data.in = data.frame(Chromosome = data.in$Chromosome,
                       Position = data.in$Position,
                       NOR.DP = data.in$File1A + data.in$File1R,
                       NOR.RD = data.in$File1R,
                       TUM.DP = data.in$File2A + data.in$File2R, 
                       TUM.RD = data.in$File2R)
  
  #' consider only postions where we have unique sequence alignments
  Y_chromosome = data.in[which(data.in$Chromosome == 'Y' & 
                                 data.in$Position >= 2654550 & 
                                 data.in$Position <= 28000000), ]
  
  #' exclude PCDH11Y and centromeric region
  Y_chromosome = Y_chromosome[with(Y_chromosome, !((Position %in% 4922131:5612269))), ]
  Y_chromosome = Y_chromosome[with(Y_chromosome, !((Position %in% 10500000:14000000))), ]
  
  ## autosomes:
  autosomes = rbind(data.in[which(data.in$Chromosome %in% c(seq(1, 22, 1), 'X')), ], Y_chromosome)
  
  ## run analysis
  data.pre = facetsY::preProcSample(rcmat = autosomes,
                                    ndepth = 30,
                                    snp.nbhd = snp.nbhd, 
                                    gbuild = 'hg19', 
                                    ndepthmax = 1000)
  data.pre = data.pre$jointseg[which(data.pre$jointseg$gcpct >= 0.4 & data.pre$jointseg$gcpct <= 0.6),
                               c('chrom', 'maploc', 'rCountN')]
  data.pre$sample = basename(data)
  data.pre
}

Normal_out = lapply(IMPACT.samples, function(x) normal_fetch(x))
message('rbindlist')
Normal_out = data.table::rbindlist(Normal_out)


## SlidingWindow approach
SlidingWindow = function(FUN, data, window, step) {
  total = length(data)
  spots = seq(from = 1, to = (total - window), by = step)
  result = vector(length = length(spots))
  for (i in 1:length(spots)) {
    result[i] = match.fun(FUN)(data[spots[i]:(spots[i] + 
                                                window - 1)])
  }
  return(result)
}

message('starting SlidingWindow now')

#' determine median read-depth across 500 bp intervals
determine_depth = function(data){
  data.selected = Normal_out[which(Normal_out$sample == data), ]
  print(basename(data))
  
  #' for loop work over chromosomes
  out_all = data.frame()
  for(j in unique(data.selected$chrom)){
    if(length(data.selected$rCountN[which(data.selected$chrom == j)]) > 500){
      median_depth = SlidingWindow(FUN = median,
                                   data = data.selected$rCountN[which(data.selected$chrom == j)],
                                   window = 500,
                                   step = 500)
    } else if(length(data.selected$rCountN[which(data.selected$chrom == j)]) <= 500 & 
              length(data.selected$rCountN[which(data.selected$chrom == j)]) > 30){
      median_depth = median(data.selected$rCountN[which(data.selected$chrom == j)])
    } else {
      median_depth = NA
    }
    
    out = data.frame(chrom = j,
                     sample = data,
                     depth = median_depth)
    out_all = rbind(out_all, out)
  }
  out_all
}

depth_out = lapply(unique(Normal_out$sample), function(x) determine_depth(x))
depth_out = data.table::rbindlist(depth_out)


write.table(depth_out, file = '/juno/home/kreitzec/Y_chromosome_loss/Mosaicism/IMPACT_coverage_bins.txt', row.names = F, sep = '\t')



#######################
#' Downstream analysis:
#' 
#' start: 08/16/2022
#' chris-kreitzer
#' bin autosomes and allosomes 
#' and calculate the read-depth ratio

#' write.table(mosaic, file = '~/Documents/MSKCC/10_MasterThesis/Data/03_Mosaicism/IMPACT_coverage_N_22320.txt', sep = '\t', quote = F, row.names = F)
mosaic = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/03_Mosaicism/IMPACT_coverage_N_22320.txt', sep = '\t')
Normal_coverage = mosaic

bins_summary = function(data){
  print(data)
  data = Normal_coverage[which(Normal_coverage$sample == data), ]
  summary_df_out = data.frame()
  if(length(unique(data$chrom)) == 24){
    for(chromo in unique(data$chrom)){
      median_genome = median(data$depth[!data$chrom %in% chromo])
      median_target = median(data$depth[which(data$chrom == chromo)])
      summary_df = data.frame(sample = unique(data$sample),
                              target = chromo,
                              ratio = median_target / median_genome,
                              log_ratio = log(median_target) / log(median_genome))
      summary_df_out = rbind(summary_df_out, summary_df)
    }

  } else {
    return()
  }
  summary_df_out
}

sample_summary = lapply(unique(Normal_coverage$sample), function(x) bins_summary(x))
sample_summary = data.table::rbindlist(sample_summary)


#' make a ploidy correction through population-wise determination of the observed (peak) to
#' expected ploidy level across chromosomes.
sample_summary$target = as.integer(as.character(sample_summary$target))
sample_summary$CN = sample_summary$log_ratio * 2
sample_summary$corrected.CN = NA

CN_correction = data.frame()
for(i in unique(sample_summary$target)){
  density.chromo = density(sample_summary$CN[which(sample_summary$target == i)], bw = 'SJ', na.rm = T)
  density.max = density.chromo$x[which.max(density.chromo$y)]
  corrected.CN = ifelse(i %in% seq(1, 22, 1), (density.max - 2), (density.max - 1))
  correction_factor = data.frame(chromosome = i,
                                 correction_factor = corrected.CN)
  CN_correction = rbind(CN_correction, correction_factor)
}

#' correct the CN values by their density as calculated above;
for(i in unique(CN_correction$chromosome)){
  sample_summary$corrected.CN[which(sample_summary$target == i)] = sample_summary$CN[which(sample_summary$target == i)] -
    CN_correction$correction_factor[which(CN_correction$chromosome == i)]
}


#' Visualization:
sample_summary$target[which(sample_summary$target == 23)] = 'X'
sample_summary$target[which(sample_summary$target == 24)] = 'Y'
write.table(sample_summary, file = '~/Documents/MSKCC/10_MasterThesis/Data/03_Mosaicism/IMPACT_mLRR-Y_summary.txt', sep = '\t', row.names = F, quote = F)
sample_summary$target = factor(sample_summary$target, levels = c(seq(1, 22, 1), 'X', 'Y'))

Germline_CN = ggplot(sample_summary, aes(x = target, y = corrected.CN)) +
  geom_boxplot() +
  geom_hline(yintercept = seq(1, 3, 1), 
             color = 'grey35', 
             linetype = 'dashed', 
             size = 0.35) +
  scale_y_continuous(expand = c(0.01, 0.05),
                     limits = c(0, 3)) +
  labs(x = 'Chromosome', 
       y = 'Germline copy number',
       title = paste0('MSK-IMPACT samples\nn=', length(unique(sample_summary$sample))))


Germline_CN

ggsave_golden(filename = '~/Documents/MSKCC/10_MasterThesis/Figures_original/IMPACT_GermlineCN.pdf', 
              plot = Germline_CN, width = 8,
              dpi = 800, device = 'pdf')



##-----------------
## assign Mosaic loss:
##-----------------
sample_summary = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/03_Mosaicism/IMPACT_mLRR-Y_summary.txt', sep = '\t')
LOY = sample_summary[which(sample_summary$target == 'Y'), ]
LOY = LOY[!is.na(LOY$corrected.CN), ]
LOY$seq = seq(1, nrow(LOY), 1)

plot(LOY$corrected.CN, LOY$seq, 
     pch = 17, 
     cex = 0.3,
     yaxt = 'n',
     xlab = '',
     ylab = '')
axis(side = 2, at = c(0, nrow(LOY)), las = 2)
abline(v = median(LOY$corrected.CN), col = 'red', lwd = 2)
#' lower 2.5% quantile
abline(v = quantile(LOY$corrected.CN, probs = c(0.025)) , col = 'grey35')
#' upper 97.5% quantile
abline(v = quantile(LOY$corrected.CN, probs = c(0.975)), col = 'grey35')
box(lwd = 2)
mtext(text = 'MSK-IMPACT observed mLRR-Y values', adj = 0, line = 0.5)
mtext(text = 'Individuals', side = 2, line = 2)
mtext(text = 'mLRR-Y', side = 1, line = 2)


##-----------------
## Association studies:
## AGE
Clinical_annotation = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/mskimpact_clinical_data_07112022.tsv', sep = '\t')
Clinical_annotation = Clinical_annotation[,c('Sample.ID', 'Age.at.Which.Sequencing.was.Reported..Years.')]
colnames(Clinical_annotation) = c('Sample.ID', 'Age')

LOY$sample.id = substr(LOY$sample, start = 17, stop = 33)
LOY = merge(LOY, Clinical_annotation, by.x = 'sample.id', by.y = 'Sample.ID', all.x = T)
LOY = LOY[!LOY$Age %in% c('0', '>90', '1', '2', '3', '4', '5', '6', '7'), ]
LOY = LOY[!is.na(LOY$Age), ]
LOY$Age = as.integer(as.character(LOY$Age))
summary(lm(LOY$corrected.CN ~ LOY$Age))





IMPACT_cohort = readRDS('Data_out/cohort_data.rds')$IMPACT.cohort
sample_summary$sample = substr(sample_summary$sample, start = 17, stop = 33)

IMPACT_age = merge(sample_summary, IMPACT_cohort[, c('SAMPLE_ID', 'AGE_AT_SEQ_REPORTED_YEARS')],
                   by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)

IMPACT_age = IMPACT_age[!IMPACT_age$AGE_AT_SEQ_REPORTED_YEARS %in% c('>90', '', '1', '2', '3', '4', '5', '6', '7'), ]
IMPACT_age$AGE_AT_SEQ_REPORTED_YEARS = as.integer(as.character(IMPACT_age$AGE_AT_SEQ_REPORTED_YEARS))
IMPACT_age = IMPACT_age[order(IMPACT_age$AGE_AT_SEQ_REPORTED_YEARS, decreasing = F), ]
IMPACT_age$Y_mosaicism = NA
IMPACT_age$Y_mosaicism[which(IMPACT_age$target == 'Y')] = ifelse(IMPACT_age$corrected.CN[which(IMPACT_age$target == 'Y')] < lower_Y_threshold[1], 'yes', 'no')
write.table(IMPACT_age, file = 'Data_out/IMPACT/GermlineCN.txt', sep = '\t', row.names = F, quote = F)


#' modelling the association between age and Y chromsome copy number
summary(lm(IMPACT_age$corrected.CN[which(IMPACT_age$target == 'Y')] ~
             IMPACT_age$AGE_AT_SEQ_REPORTED_YEARS[which(IMPACT_age$target == 'Y')]))

#' X-chromosome visualization
X_chromosome = ggplot(IMPACT_age[which(IMPACT_age$target == 'X'), ],
                      aes(x = AGE_AT_SEQ_REPORTED_YEARS, y = corrected.CN)) +
  geom_jitter(size = 0.2) +
  scale_y_continuous(limits = c(0, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.05),
                     breaks = seq(10, 90, 10)) +
  geom_smooth(method = 'lm') +
  annotate(geom = 'text',
           x = 20,
           y = 0.2,
           label = 'y = 1.08 -0.0002 x',
           size = 5.5) +
  labs(x = 'Age', y = 'Chr. X copy number')

#' Y-chromosome visualization
Y_chromosome = ggplot(IMPACT_age[which(IMPACT_age$target == 'Y'), ],
                      aes(x = AGE_AT_SEQ_REPORTED_YEARS, y = corrected.CN)) +
  geom_jitter(size = 0.2) +
  scale_y_continuous(limits = c(0, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.05),
                     breaks = seq(10, 90, 10)) +
  geom_smooth(method = 'lm') +
  annotate(geom = 'text',
           x = 20,
           y = 0.2,
           label = 'y = 0.966 -0.00052 x\np < 2e-16',
           size = 5.5) +
  labs(x = 'Age', y = 'Chr. Y copy number')


library(patchwork)
ggsave_golden(filename = 'Figures/Mosaic_AGE_IMPACT.pdf', plot = (Y_chromosome / X_chromosome), width = 7.5)


