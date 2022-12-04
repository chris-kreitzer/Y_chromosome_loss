##----------------+
## mosaic loss of the Y chromosome
## Taking MSK-IMPACT normal samples;
## Asking whether chromosome Y loss
## is due to mosaic loss of Y or truly
## due to cancer
##----------------+

## start data: 08/20/2021
## modified: 09/06/2021
## revision: 07/27/2022
## revision: 12/03/2022
## chris kreitzer


## Firstly fetch the information from the cluster
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

#' Facets accounts for mappability and GC content bias;
#' Only regions with GC content between 40 and 60% are considered;

snp.nbhd = 0
normal_fetch = function(data){
  data.in = suppressMessages(data.table::fread(data))
  data.in = data.frame(Chromosome = data.in$Chromosome,
                       Position = data.in$Position,
                       NOR.DP = data.in$File1A + data.in$File1R,
                       NOR.RD = data.in$File1R,
                       TUM.DP = data.in$File2A + data.in$File2R, 
                       TUM.RD = data.in$File2R)
  
  #' only positions with unique alignments are considered
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



##----------------+
## Sliding-Window approach
##----------------+
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
  
  #' for loop over chromosomes
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
write.table(depth_out, file = '/juno/home/kreitzec/Y_chromosome_loss/Mosaicism/IMPACT_coverage_bins.txt', 
            row.names = F, sep = '\t')



##----------------+
## Downstream analysis;
##----------------+

## start: 08/16/2022
## revision: 12/03/2022
## chris-kreitzer

#' bin autosomes and allosomes 
#' and calculate the read-depth ratio
Normal_coverage = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/03_Mosaicism/IMPACT_coverage_N_22320.txt', sep = '\t')

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

LOY$sample.id = substr(LOY$sample, start = 17, stop = 33)
LOY$LOY = ifelse(LOY$corrected.CN <= quantile(LOY$corrected.CN, probs = c(0.025)), 'yes', 'no')
write.table(LOY, file = '~/Documents/MSKCC/10_MasterThesis/Data/03_Mosaicism/IMPACT_LOY.txt', sep = '\t', row.names = F)


##-----------------
## Association studies: Age
##----------------- 
LOY = readRDS('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')
LOY = LOY$IMPACT_LOY
LOY = LOY[!LOY$Age_Sequencing %in% c('0', '1', '2', '3', '4', '5', '6', '7'), ]
LOY = LOY[!is.na(LOY$Age_Sequencing), ]
LOY$Age_Sequencing = as.integer(as.character(LOY$Age_Sequencing))
LOY$mLRR_Y = as.numeric(as.character(LOY$mLRR_Y))

#' Visualization
Y_chromosome = ggplot(LOY, aes(x = Age_Sequencing, y = mLRR_Y)) +
  geom_jitter(size = 0.2) +
  scale_y_continuous(limits = c(0.5, 1.5),
                     breaks = c(0.5, 1, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.05),
                     breaks = seq(10, 90, 10)) +
  geom_smooth(method = 'lm') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, 
                                    size = 2),
        axis.text = element_text(size = 12, 
                                 color = 'black')) +
  annotate(geom = 'text',
           x = 25,
           y = 1.40,
           label = 'y = 0.995-0.00024x\np < 2e-16',
           size = 5.5) +
  labs(x = 'Age [reported at sequencing]', y = 'mLRR-Y', title = 'MSK-IMPACT: mLRR-Y decrease with age')

Y_chromosome


##-----------------
## Sanity check for X-chromosome
##-----------------
germline = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/03_Mosaicism/IMPACT_mLRR-Y_summary.txt', sep = '\t')
germlineX = germline[which(germline$target == 'X'), ]
germlineX$sample.id = substr(germlineX$sample, start = 17, stop = 33)
cohort = readRDS('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/Cohort_07132022.rds')
cohort = cohort$IMPACT_clinicalAnnotation

germlineX = merge(germlineX, cohort[, c('SAMPLE_ID', 'AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)')],
                  by.x = 'sample.id', by.y = 'SAMPLE_ID', all.x = T)

germlineX = germlineX[!germlineX$`AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)` %in% c('0', '1', '2', '3', '4', '5', '6', '7'), ]
germlineX = germlineX[!is.na(germlineX$`AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)`), ]
germlineX$`AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)` = as.integer(as.character(germlineX$`AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)`))
colnames(germlineX)[8] = 'Age_Sequencing'

X_chromosome = ggplot(germlineX, aes(x = Age_Sequencing, y = corrected.CN)) +
  geom_jitter(size = 0.2) +
  scale_y_continuous(limits = c(0.5, 1.5),
                     breaks = c(0.5, 1, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.05),
                     breaks = seq(10, 90, 10)) +
  geom_smooth(method = 'lm') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, 
                                    size = 2),
        axis.text = element_text(size = 12, 
                                 color = 'black')) +
  annotate(geom = 'text',
           x = 25,
           y = 1.40,
           label = 'y = 1.017+0.0000327x\np = 0.06',
           size = 5.5) +
  labs(x = 'Age [reported at sequencing]', y = 'mLRR-X', title = 'MSK-IMPACT: stable X chromosome')

X_chromosome  


library(patchwork)
ggsave_golden(filename = 'Figures_original/Mosaic_AGE_IMPACT.pdf', plot = (Y_chromosome / X_chromosome), width = 7.5)

#' out