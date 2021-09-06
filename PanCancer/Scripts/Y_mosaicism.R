## Investigate Y-chromosome mosaicism:
## This should help to distinguish whether Y-chromosome loss happened
## as a physiological cause (e.g., age) or truly because of the cancer.
## Importantly, note that you are analyzing mostly blood samples - normal coverage
## and this is not directly linked to the tissue of cancer you are investigating. 
## So, mosaicism of the blood is not necessarily tissue mosaicism of the e.g. the pancreas
## 
## start data: 08/20/2021
## chris kreitzer

setwd('~/Documents/GitHub/Y_chromosome_loss/PanCancer/')

## I will start working on whole exome sequenced samples;
## Firstly fetch the information from the cluster

## install local FacetsY and pctGCdata
require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('vroom', '/juno')
library(facetsSuite)
library(dplyr)
library(data.table)
library(ggplot2)

# load WES_paths (only male samples to assess Y_chromosome mosaicism)
WES.samples = readRDS('/juno/home/kreitzec/Y_chromosome_loss/Data/cohort_data.rds')
WES.samples = as.character(WES.samples$WES_male$Facet_Countfile)

## the main hypothesis is, that the Y-chromosome DNA content is roughly half of that of autosomes 
## I claim, that I can estimate the Y-chromosome copy number from the normal sequencing results via the read depth
## I will calculate the median normal coverage among the autosomes and compare to the allosomes
## the ratio allosomes / median autosomes will inform us, whether patient has an intact Y-chromosome or not
## look for age dependence

#' mappability issues and GC content are already considered in Facets pipeline
#' I will only include regions with GC content between 45 and 55%; equal distribution across the genome.
print(length(unique(WES.samples)))

#' start processing the normal samples
snp.nbhd = 250
normal_coverage_df = data.frame()
for(i in unique(WES.samples)){
  try({
    print(i)
    data.in = suppressMessages(vroom::vroom(i))
    data.in_summary = data.frame(Chromosome = data.in$Chromosome,
                                 Position = data.in$Position,
                                 NOR.DP = data.in$File1A + data.in$File1R,
                                 NOR.RD = data.in$File1R,
                                 TUM.DP = data.in$File2A + data.in$File2R,
                                 TUM.RD = data.in$File2R)
    
    data.pre = facetsY::preProcSample(data.in_summary, 
                                      gbuild = 'hg19', 
                                      snp.nbhd = snp.nbhd)
    data.pre = data.pre$jointseg[which(data.pre$jointseg$gcpct >= 0.45 & data.pre$jointseg$gcpct <= 0.55),
                                 c('chrom', 'maploc', 'rCountN')]
    data.pre$sample = i
    normal_coverage_df = rbind(normal_coverage_df, data.pre)
    
  })
}

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


depth_out = data.frame()
for(i in unique(normal_coverage_df$sample)){
  data.selected = normal_coverage_df[which(normal_coverage_df$sample == i), ]
  print(i)
  for(j in seq_along(data.selected$chrom)){
    if(length(data.selected$rCountN[which(data.selected$chrom == j)]) > 1000){
      median_depth = SlidingWindow(FUN = median,
                                   data = data.selected$rCountN[which(data.selected$chrom == j)],
                                   window = 1000,
                                   step = 1000)
    } else if(length(data.selected$rCountN[which(data.selected$chrom == j)]) <= 1000 & 
              length(data.selected$rCountN[which(data.selected$chrom == j)]) > 100){
      median_depth = median(data.selected$rCountN[which(data.selected$chrom == j)])
    } else next
    
    out = data.frame(chrom = j,
                     sample = i,
                     depth = median_depth)
    depth_out = rbind(depth_out, out)
  }
}

write.table(normal_coverage_df, file = '/juno/home/kreitzec/Y_chromosome_loss/Mosaicism/WES_normal_coverage.txt', row.names = F, sep = '\t')
write.table(depth_out, file = '/juno/home/kreitzec/Y_chromosome_loss/Mosaicism/Normal_coverage_bins.txt', row.names = F, sep = '\t')



###############################################################################
###############################################################################
##' Downstream analysis; investigate the output from above
#' make summary over all autosomes / sample and compare to allosomes
Normal_coverage = read.csv('Mosaicism/Normal_coverage_bins.txt', sep = '\t')
Normal_coverage$sample = basename(Normal_coverage$sample)

sample_summary = data.frame() 
for(i in unique(Normal_coverage$sample)){
  if(length(unique(Normal_coverage$chrom[which(Normal_coverage$sample == i)])) == 24){
    patient = Normal_coverage[which(Normal_coverage$sample == i), ]
    for(chromo in unique(patient$chrom)){
      median_genome = median(patient$depth[!patient$chrom %in% chromo & patient$sample == i])
      median_target = median(patient$depth[which(patient$chrom == chromo & patient$sample == i)])
      summary_df = data.frame(sample = i,
                              target = chromo,
                              ratio = median_target / median_genome)
      sample_summary = rbind(sample_summary, summary_df)
    }
  } else next
}

#' make a ploidy correction through population-wise determination of the observed (peak) to 
#' expected ploidy level across chromosomes.
sample_summary$target = as.integer(as.character(sample_summary$target))
sample_summary$CN = sample_summary$ratio * 2
sample_summary$corrected.CN = NA
CN_correction = data.frame()
for(i in unique(sample_summary$target)){
  density.chromo = density(sample_summary$CN[which(sample_summary$target == i)], bw = 'SJ')
  density.max = density.chromo$x[which.max(density.chromo$y)]
  corrected.CN = ifelse(i %in% seq(1, 22, 1), (density.max - 2), (density.max - 1))
  correction_factor = data.frame(chromosome = i,
                                 correction_factor = corrected.CN)
  CN_correction = rbind(CN_correction, correction_factor)
}

#' correct the CN values by their density as calculated above;
for(i in unique(CN_correction$chromosome)){
  sample_summary$corrected.CN[which(sample_summary$target == i)] = sample_summary$CN[which(sample_summary$target == i)] - CN_correction$correction_factor[which(CN_correction$chromosome == i)]
}


#' Visualization:
sample_summary$target = factor(sample_summary$target, levels = seq(1, 24, 1))
Germline_CN = ggplot(sample_summary, aes(x = target, y = corrected.CN)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = 'chromosome', y = 'corrected copy Number', 
       title = paste0('Estimated germline copy number\ndetermined from read depth data\nn=', length(unique(sample_summary$sample)), ' WES samples'))

ggsave_golden(plot = Germline_CN, filename = 'Figures/Germline_CN.pdf', width = 8)



#' define the LOY mosaicism threshold based on the 99% CI
normConfInt = function(x, alpha = 0.05){
  mean(x) + qt(1 - alpha / 2, length(x) - 1) * sd(x) / sqrt(length(x)) * c(-1, 1)
}


lower_Y_threshold = normConfInt(x = sample_summary$corrected.CN[which(sample_summary$target == 24)])
length(sample_summary$sample[which(sample_summary$target == 24 & sample_summary$corrected.CN < 0.8189)])


#' out