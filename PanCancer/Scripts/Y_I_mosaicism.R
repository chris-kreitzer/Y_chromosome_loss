## Investigate Y-chromosome mosaicism: Here I am looking into IMPACT samples
## This should help to distinguish whether Y-chromosome loss happened
## as a physiological cause (e.g., age) or truly because of the cancer.
## 
## start data: 08/20/2021
## modified: 09/06/2021
## chris kreitzer

## I will start working on whole exome sequenced samples;
## Firstly fetch the information from the cluster

## install local FacetsY and pctGCdata
require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('vroom', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('data.table', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
library(facetsSuite)
library(dplyr)
library(data.table)

# load FACETS countsfiles for Prostate WES samples
IMPACT.samples = readRDS('/juno/home/kreitzec/Y_chromosome_loss/Data/cohort_data.rds')
IMPACT.samples = as.character(IMPACT.samples$IMPACT.cohort$counts_file)

#' mappability issues and GC content are already considered in Facets pipeline
#' I will only include regions with GC content between 45 and 55%; equal distribution across the genome.
print(length(unique(IMPACT.samples)))


#' start processing the normal samples
snp.nbhd = 0
normal_fetch = function(data){
  data.in = suppressMessages(vroom::vroom(data))
  data.in = data.frame(Chromosome = data.in$Chromosome,
                       Position = data.in$Position,
                       NOR.DP = data.in$File1A + data.in$File1R,
                       NOR.RD = data.in$File1R,
                       TUM.DP = data.in$File2A + data.in$File2R,
                       TUM.RD = data.in$File2R)
  data.pre = facetsY::preProcSample(rcmat = data.in,
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



###############################################################################
#' Downstream analysis:
#' bin autosomes and allosomes and calculate the read-depth ratio
data = read.csv('Mosaicism/IMPACT_coverage_bins.txt', sep = '\t')
Normal_coverage = data

bins_summary = function(data){
  data = Normal_coverage[which(Normal_coverage$sample == data), ]
  summary_df_out = data.frame()
  if(length(unique(data$chrom)) == 24){
    for(chromo in unique(data$chrom)){
      median_genome = median(data$depth[!data$chrom %in% chromo])
      median_target = median(data$depth[which(data$chrom == chromo)])
      summary_df = data.frame(sample = unique(data$sample),
                              target = chromo,
                              ratio = median_target / median_genome)
      summary_df_out = rbind(summary_df_out, summary_df)
    }
    
  } else {
    return()
  }
  summary_df_out
}

sample_summary = lapply(unique(Normal_coverage$sample), function(x) bins_summary(x))
sample_summary = data.table::rbindlist(sample_summary)







