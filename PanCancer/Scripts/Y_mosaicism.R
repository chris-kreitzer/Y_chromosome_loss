## Investigate Y-chromosome mosaicism:
## This should help to distinguish whether Y-chromosome loss happened
## as a physiological cuase (age) or truly because of the cancer.
## 
## start data: 08/20/2021
## chris kreitzer

## I will start working on whole exome sequenced samples;
## Firstly fetch the information from the cluster

## install local FacetsY and pctGCdata
require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
source('/juno/home/kreitzec/WES_Prostate/FacetsQC.R')
library(facetsSuite)
library(dplyr)
library(data.table)

# load FACETS countsfiles for Prostate WES samples
IMPACT.prostate.samples = read.csv('/juno/home/kreitzec/WES_Prostate/Panel.prostate.samplepath.txt', header = F)
IMPACT.prostate.samples = as.character(IMPACT.prostate.samples$V1)

#' mappability issues and GC content are already considered in Facets pipeline
#' I will only include regions with GC content between 45 and 55%; equal distribution across the genome.
snp.nbhd = 0
normal_coverage_df_prostate = data.frame()
for(i in unique(IMPACT.prostate.samples)){
  try({
    print(i)
    data.in = facetsY::readSnpMatrix(i)
    data.pre = facetsY::preProcSample(data.in, 
                                      gbuild = 'hg19', 
                                      snp.nbhd = snp.nbhd)
    data.pre = data.pre$jointseg[which(data.pre$jointseg$gcpct >= 0.45 & data.pre$jointseg$gcpct <= 0.55),
                                 c('chrom', 'maploc', 'rCountN', 'gcpct')]
    data.pre$sample = i
    
    normal_coverage_df_prostate = rbind(normal_coverage_df_prostate , data.pre)
    
  })
}

write.table(normal_coverage_df_prostate, file = '/juno/home/kreitzec/Y_chromosome_loss/Mosaicism/WES_normal_coverage_Prostate.txt', row.names = F, sep = '\t')



## Start with the downstream analysis;
## the main hypothesis is, that the Y-chromosome DNA content is roughly half of that of autosomes 
## I claim that I can estimate the Y-chromosome copy number from the normal sequencing results via the read depth
## I will calculate the median normal coverage among the autosomes and compare to the allosomes
## the ratio allosomes / median autosomes will inform us, whether patient has an intact Y-chromosome or not
## look for age dependence

#' sliding window approach with evoBiR (1kb windows; median depth of sequencing in the normals)
library(evobiR)

## Input
Normals_Prostate = read.csv('../Data_out/WES_normal_coverage_Prostate.txt', sep = '\t')
Normals_Prostate$sample = basename(Normals_Prostate$sample)
Annotation = read.table('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/WES_clinical_annotation.txt', sep = '\t', quote = "", header = T)
WES_paths = read.csv('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_in/WES_paths.txt', header = F)
WES_paths = paste0('/juno/work/tempo/wes_repo/Results/v1.3.x/cohort_level/MSKWESRP/somatic/', WES_paths$V1, '/facets/', WES_paths$V1, '/', WES_paths$V1, '.snp_pileup.gz')
write.table(WES_paths, file = '~/Documents/GitHub/Y_chromosome_loss/PanCancer/Data_out/WES_facetspaths.txt', col.names = F, row.names = F, quote = F)


depth_out = data.frame()
for(i in unique(Normals_Prostate$sample)){
  data.selected = Normals_Prostate[which(Normals_Prostate$sample == i), ]
  print(i)
  for(j in seq_along(data.selected$chrom)){
    if(length(data.selected$rCountN[which(data.selected$chrom == j)]) > 1000){
      median_depth = evobiR::SlidingWindow(FUN = median,
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


#' make summary over all autosomes / sample and compare to allosomes
chromo_out = data.frame()

for(i in unique(depth_out$sample)){
  if(length(unique(depth_out$chrom[which(depth_out$sample == i)])) == 24){
    print('yes')
    patient = depth_out[which(depth_out$sample == i), ]
    for(chromo in unique(patient$chrom)){
      median_genome = median(patient$depth[!patient$chrom %in% chromo & patient$sample == i])
      median_target = median(patient$depth[which(patient$chrom == chromo & patient$sample == i)])
      summary_df = data.frame(sample = i,
                              target = chromo,
                              ratio = median_target / median_genome)
      chromo_out = rbind(chromo_out, summary_df)
    }
  } else next
}

#' make a ploidy correction through population-wise determination of the observed (peak) to 
#' expected ploidy level across chromosomes.

chromo_out$CN = chromo_out$ratio * 2
chromo_out$corrected.CN = NA
for(i in unique(chromo_out$target)){
  density.chromo = density(chromo_out$CN[which(chromo_out$target == i)], bw = 'SJ')
  print(density.chromo$x[which.max(density.chromo$y)])
  density.max = density.chromo$x[which.max(density.chromo$y)]
  corrected.CN = ifelse(i <= 22, (density.max - 2), (density.max - 1))
  print(corrected.CN)
  for(patient in unique(chromo_out$sample)){
    chromo_out$corrected.CN[which(chromo_out$target == i & chromo_out$sample == patient)] = chromo_out$CN[which(chromo_out$sample == patient & chromo_out$target == i)] - corrected.CN 
  }
}

chromo_out$target = factor(chromo_out$target, levels = seq(1, 24, 1))
library(ggplot2)
ggplot(chromo_out, aes(x = target, y = corrected.CN)) + geom_boxplot()

## calculate the CI interval
normConfInt = function(x, alpha = 0.05){
  mean(x) + qt(1 - alpha / 2, length(x) - 1) * sd(x) / sqrt(length(x)) * c(-1, 1)
}
  
normConfInt(x = chromo_out$corrected.CN[which(chromo_out$target == 24)], 0.1)
summary(chromo_out$corrected.CN[which(chromo_out$target == 24)])

hist(chromo_out$corrected.CN[which(chromo_out$target == 24)], nclass = 50)
lines(density(chromo_out$corrected.CN[which(chromo_out$target == 24)]))





