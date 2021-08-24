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
Y_mosaic_df = data.frame()
for(i in unique(depth_out$sample)){
  if(length(depth_out$depth[which(depth_out$sample == i & depth_out$chrom == 24)]) != 0){
    autosomes.median = median(depth_out$depth[which(depth_out$chrom %in% seq(1, 23, 1) & depth_out$sample == i)])
    allosomes.median = median(depth_out$depth[which(depth_out$chrom == 24 & depth_out$sample == i)])
    Y_mosaic = data.frame(Sample = i,
                          autosomes.median = autosomes.median,
                          allosomes.median = allosomes.median)
  } else next
  
  Y_mosaic_df = rbind(Y_mosaic_df, Y_mosaic)
}


#' investigate the results:
#' the main hypothesis is the following:
#' The older patient gets, the more likely the suffer from a Y-chromosome loss.
#' We use the median read depth ratio as proxy
Y_mosaic_df$Sample = substr(x = Y_mosaic_df$Sample, start = 17, stop = 33)
Y_mosaic_df$ratio = (Y_mosaic_df$allosomes.median / Y_mosaic_df$autosomes.median) * 2
Y_ratio_density = density(Y_mosaic_df$ratio, bw = 'SJ')
Y_max = Y_ratio_density$x[which.max(Y_ratio_density$y)]
Y_mosaic_df$ratio.corrected = Y_mosaic_df$ratio - Y_max

#' merge with clinical data:
x = merge(Y_mosaic_df, Annotation[c('Sample.ID', 'Age.at.Which.Sequencing.was.Reported..Years.')],
          by.x = 'Sample', by.y = 'Sample.ID', all.x = T)

x = x[order(x$Age.at.Which.Sequencing.was.Reported..Years., decreasing = T), ]
View(x)






plot(density(Y_mosaic_df$ratio.corrected))
head(Y_mosaic_df)


a = Normals_Prostate[which(Normals_Prostate$Sample == 'countsMerged____P-0021240-T01-IM6_P-0021240-N01-IM6.dat.gz'), ]

al = data.frame()
for(i in unique(a$chrom)){
  b = data.frame(chrom = i,
                 median = median(a$rCountN[which(a$chrom == i)]))
  al = rbind(al, b)
}

plot(al$chrom, al$median)
abline(h = 383)
median(al$median[which(al$chrom %in% seq(1, 23, 1))])
102/383 * 2



Y_mosaic_df$ratio = 
plot(Y_mosaic_df$ratio)

confint(Y_mosaic_df$ratio)

abline(h = 0.4)
plot(density(Y_mosaic_df$ratio))
mad(Y_mosaic_df$ratio)







