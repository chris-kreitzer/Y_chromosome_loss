## PanCancer Analysis on Y-chromosome loss:
## 
## look into the average amount of protein coding genes which are 
## directly covered by WES- and IMPACT-Seq.
## 
## Start: 14/04/2021


rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')


## Libraries
require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
library(dplyr)
library(data.table)


## Input
Features.Y_chromosome = read.csv('Data_out/Y_chromosome_genes.txt', sep = '\t')
Sample.path = read.csv('/juno/home/kreitzec/WES_Prostate/Panel.prostate.samplepath.txt', header = F)
Sample.path = as.character(Sample.path$V1)
print(Sample.path)

## Processing:
Protein.coding.genes = Features.Y_chromosome[which(Features.Y_chromosome$gene_biotype == 'protein_coding'), 
                                         c('hgnc_symbol', 'start_position', 'end_position', 'external_gene_name')]



## Analysis:
number.protein.coding.genes = nrow(Protein.coding.genes)

cval.preprocess = 75
cval.postprocess = 100
snp.nbhd = 250

Features.covered_out = data.frame()
Features.covered_summary = data.frame()

for(i in unique(Sample.path)){
  try({
    print(i)
    ID = substr(i, start = 61, stop = 90)
    data.in = facetsY::readSnpMatrix(i)
    data.pre = facetsY::preProcSample(data.in, 
                                      cval = cval.preprocess, 
                                      gbuild = 'hg19', 
                                      snp.nbhd = snp.nbhd)
    data.process = facetsY::procSample(data.pre, 
                                       cval = cval.postprocess)
    
    data.analysis = data.process$jointseg
    data.analysis = data.analysis[which(data.analysis$chrom == 24), ]
    
    #' look into Feature coverage
    for(i in 1:nrow(Protein.coding.genes)){
      if(any(between(x = data.analysis$maploc, lower = Protein.coding.genes$start_position[i], 
                     upper = Protein.coding.genes$end_position[i]))){
        out = data.frame(Gene = Protein.coding.genes$external_gene_name[i],
                         Match = 1,
                         number = sum(between(x = data.analysis$maploc, 
                                              lower = Protein.coding.genes$start_position[i], 
                                              upper = Protein.coding.genes$end_position[i])),
                         ID = ID)
        
      } else next
      Features.covered_out = rbind(Features.covered_out, out)
    }
    
    #' summary stats per sample
    sample.summary = data.frame(Sample = ID,
                                Features_covered = sum(out$Match) / number.protein.coding.genes)
    Features.covered_summary = rbind(Features.covered_summary, sample.summary)
  })
    
}

write.table(Features.covered_out, file = '/juno/home/kreitzec/WES_Prostate/IMPACT.Features.covered.out.txt', row.names = F, sep = '\t')
write.table(Features.covered_summary, file = '/juno/home/kreitzec/WES_Prostate/IMPACT_Features.covered.summary.txt', sep = '\t', row.names = F)
