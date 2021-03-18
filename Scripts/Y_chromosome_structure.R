## PanCancer analysis of Y chromosome loss:
## Structure of the Y_chromosome:
## 
## The initial analysis started in fall 2020; 
## Initial scripts and pre-analysis were run by Bastien and Subhi (back in 2019 and 2020)
## 
## Here, I will revise some concepts and provide technical validation of Y_chromosome_loss calls
## 
## Start (revision): 11/03/2021:
## Scripts will be hosted on Github: https://github.com/chris-kreitzer/Y_chromosome_loss.git

rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')


## Libraries
library(pctGCdata)
devtools::install_github('https://github.com/BastienNguyen/facetsY.git', 
                         build_vignettes = F, force = T)
library(facetsY)
library(tidyverse)
library(cowplot)
library(statebins)
library(tempR)

## Input data:

# Pileup = vroom::vroom('~/Documents/MSKCC/04_Y_chromo_loss/Data/FacetsPileupFile', delim = '\t', skip = 153)
# colnames(Pileup)[1] = 'Chrom'
# Pileup.Y_df = Pileup[which(Pileup$Chrom == 'chrY'), ]
# write.table(Pileup.Y_df, file = '~/Documents/MSKCC/04_Y_chromo_loss/Data/Y_chromosome_Facets_Pileup.txt', sep = '\t', row.names = F)

# MSK.WES_df = read.csv('~/Documents/MSKCC/04_Y_chromo_loss/s_C_740W28_M001_d__s_C_740W28_N001_d.snp_pileup.gz', sep = ',')
# MSK.WES_Y_df = MSK.WES_df[which(MSK.WES_df$Chromosome == 'Y'), ]
# write.table(MSK.WES_Y_df, file = '~/Documents/MSKCC/04_Y_chromo_loss/Data/Y_chromosome_IMPACT_WES.example.txt', sep = '\t', row.names = F)

# MSK.Panel_df = read.csv('~/Desktop/mnt/ATMcountdata/countsMerged____P-0017892-T02-IM6_P-0017892-N01-IM6.dat.gz', sep = ',')
# MSK.Panel_Y_df = MSK.Panel_df[which(MSK.Panel_df$Chromosome == 'Y'), ]
# write.table(MSK.Panel_Y_df, file = 'Data/Y_chromosome_IMPACT_Panel.example.txt', sep = '\t', row.names = F)

Pileup.Y_df = read.csv('Data/Y_chromosome_Facets_Pileup.txt', sep = '\t')
WES.Y_df = read.csv('Data/Y_chromosome_IMPACT_WES.example.txt', sep = '\t')
Panel.Y_df = read.csv('Data/Y_chromosome_IMPACT_Panel.example.txt', sep = '\t')
Coding.genes_df = read.csv(file = '~/Documents/MSKCC/00_Data/Gene.Model.Biomart.txt', sep = '\t')

## Functions:
#' split sequence in n parts
chunk.function = function(x, n){
  split(x, cut(seq_along(x), n, labels = F))
}

#' min and max position per bin
min.function = function(data, bin){min(data$loc[which(data$bin == bin)])}
max.function = function(data, bin){max(data$loc[which(data$bin == bin)])}

#' select whether Tumor OR Normal Coverage per bin should be calculated
ave.function = function(data, type, bin){
  message('type = mean, median, max, min')
  message('TumorCoverage is used as default')
  if(type == 'mean'){
    mean(data$Tumor_Coverage[which(data$bin == bin)])
    } else if (type == 'median') {
    median(data$Tumor_Coverage[which(data$bin == bin)])
  } else {
    max(data$Tumor_Coverage)[which(data$bin == bin)]
  }
}

#' select probes which overcome default Facets parameter; see github/facets
Facets.preprocess = function(data){
  # default input parameters
  err.thresh = Inf 
  del.thresh = Inf
  ndepth = 35
  ndepthmax = 1000
  
  data$Normal_Coverage = data$File1A + data$File1R
  data$Tumor_Coverage = data$File2A + data$File2R
  
  # discard non-fintite values
  ii = which(data$File1E <= err.thresh & 
               data$File1D <= del.thresh & 
               data$File2E <= err.thresh & 
               data$File2D <= del.thresh)
  
  data = data[ii, ]
  
  # discard low coverage positions
  depthN.keep = data$Normal_Coverage >= ndepth & 
    data$Normal_Coverage < ndepthmax
  
  data = data[depthN.keep, ]
}  

#' change color saturation; alpha
t_col = function(color, percent = 50, name = NULL) {
  rgb.val = col2rgb(color)
  t.col = rgb(rgb.val[1], rgb.val[2], rgb.val[3],
              max = 255,
              alpha = (100 - percent) * 255 / 100,
              names = name)
  
  return(t.col)
}

## Analysis: 
#' create Y_chromosome Master Structure File
bin.number = 1000  # one bin ~ 59 kbp long (59,000 bps)
start.chromosome = 10000
end.chromosome = 59000000

Y.chromo.MasterFile_df = data.frame(bin = NA,
                                    loc = seq(start.chromosome, end.chromosome, 1))
                        
bins.out = lapply(seq(1, bin.number, 1), 
                  function(x) rep(paste0('bin:', x), 
                                  nrow(Y.chromo.MasterFile_df) / bin.number))

bins.out = unlist(bins.out)
bins.out = c(bins.out, bins.out[length(bins.out)])
Y.chromo.MasterFile_df$bin = bins.out


#' merge Y_chromosome MasterFile with Panel Count matrix (Panel 468 in this example)
MSK.Panel_Y_df = Facets.preprocess(data = MSK.Panel_Y_df)
Y.panel.merged_df = merge(Y.chromo.MasterFile_df, 
                          MSK.Panel_Y_df[,c('Position', 'Normal_Coverage', 'Tumor_Coverage')], 
                          by.x = 'loc', 
                          by.y = 'Position', 
                          all.x = T)

#' replace NA in those bins without coverage
Y.panel.merged_df[is.na(Y.panel.merged_df)] = 0


#' bin summary for IMPACT-Panel example
Panel.summary.out = data.frame(bin = unique(Y.panel.merged_df$bin),
                               start = unlist(lapply(unique(Y.panel.merged_df$bin), function(x) min.function(data = Y.panel.merged_df, bin = x))),
                               end = unlist(lapply(unique(Y.panel.merged_df$bin), function(x) max.function(data = Y.panel.merged_df, bin = x))),
                               ave.coverage = unlist(lapply(unique(Y.panel.merged_df$bin), function(x) ave.function(data = Y.panel.merged_df, type = 'mean', bin = x))))

Panel.summary.out$sort = extract_numeric(Panel.summary.out$bin)
Panel.summary.out$sort = factor(Panel.summary.out$sort, levels = order(Panel.summary.out$sort))


#' bin summary for IMPACT WES example:
WES.Y_df = Facets.preprocess(data = WES.Y_df)

## merge with MasterFile:
Y.WES.merged_df = merge(Y.chromo.MasterFile_df, 
                        WES.Y_df[,c('Position', 'Normal_Coverage', 'Tumor_Coverage')], 
                        by.x = 'loc', 
                        by.y = 'Position', 
                        all.x = T)

Y.WES.merged_df[is.na(Y.WES.merged_df)] = 0

#' bin-wide summary stats for WES
WES.summary_out = data.frame(bin = unique(Y.WES.merged_df$bin),
                             start = unlist(lapply(unique(Y.WES.merged_df$bin), function(x) min.function(data = Y.WES.merged_df, bin = x))),
                             end = unlist(lapply(unique(Y.WES.merged_df$bin), function(x) max.function(data = Y.WES.merged_df, bin = x))),
                             ave.coverage = unlist(lapply(unique(Y.WES.merged_df$bin), function(x) ave.function(data = Y.WES.merged_df, bin = x, type = 'mean'))))


WES.summary_out$sort = extract_numeric(WES.summary_out$bin)
WES.summary_out$sort = factor(WES.summary_out$sort, levels = order(WES.summary_out$sort))


#' genomic annotations of the Y chromosome
Y.genes_df = Coding.genes_df[which(Coding.genes_df$Chromosome.scaffold.name == 'Y'), ]
Y.genes_df = Y.genes_df[!duplicated(Y.genes_df$Gene.stable.ID.version), ]
Y.genes_df = Y.genes_df[, c('Gene.name', 'Gene.type', 'Gene.start..bp.', 'Gene.end..bp.')]
colnames(Y.genes_df) = c('Gene', 'Gene.type', 'Start', 'End')

#' fetch bin where Y-chromosome gene(s) map
Y.start = min(Y.chromo.MasterFile_df$loc)
Y.end = max(Y.chromo.MasterFile_df$loc)

Y.chromo.sequence = unlist(apply(data.frame(start = Y.start, end = Y.end), 1, Reduce, f = seq))

Y.gene.coverage_df = data.frame()
for(i in 1:nrow(Y.genes_df)){
  print(i)
  if(apply(Y.genes_df[i, c(3, 4)], 1, function(x) any(seq(x[1], x[2]) %in% Y.chromo.sequence))){
    gene.out = data.frame(bin = Y.chromo.MasterFile_df$bin[which(Y.chromo.MasterFile_df$loc == Y.genes_df[i, 'Start'])],
                          gene = Y.genes_df[i, 'Gene'])
  } else next
  Y.gene.coverage_df = rbind(Y.gene.coverage_df, gene.out)
}

row.names(Y.gene.coverage_df) = NULL


#' if more genes per bin; select protein coding
Y.chromo.genes_df = merge(Y.gene.coverage_df, 
                          Y.genes_df, 
                          by.x = 'gene', 
                          by.y = 'Gene',
                          all.x = T)
Y.chromo.genes_df = Y.chromo.genes_df[!duplicated(Y.chromo.genes_df$gene), ]

bins.plot = data.frame()
for(i in unique(Y.chromo.genes_df$bin)){
  sub.bin = Y.chromo.genes_df[which(Y.chromo.genes_df$bin == i), ]
  
  if(nrow(sub.bin) == 1){
    bin.out = sub.bin
    
  } else if(nrow(sub.bin) > 1 & 'protein_coding' %in% sub.bin$Gene.type){
    bin.out = sub.bin[which(sub.bin$Gene.type == 'protein_coding'), ]
    bin.out = bin.out[1, ]
    
  } else if (nrow(sub.bin) > 1 & !is.element('protein_coding', sub.bin$Gene.type)) {
    bin.out = sub.bin[1, ]
  }
  
  bins.plot = rbind(bins.plot, bin.out)
}

#' backfill every bin (1 through 1000)
missing.bins = setdiff(unique(Y.chromo.MasterFile_df$bin), bins.plot$bin)
missing.bins = data.frame(gene = NA,
                          Gene.type = NA,
                          bin = missing.bins)

Feature.Y.plot_df = rbind(bins.plot[,c(1,3,2)], missing.bins)
Feature.Y.plot_df$sort = extract_numeric(Feature.Y.plot_df$bin)
Feature.Y.plot_df$category = ifelse(Feature.Y.plot_df$Gene.type == 'protein_coding', 'protein_coding',
                                    ifelse(Feature.Y.plot_df$Gene.type %in% c("processed_pseudogene", "unprocessed_pseudogene", 
                                                              "transcribed_unprocessed_pseudogene", "transcribed_processed_pseudogene", 
                                                              'rRNA_pseudogene'), 'Pseudogene',
                                           ifelse(Feature.Y.plot_df$Gene.type %in% c('misc_RNA', 'lncRNA', 'snoRNA', 'snRNA'), 'nonCodingRNA',
                                                  ifelse(is.na(Feature.Y.plot_df$Gene.type), '', NA))))

Feature.Y.plot_df$sort = factor(Feature.Y.plot_df$sort, levels = seq(1, bin.number, 1))
rm(missing.bins)

#' ubiquitiously expressed genes in every tissue
#' see own script for GTEX analysis (heatmap) and https://www.karger.com/Article/FullText/508564
ubiquit.genes = c("RPS4Y1",
                  "DDX3Y",
                  "KDM5D",
                  "LINC00278",
                  "EIF1AY",
                  "USP9Y",
                  "PRKY",
                  "UTY",
                  "ZFY",
                  "TTTY14",
                  "TMSB4Y",
                  "NLGN4Y",
                  "USP9Y",
                  "RPS4Y2")

Ubiquit.genes_df = Y.chromo.genes_df[which(Y.chromo.genes_df$gene %in% ubiquit.genes),, drop = F]
missing.bins = setdiff(unique(Y.chromo.MasterFile_df$bin), Ubiquit.genes_df$bin)
missing.bins = data.frame(gene = NA,
                          Gene.type = NA,
                          bin = missing.bins)

Ubi.genes.Y.plot_df = rbind(Ubiquit.genes_df[,c(1,3,2)], missing.bins)
Ubi.genes.Y.plot_df$sort = extract_numeric(Ubi.genes.Y.plot_df$bin)
Ubi.genes.Y.plot_df$sort = factor(Ubi.genes.Y.plot_df$sort, levels = seq(1, bin.number, 1))
rm(missing.bins)

#' Y-specific genes; no homologous on X
Y_specific = c('KDM5D',
               'DDX3Y',
               'SRY',
               'XKR1',
               'CDY2B',
               'CDY2A',
               'XKRY2',
               'HSFY1',
               'HSFY2',
               'PRORY',
               'PRY2',
               'PRY',
               'BPY2',
               'DAZ1',
               'DAZ2',
               'BPY2B',
               'DAZ3',
               'DAZ4',
               'BPY2C')

Y_specific_df = Y.chromo.genes_df[which(Y.chromo.genes_df$gene %in% Y_specific & 
                                          Y.chromo.genes_df$Gene.type == 'protein_coding'), ]

missing.bins = setdiff(unique(Y.chromo.MasterFile_df$bin), Y_specific_df$bin)
missing.bins = data.frame(gene = NA,
                          Gene.type = NA,
                          bin = missing.bins)

Y_specific.plot_df = rbind(Y_specific_df[, c(1,3,2)], missing.bins)
Y_specific.plot_df$sort = extract_numeric(Y_specific.plot_df$bin)
Y_specific.plot_df$sort = factor(Y_specific.plot_df$sort, levels = seq(1, bin.number, 1))


#' Seq. coverage of specific platforms and groups: ALL features (n = 521): BIOMART
Y.chromo.genes_df # contain all features found on Y chromosome and respective bins:
WES.summary_out # bin-wise summary of average Tumor coverage (from sequencing)

Y.chromo.genes_df$cov.WES = NA
for(i in 1:nrow(Y.chromo.genes_df)){
  Y.chromo.genes_df$cov.WES[i] = ifelse(WES.summary_out$ave.coverage[which(WES.summary_out$bin == Y.chromo.genes_df$bin[i])] != 0, 1, 0)
}

#' Seq. coverage of PANEL sequencing: ALL features (n = 521): BIOMART
Y.chromo.genes_df$cov.Panel = NA
for(i in 1:nrow(Y.chromo.genes_df)){
  Y.chromo.genes_df$cov.Panel[i] = ifelse(Panel.summary.out$ave.coverage[which(Panel.summary.out$bin == Y.chromo.genes_df$bin[i])] != 0, 1, 0)
}


#' Seq. coverage of ubiquitously expressed genes
Ubiquit.genes_df$cov.WES = NA
for(i in 1:nrow(Ubiquit.genes_df)){
  Ubiquit.genes_df$cov.WES[i] = ifelse(WES.summary_out$ave.coverage[which(WES.summary_out$bin == Ubiquit.genes_df$bin[i])] != 0, 1, 0)
}

Ubiquit.genes_df$cov.Panel = NA
for(i in 1:nrow(Ubiquit.genes_df)){
  Ubiquit.genes_df$cov.Panel[i] = ifelse(Panel.summary.out$ave.coverage[which(Panel.summary.out$bin == Ubiquit.genes_df$bin[i])] != 0, 1, 0)
}

#' Seq. coverage Y specific (protein coding genes):
Y_specific.plot_df$cov.WES = NA
Y_specific.plot_df$cov.Panel = NA

for(i in 1:nrow(Y_specific.plot_df)){
  Y_specific.plot_df$cov.WES[i] = ifelse(WES.summary_out$ave.coverage[which(WES.summary_out$bin == Y_specific.plot_df$bin[i])] != 0, 1, 0)
}

for(i in 1:nrow(Y_specific.plot_df)){
  Y_specific.plot_df$cov.Panel[i] = ifelse(Panel.summary.out$ave.coverage[which(Panel.summary.out$bin == Y_specific.plot_df$bin[i])] != 0, 1, 0)
}

Y_specific.plot_df = Y_specific.plot_df[!is.na(Y_specific.plot_df$gene), ]

#' prepare data.frame for plotting:
Y.feature.coverage_df = data.frame(Platform = c(rep('WES', 2), rep('Panel', 2)),
                                   Feature = c(rep(c('All', 'Y.specifc'), 2)),
                                   value = c(table(Y.chromo.genes_df$cov.WES)[['1']] / sum(table(Y.chromo.genes_df$cov.WES)[['1']], table(Y.chromo.genes_df$cov.WES)[['0']]),
                                             table(Y_specific.plot_df$cov.WES)[['1']] / sum(table(Y_specific.plot_df$cov.WES)[['1']], table(Y_specific.plot_df$cov.WES)[['0']]),
                                             table(Y.chromo.genes_df$cov.Panel)[['1']] / sum(table(Y.chromo.genes_df$cov.Panel)[['1']], table(Y.chromo.genes_df$cov.Panel)[['0']]),
                                             table(Y_specific.plot_df$cov.Panel)[['1']] / sum(table(Y_specific.plot_df$cov.Panel)[['1']], table(Y_specific.plot_df$cov.Panel)[['0']])))

Y.feature.coverage_df

## Plot analysis results:
#' the general structure of the Y chromosome (scaled to bin size)
#' average (tumor) read coverage from PANEL sequencing (468)
#' average (tumor) read coverage from WES sequencing
#' protein coding (highlighting and the other groups in lighter color)
#' genes ubiquitinously expressed in every gene (GTEX)


## Y chromosome structure:
#' Y chromosome, split into n = bin.numer (input above)
full.y.chromo = ggplot() + 
  scale_x_continuous(limits = c(0, bin.number), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 1.5), expand = c(0, 0)) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'cm'))

## add PAR1 region: 
#' PAR1 start = 10,001 == bin:1
#' PAR1 end = 2649520 == Y.chromo.Masterfile_df$bin[which(Y.chromo.MasterFile_df$loc == 2649520)]
PAR1.start = 1
PAR1.end = unique(Y.chromo.MasterFile_df$bin[which(Y.chromo.MasterFile_df$loc == 2649520)])
PAR1.end = as.numeric(extract_numeric(PAR1.end))

full.y.chromo = full.y.chromo + statebins:::geom_rrect(mapping = aes(xmin = PAR1.start,
                                                                     xmax = PAR1.end,
                                                                     ymin = -1,
                                                                     ymax = 1),
                                                       color = NA,
                                                       fill = '#527ba4',
                                                       radius = grid:::unit(10, 'pt')) +
  annotate('text', x = (PAR1.start + PAR1.end) / 2, 
           y = 1.0, 
           label = 'PAR1',
           hjust = -0.5,
           vjust = 0.5,
           size = 5,
           angle = 90)

## add MSY region:
MSY.start = PAR1.end
MSY.end = unique(Y.chromo.MasterFile_df$bin[which(Y.chromo.MasterFile_df$loc == 59034049)])
MSY.end = bin.number - 1

full.y.chromo = full.y.chromo + 
  statebins:::geom_rrect(mapping = aes(xmin = MSY.start,
                                       xmax = MSY.end,
                                       ymin = -1,
                                       ymax = 1),
                         color = NA,
                         fill = '#e0a347',
                         radius = grid:::unit(10, 'pt')) +
  
  annotate('text', x = (MSY.start + MSY.end) / 2, 
           y = 1.0, 
           label = 'MSY',
           hjust = 0.5,
           vjust = -1,
           size = 5)

## add PAR2 region:
PAR2.start = MSY.end
PAR2.end = bin.number

full.y.chromo = full.y.chromo + statebins:::geom_rrect(mapping = aes(xmin = PAR2.start,
                                                                     xmax = PAR2.end,
                                                                     ymin = -1,
                                                                     ymax = 1),
                                                       color = NA,
                                                       fill = '#527ba4',
                                                       radius = grid:::unit(2, 'pt')) +
  annotate('text', x = (PAR2.start + PAR2.end) / 2, 
           y = 1.0, 
           label = 'PAR2',
           hjust = -0.5,
           vjust = 0,
           size = 5,
           angle = 90)


#' Panel coverage
panel.y.plot = ggplot(Panel.summary.out, 
                      aes(x = sort, 
                          y = ave.coverage, 
                          group = 1)) + 
  geom_line(size = 0.7) +
  theme(axis.text = element_blank(),
        panel.border = element_rect(size = 1.2, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'cm')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = '', y = 'ave. sequencing coverage [a.u.]') +
  annotate('text',
           x = 5,
           y = max(Panel.summary.out$ave.coverage) * 19/20,
           label = 'IMPACT-Panel',
           hjust = 0,
           size = 8)

#' WES coverage
wes.y.plot = ggplot(WES.summary_out, 
                    aes(x = sort, 
                        y = ave.coverage, 
                        group = 1)) + 
  scale_x_discrete(breaks = c(1, 200, 400, 600, 800, 1000),
                   labels = c(paste0(round(WES.summary_out$start[which(WES.summary_out$sort == 1)] / 1e6), ' Mb'),
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == 200)] / 1e6), ' Mb'),
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == 400)] / 1e6), ' Mb'),
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == 600)] / 1e6), ' Mb'),
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == 800)] / 1e6), ' Mb'),
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == 1000)] / 1e6), ' Mb'))) +
  geom_line(size = 0.7) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, face = 'bold', color = 'black'),
        panel.border = element_rect(size = 1.2, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'cm')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = '', y = 'ave. sequencing coverage [a.u.]') +
  annotate('text',
           x = 5, 
           y = max(WES.summary_out$ave.coverage) * 19/20,
           label = 'IMPACT-WES',
           hjust = 0,
           size = 8)
wes.y.plot


#' Y chromosome Feature annotation;
#' Features derived from Biomart: only one feature/bin is shown
c.nonCodingRNA = t_col('#005155', percent = 60)
c.Pseudogene = t_col('#06C1CF', percent = 60)
c.proteinCoding = '#8F1918'

# currently legend is disabled:  
gene_coverage.plot = ggplot(Feature.Y.plot_df, 
                            aes(x = sort, 
                                y = NA, 
                                fill = category)) + 
  geom_tile() +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0,0,0,0), 'cm'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, size = 1.2),
        panel.background = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank()) +
  
  scale_fill_manual(values = c('NA' = 'white',
                      'nonCodingRNA' = c.nonCodingRNA,
                      'Pseudogene' = c.Pseudogene,
                      'protein_coding' = c.proteinCoding),
                    name = 'Type',
                    labels = c('nonCodingRNA', 'Protein Coding', 'Pseudogene', '')) +
  labs(x = '', y = '')

gene_coverage.plot


#' annotate genes which are ubiquitiniously expressed in every tissues (GTEX)
ubi.genes.y.plot = ggplot(Ubi.genes.Y.plot_df,
                          aes(x = sort, 
                              y = NA,
                              fill = Gene.type)) + 
  geom_tile() +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0,0,0,0), 'cm'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, size = 1.2),
        panel.background = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank()) +
  scale_fill_manual(values = c('NA' = 'white',
                               'protein_coding' = c.proteinCoding,
                               'lncRNA' = c.nonCodingRNA),
                    name = 'Type',
                    labels = c('long non-coding RNA', 'Protein Coding', '')) +
  labs(x = '', y = '')

ubi.genes.y.plot

#' annotate Y-specific genes which are protein coding:
breaks.x = Y_specific.plot_df$sort[!is.na(Y_specific.plot_df$gene)]
labels.x = Y_specific.plot_df$gene[match(Y_specific.plot_df$sort, breaks.x)]
labels.x = labels.x[!is.na(labels.x)]

Y_specific_y.plot = ggplot(Y_specific.plot_df,
                          aes(x = sort, 
                              y = NA, fill = Gene.type)) + 
  geom_tile() +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0, 0),
                   breaks = breaks.x,
                   labels = labels.x,
                   guide = guide_axis(check.overlap = T)) +
  
  theme(plot.margin = unit(c(0,0,0,0), 'cm'),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, color = 'black', size = 10),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, size = 1.2),
        panel.background = element_blank(),
        legend.position = 'bottom',
        panel.grid = element_blank()) +
  
  scale_fill_manual(values = c('NA' = 'white',
                               'protein_coding' = c.proteinCoding),
                    name = 'Y-specific') +
  labs(x = '', y = '')


Y_specific_y.plot


#' Y feature coverage plot, across platforms:
c.WES = '#448C82'
c.Panel = '#C35E34'
Y.feature.coverage_df$Platform = factor(Y.feature.coverage_df$Platform, levels = c('WES', 'Panel'))
Y.feature.coverage_df$Feature = ifelse(Y.feature.coverage_df$Feature == 'All', paste0('func.-/non-func.\nelements'),
                                       paste0('Protein Coding\nY_specific'))

Y.coverage_plot = ggplot(Y.feature.coverage_df, 
                         aes(x = Platform,
                             y = value * 100,
                             fill = Platform,
                             alpha = Feature)) +
  
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.9) +
  
  scale_alpha_manual(values = c('Protein Coding\nY_specific' = 1,
                                 'func.-/non-func.\nelements' = 0.5)) +
  
  scale_fill_manual(values = c('WES' = c.WES,
                               'Panel' = c.Panel)) +
  
  geom_text(aes(label = Feature, y = 0.5), 
            position = position_dodge(width = 0.9), 
            color = 'white', 
            angle = 90,
            vjust = 0.5,
            hjust = 0,
            size = 5) +
  
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.5, 0)) +
  
  theme(panel.background = element_blank(),
        legend.position = 'none',
        aspect.ratio = 2,
        axis.line.y = element_line(color = 'black', size = 0.6),
        axis.line.x = element_line(color = 'black', size = 0.6),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10, color = 'black')) +
  labs(x = '', y = 'Percentage elements covered by sequencing')
  

Y.coverage_plot

ggsave(filename = 'Features.covered.sequencing.pdf',
       plot = Y.coverage_plot,
       device = 'pdf',
       width = 4,
       height = 8,
       units = 'in')


structure_grid_plot = plot_grid(full.y.chromo, 
                                panel.y.plot, 
                                wes.y.plot,
                                gene_coverage.plot,
                                ubi.genes.y.plot,
                                Y_specific_y.plot,
                                nrow = 6, 
                                align = 'v',
                                rel_heights = c(0.5, 3, 3, 1, 1, 2))

structure_grid_plot

ggsave(filename = 'Y_chromosome_structure.pdf',
       plot = structure_grid_plot,
       device = 'pdf',
       width = 24,
       height = 12,
       units = 'in')


