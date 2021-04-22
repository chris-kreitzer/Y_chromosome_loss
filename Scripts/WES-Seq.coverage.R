## PanCancer analysis of Y chromosome loss:
## Y_chromosome structure and FacetsY validation:
##
## WES provide us with a granular view on the Y-chromosome
## Is the utilization of FacetsY justified
## 
## The initial analysis started in fall 2020; 
## Initial scripts and pre-analysis were run by Bastien and Subhi (back in 2019 and 2020)
## 
## Here, I will revise some concepts and provide technical validation of Y_chromosome_loss calls
## 
## Start (revision): 11/03/2021:
## Start (updates): 23/03/2021:
## updates: 24/03/2021:
## updates: 14/04/2021:
## 
## Scripts will be hosted on Github: https://github.com/chris-kreitzer/Y_chromosome_loss.git
## Scripts will be hosted equally on GitHub enterprise: https://github.mskcc.org/kreitzeC/Y-chromosome-loss

rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')
source('Y_chromosome_loss/Plotting_theme.R')

## Libraries
library(pctGCdata)
library(devtools)
if('facetsY' %in% rownames(installed.packages())){
  message('FacetsY is already installed on machine')
  library(facetsY)
} else{
  message('FacetsY will be installed from GitHub')
  devtools::install_github('https://github.com/BastienNguyen/facetsY.git', 
                           build_vignettes = F, force = T)
}
library(tidyverse)
library(cowplot)
library(statebins)
library(tempR)
library(dplyr)
library(patchwork)
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(ggplot2))


## Input:
WES.coverage_out = read.csv('Data_out/tmp.coverageWES.txt', sep = '\t') #' Seq-Coverage of WES (n = 3)
WES.breakpoints_df = read.csv('Data_out/WES.Breakpoints.Y.txt', sep = '\t')
WES.features.covered_df = read.csv('Data_out/WES.Features.covered.out.txt', sep = '\t')
Y.chromosome.genes_df = read.csv('Data_out/Y_chromosome_genes.txt', sep = '\t')

#' retrieve hg19 genes (as FacetsY was run on hg19)
# ensembl = useEnsembl(biomart = 'ensembl', 
#                      dataset = 'hsapiens_gene_ensembl', 
#                      GRCh = 37)
# Y_genes = getBM(attributes = c('ensembl_gene_id', 
#                                'hgnc_symbol', 
#                                'chromosome_name', 
#                                'start_position', 
#                                'end_position',
#                                'external_gene_name',
#                                'gene_biotype'), 
#                 filters = 'chromosome_name', 
#                 values = "Y", 
#                 mart = ensembl)

## Functions
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

#' count how many positions (bases) or covered > threshold in given bin
coverage.count.function = function(data, bin){
  n = length(data$Tumor_Coverage[which(data$Tumor_Coverage != 0 & data$bin == bin)])
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



## Processing
#' summary for Seq-Coverage across the Y-chromosome
#' here I have just one sample run:
WES.coverage_out = WES.coverage_out[order(WES.coverage_out$sort), ]


#' look into how many feature are covered by WESeq
number.protein.coding.genesY = length(Y.chromosome.genes_df$external_gene_name[which(Y.chromosome.genes_df$gene_biotype == 'protein_coding')])

WES.coverage.summary_df = data.frame()
for(i in unique(WES.features.covered_df$ID)){
  rel.amount = data.frame(ID = i,
                          rel.covered = sum(WES.features.covered_df$Match[which(WES.features.covered_df$ID == i)]) / number.protein.coding.genesY)
  WES.coverage.summary_df = rbind(WES.coverage.summary_df, rel.amount)
}

#' Features covered by WES;
WES.covered = Y.chromosome.genes_df[which(Y.chromosome.genes_df$external_gene_name %in% unique(WES.features.covered_df$Gene)), ]


## Analysis
#' Y-chromosome coverage:
bin.number = 5900  # one bin 10 kbp long (50,000 bps)
start.chromosome = 1
end.chromosome = 59000000
Y.chromo.MasterFile_df = data.frame(bin = NA,
                                    loc = seq(start.chromosome, end.chromosome, 1))
bins.out = lapply(seq(1, bin.number, 1), 
                  function(x) rep(paste0('bin:', x), 
                                  nrow(Y.chromo.MasterFile_df) / bin.number))
bins.out = unlist(bins.out)
Y.chromo.MasterFile_df$bin = bins.out

##-----------------------------------------------------------------------------
#' merge with MasterFile:
#' Y.WES.merged_df = merge(Y.chromo.MasterFile_df, 
#'                         WES.Y_df[, c('Position', 'Normal_Coverage', 'Tumor_Coverage')], 
#'                         by.x = 'loc', 
#'                         by.y = 'Position', 
#'                         all.x = T)
#' 
#' Y.WES.merged_df[is.na(Y.WES.merged_df)] = 0
#' 
#' #' bin-wide summary stats for WES
#' WES.summary_out = data.frame(bin = unique(Y.WES.merged_df$bin),
#'                              start = unlist(lapply(unique(Y.WES.merged_df$bin), function(x) min.function(data = Y.WES.merged_df, bin = x))),
#'                              end = unlist(lapply(unique(Y.WES.merged_df$bin), function(x) max.function(data = Y.WES.merged_df, bin = x))),
#'                              SNPs = unlist(lapply(unique(Y.WES.merged_df$bin), function(x) coverage.count.function(data = Y.WES.merged_df, bin = x))))
#' 
#' 
#' WES.summary_out$sort = extract_numeric(WES.summary_out$bin)
#' WES.summary_out$sort = factor(WES.summary_out$sort, levels = order(WES.summary_out$sort))
# write.table(WES.summary_out, file = 'Data_out/tmp.coverageWES.txt', sep = '\t', row.names = F)
##-----------------------------------------------------------------------------

#' genomic annotation of Y-chromosome features
#' here I use the ENSEMBL coordinates of hg19 
WES.covered = WES.covered[, c('hgnc_symbol', 'gene_biotype', 'start_position', 'end_position')]
colnames(WES.covered) = c('Gene', 'Gene.type', 'Start', 'End')

#' assign bins to respective genes (mapping to chromosome)
Y.chromo.start = min(Y.chromo.MasterFile_df$loc)
Y.chromo.end = max(Y.chromo.MasterFile_df$loc)

Y.chromo.sequence = unlist(apply(data.frame(start = Y.chromo.start, end = Y.chromo.end), 1, Reduce, f = seq))

WES.gene.coverage_df = data.frame()
for(i in 1:nrow(WES.covered)){
  print(i)
  if(apply(WES.covered[i, c(3, 4)], 1, 
           function(x) any(seq(x[1], x[2]) %in% Y.chromo.sequence))){
    gene.out = data.frame(bin = Y.chromo.MasterFile_df$bin[which(Y.chromo.MasterFile_df$loc == WES.covered[i, 'Start'])],
                          gene = WES.covered[i, 'Gene'])
  } else next
  WES.gene.coverage_df = rbind(WES.gene.coverage_df, gene.out)
}

row.names(WES.gene.coverage_df) = NULL

#' backfill with missing bins
missing.bins = data.frame(gene = NA,
                          bin = setdiff(unique(Y.chromo.MasterFile_df$bin), WES.gene.coverage_df$bin))
WES.Feature.Y.plot_df = rbind(WES.gene.coverage_df[, c(1,2)], missing.bins)
WES.Feature.Y.plot_df$sort = extract_numeric(WES.Feature.Y.plot_df$bin)
WES.Feature.Y.plot_df = WES.Feature.Y.plot_df[!duplicated(WES.Feature.Y.plot_df$bin), ]
WES.Feature.Y.plot_df = WES.Feature.Y.plot_df[order(WES.Feature.Y.plot_df$sort), ]



## Visualization:
#' Seq Coverage across the Y-chromosome
breaks.max = max(WES.coverage_out$sort)
breaks.seq = seq(0, breaks.max, length.out = 6)
breaks.seq[1] = 1
y.max = max(WES.coverage_out$SNPs) + 1

WES.Y.coverage_plot = ggplot(WES.coverage_out, 
                                aes(x = sort, 
                                    y = SNPs, 
                                    group = 1)) + 
  
  geom_line(size = 0.7) +
  
  annotate('segment',
           x = 5000, xend = 5500,
           y = 55, yend = 55,
           size = 1.7) +
  
  annotate('text',
           x = 5000, 
           y = 55,
           label = '≈10 Markers',
           hjust = 0,
           vjust = -1) +
  
  scale_x_continuous(expand = c(0.01, 0),
                     breaks = breaks.seq,
                     labels = c('0',
                                paste0('~', round(IMPACT.coverage.summary$start[which(IMPACT.coverage.summary$sort == breaks.seq[[2]])] / 1e6), ' Mb'),
                                paste0('~', round(IMPACT.coverage.summary$start[which(IMPACT.coverage.summary$sort == breaks.seq[[3]])] / 1e6), ' Mb'),
                                paste0('~', round(IMPACT.coverage.summary$start[which(IMPACT.coverage.summary$sort == breaks.seq[[4]])] / 1e6), ' Mb'),
                                paste0('~', round(IMPACT.coverage.summary$start[which(IMPACT.coverage.summary$sort == breaks.seq[[5]])] / 1e6), ' Mb'),
                                paste0('~', round(IMPACT.coverage.summary$start[which(IMPACT.coverage.summary$sort == breaks.seq[[6]])] / 1e6), ' Mb'))) +
  
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, y.max)) +

theme_Y(base_size = 16, base_rect_size = 0.1, base_line_size = 0.1) +
  
  theme(axis.text.y = element_blank(),
        panel.border = element_rect(size = 1.0, fill = NA),
        axis.ticks = element_blank(),
        plot.subtitle = element_text(size = 13)) +
  
  labs(title = 'WE-Seq resolution across the Y-chromosome', 
       subtitle = 'Peaks depict the number of markers in given bins (1bin = 1Kb) which retain for segmentation',
       x = '', y = 'Average number of markers')


#' geom_tile for gene annotation along the Y-chromosome
#' annotate genes which are directly covered by FacetsY
breaks.x = WES.Feature.Y.plot_df$sort[!is.na(WES.Feature.Y.plot_df$gene)]
labels.x = c()
for(i in 1:length(breaks.x)){
  label = WES.Feature.Y.plot_df$gene[which(WES.Feature.Y.plot_df$sort == breaks.x[i])]
  labels.x = c(labels.x, label)
}


color.out = c('SRY' = 'red',
'RPS4Y1' = 'red',
'ZFY' = 'red',
'TGIF2LY' = 'red',
'PCDH11Y' = 'red',
'TSPY2' = 'red',
'AMELY' = 'red',
'TBL1Y' = 'red',
'TSPY4' = 'red',
'TSPY8' = 'red',
'TSPY3' = 'red',
'TSPY6P' = 'red',
'TSPY10' = 'red',
'FAM197Y1' = 'red',
'SLC9B1P1' = 'red',
'USP9Y' = 'red',
'DDX3Y' = 'red',
'UTY' = 'red',
'TMSB4Y' = 'red',
'NLGN4Y' = 'red',
'KDM5D' = 'red',
'EIF1AY' = 'red',
'RPS4Y2' = 'red',
'CYorf17' = 'red',
'RBMY1A1' = 'red',
'RBMY1E' = 'red',
'RBMY1F' = 'red',
'RBMY1J' = 'red',
'DAZ1' = 'red',
'DAZ2' = 'red',
'DAZ3' = 'red',
'DAZ4' = 'red',
'BPY2C' = 'red',
'MISSING' = 'red')

WES.coverage.plot = ggplot(WES.Feature.Y.plot_df , 
                              aes(x = sort,
                                  y = NA,
                                  fill = gene)) + 
  geom_tile(na.rm = T) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0.01, 0),
                   breaks = breaks.x,
                   labels = labels.x,
                   guide = guide_axis(check.overlap = T)) +
  scale_fill_manual(values = color.out) +
  theme_Y(base_size = 16) +
  theme(plot.margin = unit(c(0,0,0,0), 'cm'),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1,
                                   color = 'black', 
                                   size = 10),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA, size = 0.8),
        panel.background = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = 'Genomic position', y = '')
  
  

Fig1A = WES.Y.coverage_plot / WES.coverage.plot +
  plot_layout(widths = unit(c(25, 25), c('cm', 'cm')),
              heights = unit(c(10, 1), c('cm', 'cm'))) 

  
ggsave_golden(filename = 'Figures/WES_Seq_coverage_full.pdf', plot = Fig1A, width = 16)
