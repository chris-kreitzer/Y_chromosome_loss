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
## updates: 24/03/2021
## 
## Scripts will be hosted on Github: https://github.com/chris-kreitzer/Y_chromosome_loss.git
## Scripts will be hosted equally on GitHub enterprise: https://github.mskcc.org/kreitzeC/Y-chromosome-loss

rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')


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


## Input data:

WES.Y_df = read.csv('Data/Y_chromosome_IMPACT_WES.example.txt', sep = '\t')
Panel.Y_df = read.csv('Data/Y_chromosome_IMPACT_Panel.example.txt', sep = '\t')
WES.breakpoints_df = read.csv('Data_out/breakpoints_WES.txt', sep = '\t')
Panel.coverage.3samples_df = read.csv('Data_out/Panel.coverage.txt', sep = '\t')
Ubi.expressed.genes = read.csv('../Data_out/ubiquitinously.expressed.genes.txt', sep = '\t')
WES.summary_out = read.csv('../Data_out/tmp.coverageWES.txt', sep = '\t')

#' retrieve hg19 genes (as FacetsY was run on hg19)
ensembl = useEnsembl(biomart = 'ensembl', 
                     dataset = 'hsapiens_gene_ensembl', 
                     GRCh = 37)
Y_genes = getBM(attributes = c('ensembl_gene_id', 
                               'hgnc_symbol', 
                               'chromosome_name', 
                               'start_position', 
                               'end_position',
                               'external_gene_name',
                               'gene_biotype'), 
                filters = 'chromosome_name', 
                values = "Y", 
                mart = ensembl)


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

#' apply quality metrics (default by Facets when analysis)
#' Panel sequencing example
MSK.Panel_Y_df = Facets.preprocess(data = Panel.Y_df)
#' WES example
WES.Y_df = Facets.preprocess(data = WES.Y_df)


## Analysis

#' quality check: 
#' average positions (breakpoints) kept when using FacetsY on n = 55 samples (PANEL sequencing)
#' this example runs on n = 55 prostate IMPACT panel samples (locally stored)
Panel.sequencing.examples = list.files('~/Desktop/mnt/ATMcountdata/', full.names = T)
breakpoints.panel.out_df = data.frame()
for(i in sample(Panel.sequencing.examples, 55, replace = T)){
  input = vroom::vroom(i)
  input = input[which(input$Chromosome == 'Y'), ]
  input.raw = nrow(input)
  input.modi = nrow(Facets.preprocess(input))
  ave.breaks = data.frame(id = str_extract(i, pattern = 'P-.*'),
                          seq.accessed = input.raw,
                          breakpoints.used_facets = input.modi)
  
  breakpoints.panel.out_df = rbind(breakpoints.panel.out_df, ave.breaks)
}

#' breakpoints from WES seq.
#' The data obtained from this script was run on the xjuno cluster
WES.breakpoints.density = density(WES.breakpoints_df$breakpoints.used_facets, bw = 600) # pre-defined bandwith


#' Y_chromosome master structure file:
#' The Y-chromosome is split into 5,900 bins (1 bin = 10 Kbp)
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


#' Y_chromosome MasterFile will be merged with Panel data (where reads are available) (Panel 468 in this example)
Y.panel.merged_df = merge(Y.chromo.MasterFile_df, 
                          MSK.Panel_Y_df[, c('Position', 'Normal_Coverage', 'Tumor_Coverage')], 
                          by.x = 'loc', 
                          by.y = 'Position',
                          all.x = T)

#' replace NA in those bins without coverage
Y.panel.merged_df[is.na(Y.panel.merged_df)] = 0

#' bin summary for IMPACT-Panel example
Panel.summary.out = data.frame(bin = unique(Y.panel.merged_df$bin),
                               start = unlist(lapply(unique(Y.panel.merged_df$bin), function(x) min.function(data = Y.panel.merged_df, bin = x))),
                               end = unlist(lapply(unique(Y.panel.merged_df$bin), function(x) max.function(data = Y.panel.merged_df, bin = x))),
                               ave.coverage = unlist(lapply(unique(Y.panel.merged_df$bin), function(x) ave.function(data = Y.panel.merged_df, type = 'mean', bin = x))),
                               SNPs.covered = unlist(lapply(unique(Y.panel.merged_df$bin), function(x) coverage.count.function(data = Y.panel.merged_df, bin = x))))

Panel.summary.out$sort = extract_numeric(Panel.summary.out$bin)
Panel.summary.out$sort = factor(Panel.summary.out$sort, levels = order(Panel.summary.out$sort))

#' as the above example exclusively provides us with one example, I repeated this analysis on the juno cluster
#' to fetch the information on 3 independent samples (better overview of Y chromosome coverage)

#' make a bin-wide summary for the 3 independent samples, where the above functions were run:
Panel.Y.coverage.summary = Panel.coverage.3samples_df %>%
  group_by(bin) %>%
  summarize(mean.SNPs = mean(SNPs.covered))

Panel.Y.coverage.summary$sort = extract_numeric(Panel.Y.coverage.summary$bin)
Panel.Y.coverage.summary$sort = factor(Panel.Y.coverage.summary$sort, levels = seq(1, 5900, 1))


#' repeat the same analysis with WES samples
#' merge with MasterFile:
Y.WES.merged_df = merge(Y.chromo.MasterFile_df, 
                        WES.Y_df[, c('Position', 'Normal_Coverage', 'Tumor_Coverage')], 
                        by.x = 'loc', 
                        by.y = 'Position', 
                        all.x = T)

Y.WES.merged_df[is.na(Y.WES.merged_df)] = 0

#' bin-wide summary stats for WES
WES.summary_out = data.frame(bin = unique(Y.WES.merged_df$bin),
                             start = unlist(lapply(unique(Y.WES.merged_df$bin), function(x) min.function(data = Y.WES.merged_df, bin = x))),
                             end = unlist(lapply(unique(Y.WES.merged_df$bin), function(x) max.function(data = Y.WES.merged_df, bin = x))),
                             SNPs = unlist(lapply(unique(Y.WES.merged_df$bin), function(x) coverage.count.function(data = Y.WES.merged_df, bin = x))))


WES.summary_out$sort = extract_numeric(WES.summary_out$bin)
WES.summary_out$sort = factor(WES.summary_out$sort, levels = order(WES.summary_out$sort))
# write.table(WES.summary_out, file = 'Data_out/tmp.coverageWES.txt', sep = '\t', row.names = F)


#' genomic annotation of Y-chromosome features
#' here I use the ENSEMBL coordinates of hg19 
Y_genes = Y_genes[!duplicated(Y_genes$ensembl_gene_id), ]
Y_genes = Y_genes[, c('hgnc_symbol', 'gene_biotype', 'start_position', 'end_position')]
colnames(Y_genes) = c('Gene', 'Gene.type', 'Start', 'End')

#' assign bins to respective genes (mapping to chromosome)
Y.chromo.start = min(Y.chromo.MasterFile_df$loc)
Y.chromo.end = max(Y.chromo.MasterFile_df$loc)

Y.chromo.sequence = unlist(apply(data.frame(start = Y.chromo.start, end = Y.chromo.end), 1, Reduce, f = seq))

Y.gene.coverage_df = data.frame()
for(i in 1:nrow(Y_genes)){
  print(i)
  if(apply(Y_genes[i, c(3, 4)], 1, 
           function(x) any(seq(x[1], x[2]) %in% Y.chromo.sequence))){
    gene.out = data.frame(bin = Y.chromo.MasterFile_df$bin[which(Y.chromo.MasterFile_df$loc == Y_genes[i, 'Start'])],
                          gene = Y_genes[i, 'Gene'])
  } else next
  Y.gene.coverage_df = rbind(Y.gene.coverage_df, gene.out)
}

row.names(Y.gene.coverage_df) = NULL

#' merge genes with gene biotype
Y_genes = Y_genes[which("" != Y_genes$Gene), ]
Y_genes$bin = NA

for(i in 1:nrow(Y_genes)){
  if(Y_genes$Gene[i] %in% Y.gene.coverage_df$gene){
    Y_genes$bin[i] = Y.gene.coverage_df$bin[which(Y.gene.coverage_df$gene == Y_genes$Gene[i])]
  } else next
}


#' display only genes which are protein coding and those
#' features which are ubiquitinously expressed
gene.plot = c(as.character(unique(Ubi.expressed.genes$Features)), 
              as.character(unique(Y_genes$Gene[which(Y_genes$Gene.type == 'protein_coding')])))
gene.plot = unique(gene.plot)

Genes.plot_df = Y_genes[which(Y_genes$Gene %in% gene.plot),, drop = F]
Genes.plot_df = Genes.plot_df[!duplicated(Genes.plot_df$bin), ]

#' backfill every bin (1 through 5900): needed for plot
missing.bins = setdiff(unique(Y.chromo.MasterFile_df$bin),Genes.plot_df$bin)
missing.bins = data.frame(Gene = NA,
                          Gene.type = NA,
                          bin = missing.bins)

Feature.Y.plot_df = rbind(Genes.plot_df[, c(1,2,5)], missing.bins)
Feature.Y.plot_df$sort = extract_numeric(Feature.Y.plot_df$bin)
Feature.Y.plot_df$sort = factor(Feature.Y.plot_df$sort, levels = seq(1, bin.number, 1))
Feature.Y.plot_df$Gene.type = factor(Feature.Y.plot_df$Gene.type, levels = c('antisense', 'lincRNA', 'pseudogene', 'protein_coding'))
Feature.Y.plot_df$sort.new = as.character(as.factor(Feature.Y.plot_df$sort))

#' check how many features are DIRECTLY covered by WES Seq.
Feature.Y.plot_df$cov.WES = NA
for(i in 1:nrow(Feature.Y.plot_df)){
  Feature.Y.plot_df$cov.WES[i] = ifelse(WES.summary_out$SNPs[which(WES.summary_out$bin == Feature.Y.plot_df$bin[i])] != 0, 1, 0)
}


## Visualizations:
#' WES coverage (SNPs along the Y-chromosome):
#' average (tumor) read coverage from WES sequencing

#' Y chromosome coverage of WES seq
breaks.max = bin.number
breaks.seq = seq(0, breaks.max, length.out = 6)
breaks.seq[1] = 1

wes.y.plot = ggplot(WES.summary_out, 
                    aes(x = sort, 
                        y = SNPs, 
                        group = 1)) + 
  
  geom_line(size = 0.7) +
  
  scale_x_discrete(expand = c(0.01, 0),
                   breaks = breaks.seq,
                   labels = c('0',
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == breaks.seq[[2]])] / 1e6), ' Mb'),
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == breaks.seq[[3]])] / 1e6), ' Mb'),
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == breaks.seq[[4]])] / 1e6), ' Mb'),
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == breaks.seq[[5]])] / 1e6), ' Mb'),
                              paste0('~', round(WES.summary_out$start[which(WES.summary_out$sort == breaks.seq[[6]])] / 1e6), ' Mb'))) +
  
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, max(WES.summary_out$SNPs + 1))) +
  
  theme_Y(base_size = 14, base_rect_size = 0.1, base_line_size = 0.1) +
  
  theme(axis.text.y = element_blank(),
        panel.border = element_rect(size = 1.0, fill = NA),
        axis.ticks = element_blank()) +
  
  labs(title = 'IMPACT-WES Seq. resolution on Y-chromosome', 
       subtitle = 'Peaks depict bins (1bin = 1Kb), where informative SNPs retain for segmentation',
       x = '', y = 'Coverage')



#' Breakpoint analysis of 
WES.density = data.frame(x = WES.breakpoints.density$x,
                         y = WES.breakpoints.density$y)

wes.density = ggplot(WES.density,
                     aes(x = x, 
                         y = y)) + 
  
  geom_line(size = 1) +
  
  geom_vline(xintercept = round(WES.density$x[which.max(WES.density$y)]),
             color = 'lightblue', 
             linetype = 'dashed', 
             size = 1) +
  
  scale_x_continuous(breaks = c(round(WES.density$x[which.max(WES.density$y)])),
                     labels = c(round(WES.density$x[which.max(WES.density$y)]))) +
  
  scale_y_continuous(limits = c(0, max(WES.density$y) + max(WES.density$y) / 20)) +
  
  coord_cartesian(expand = F, 
                  clip = 'off') +
  
  annotate('curve',
           x = 1500, 
           xend = 1053,
           y = 0.00007, yend = 0.000005,
           curvature = 0.25,
           size = 0.6,
           arrow = arrow(length = unit(0.1, 'inches'), type = 'closed')) +
  
  annotate('text', 
           x = 1510,
           y = 0.000070,
           label = 'Median',
           family = "RobotoCondensed-Regular",
           size = 3.3, 
           color = "grey55",
           fontface = "bold.italic", 
           lineheight = .85,
           hjust = 0,
           vjust = 0.5) +
  
  theme_Y(base_size = 12) +
  
  theme(axis.text.y = element_blank(),
        panel.border = element_rect(size = 0.8, fill = NA),
        axis.ticks.y = element_blank()) +
  
  
  labs(title = '', y = 'Density', x = '',
       subtitle = '# of base positions (breakpoints) for FacetsY')


#' make the inset (WES coverage + WES density (breakpoints))
WES_inset = wes.y.plot + 
  annotation_custom(grob = ggplotGrob(wes.density),
                    xmin = 2700, xmax = 5800,
                    ymin = 20, 
                    ymax = 62)

ggsave_golden(filename = 'Figures/WES_Seq_Y_coverage.pdf', plot = WES_inset, width = 12)


#' geom_tile for gene annotation along the Y-chromosome
#' annotate genes which are directly covered by FacetsY
breaks.x = Feature.Y.plot_df$sort.new[which(Feature.Y.plot_df$cov.WES == 1 &!is.na(Feature.Y.plot_df$Gene))]
breaks.x = as.character(as.factor(breaks.x))
labels.x = c()
for(i in 1:length(breaks.x)){
  label = Feature.Y.plot_df$Gene[which(Feature.Y.plot_df$sort == breaks.x[i])]
  labels.x = c(labels.x, label)
}

#' make the geom_tile plot
c.proteinCoding = '#E61F2B'
c.nonCodingRNA = '#5EA1C3'
c.Pseudogene = '#006AA0'
c.antisense = '#7FC9C7'

gene_coverage.plot = ggplot(Feature.Y.plot_df, 
                            aes(x = sort, 
                                y = NA, 
                                fill = Gene.type)) + 
  geom_tile(na.rm = T) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0, 0),
                   breaks = breaks.x,
                   labels = labels.x,
                   guide = guide_axis(check.overlap = T)) +
  theme_Y() +
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
        legend.position = 'bottom',
        panel.grid = element_blank()) +
  
  scale_fill_manual(values = c('NA' = 'white',
                               'lincRNA' = c.nonCodingRNA,
                               'pseudogene' = c.Pseudogene,
                               'protein_coding' = c.proteinCoding,
                               'antisense' = c.antisense),
                    name = 'Type') +
  labs(x = 'Genomic position', y = '')

Fig1A = WES_inset / gene_coverage.plot + 
  plot_layout(widths = unit(c(35, 35), c('cm', 'cm')),
              heights = unit(c(10, 1), c('cm', 'cm'))) 

ggsave_golden(filename = 'Figures/WES_Seq_coverage_full.pdf', plot = Fig1A, width = 8)



## Repeat the same analysis for Panel data;
## to be done: 03/28/2021

#' informative (density) of postions used for Facets Clustering
density.SNPs.Panel = density(breakpoints.panel.out_df$breakpoints.used_facets,
                             bw = 40)
panel.density.plot = data.frame(x = density.SNPs.Panel$x,
                                y = density.SNPs.Panel$y)

panel.density = ggplot(panel.density.plot, 
                       aes(x = x, y = y)) + 
  geom_line(size = 1) +
  theme_Y(base_size = 12) +
  coord_cartesian(expand = F, clip = 'off') +
  scale_x_continuous(breaks = c(0, round(panel.density.plot$x[which.max(panel.density.plot$y)]), round(max(panel.density.plot$x))),
                     labels = c(0, round(panel.density.plot$x[which.max(panel.density.plot$y)]), round(max(panel.density.plot$x)))) +
  theme(axis.text.y = element_blank(),
        panel.border = element_rect(size = 0.8, fill = NA),
        axis.ticks.y = element_blank()) +
  geom_vline(xintercept = round(panel.density.plot$x[which.max(panel.density.plot$y)]),
             color = 'lightblue', linetype = 'dashed', size = 1) +
  labs(title = 'IMPACT-Panel Seq.', y = 'Density', x = '# informative SNPs')


#' Y chromosome coverage of IMPACT Panel Seq. (informative Positions)
#' working with Y_coverage_summary: bin summary of 15 independent Panel samples
#' panel.y.plot = ggplot(data = Y_coverage_summary, 
#'                       aes(x = sort, 
#'                           y = mean.SNPs, 
#'                           group = 1)) + 
#'   
#'   geom_line(size = 0.7) + 
#'   
#'   scale_y_continuous(expand = c(0, 0),
#'                      limits = c(0, max(Y_coverage_summary$mean.SNPs) + 1)) +
#'   
#'   theme_Y(base_size = 12, base_line_size = 0.1, base_rect_size = 0.1) +
#'   
#'   theme(panel.border = element_rect(size = 1.0, fill = NA),
#'         axis.text.x = element_blank(),
#'         axis.ticks = element_blank(),
#'         axis.text.y = element_blank()) +
#'   
#'   labs(x = '', y = 'Seq. coverage')
#' 
#' 
#' #' add inset plot to main plot
#' Panel_inset = panel.y.plot + 
#'   annotation_custom(grob = ggplotGrob(panel.density),
#'                     xmin = 1, xmax = 2300,
#'                     ymin = 8, ymax = 20)

