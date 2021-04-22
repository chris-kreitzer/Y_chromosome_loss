## PanCancer analysis of Y chromosome loss:
## Y_chromosome structure and coverage by IMPACT data
##
## Start (revision): 14/04/2021:

##-----------------------------------------------------------------------------
## To obtain the data for the Seq-Coverage across the Y-chromosome for 
## IMPACT Seq samples - see script: IMPACT-Y.coverage.R (run on Juno)
##
## IMPACT.features.covered comes from: IMPACT.Feature.Coverage.R (Juno)
##-----------------------------------------------------------------------------

rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')
source('Y_chromosome_loss/Plotting_theme.R')

## Libraries:
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggdist)

## Input:
IMPACT.coverage_df = read.csv('Data_out/IMPACT-Y.coverage.txt', sep = '\t') #' 10 random samples
IMPACT.features.covered_df = read.csv('Data_out/IMPACT.Features.covered.out.txt', sep = '\t')
IMPACT.breakpoints_df = read.csv('Data_out/IMPACT.breakpoints_out.txt', sep = '\t')
WES.breakpoints_df = read.csv('Data_out/WES.Breakpoints.Y.txt', sep = '\t')
Y.genes_df = read.csv('Data_out/Y_chromosome_genes.txt', sep = '\t')


## Processing:
#' summarise the IMPACT coverage data from 10 independent samples
IMPACT.coverage.summary = data.frame(tapply(IMPACT.coverage_df$SNPs.covered, 
                                            IMPACT.coverage_df$bin, mean))

IMPACT.coverage.summary$bin = row.names(IMPACT.coverage.summary)
colnames(IMPACT.coverage.summary) = c('Markers', 'bin')
row.names(IMPACT.coverage.summary) = NULL
IMPACT.coverage.summary$start = NA

for(i in 1:nrow(IMPACT.coverage.summary)){
  IMPACT.coverage.summary$start[i] = unique(IMPACT.coverage_df$start[which(IMPACT.coverage_df$bin == IMPACT.coverage.summary$bin[i])])
}

IMPACT.coverage.summary$sort = extract_numeric(IMPACT.coverage.summary$bin)
IMPACT.coverage.summary = IMPACT.coverage.summary[order(IMPACT.coverage.summary$sort), ]

#' how many protein_coding genes are directly covered by IMPACT
#' fetch the bin and make data frame for display underneath coverage plot
protein.coding.genes.covered = unique(IMPACT.features.covered_df$Gene)
Genes.selected.print = Y.genes_df[which(Y.genes_df$external_gene_name %in% protein.coding.genes.covered), ]
Genes.selected.print = Genes.selected.print[, c('hgnc_symbol', 'gene_biotype', 'start_position', 'end_position')]
colnames(Genes.selected.print) = c('Gene', 'Gene.type', 'Start', 'End')

#' assign bins to respective genes (mapping to chromosome)
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
Y.chromo.start = min(Y.chromo.MasterFile_df$loc)
Y.chromo.end = max(Y.chromo.MasterFile_df$loc)
Y.chromo.sequence = unlist(apply(data.frame(start = Y.chromo.start, end = Y.chromo.end), 1, Reduce, f = seq))

Y.gene.coverage_df = data.frame()
for(i in 1:nrow(Genes.selected.print)){
  print(i)
  if(apply(Genes.selected.print[i, c(3, 4)], 1, 
           function(x) any(seq(x[1], x[2]) %in% Y.chromo.sequence))){
    gene.out = data.frame(bin = Y.chromo.MasterFile_df$bin[which(Y.chromo.MasterFile_df$loc == Genes.selected.print[i, 'Start'])],
                          gene = Genes.selected.print[i, 'Gene'])
  } else next
  Y.gene.coverage_df = rbind(Y.gene.coverage_df, gene.out)
}

row.names(Y.gene.coverage_df) = NULL
colnames(Y.gene.coverage_df) = c('bin', 'Gene')

#' backfill every bin (1 through 5900): needed for plot
missing.bins = data.frame(Gene = NA,
                          bin = setdiff(unique(Y.chromo.MasterFile_df$bin), Y.gene.coverage_df$bin))
IMPACT.Feature.Y.plot_df = rbind(Y.gene.coverage_df[, c(1,2)], missing.bins)
IMPACT.Feature.Y.plot_df$sort = extract_numeric(IMPACT.Feature.Y.plot_df$bin)
IMPACT.Feature.Y.plot_df = IMPACT.Feature.Y.plot_df[order(IMPACT.Feature.Y.plot_df$sort), ]


#' How many protein coding genes are directly covered by IMPACT-Seq?
number.protein.coding.genesY = length(Y.genes_df$external_gene_name[which(Y.genes_df$gene_biotype == 'protein_coding')])

IMPACT.coverage.summary_df = data.frame()
for(i in unique(IMPACT.features.covered_df$ID)){
  rel.amount = data.frame(ID = i,
                          rel.covered = sum(IMPACT.features.covered_df$Match[which(IMPACT.features.covered_df$ID == i)]) / number.protein.coding.genesY)
  IMPACT.coverage.summary_df = rbind(IMPACT.coverage.summary_df, rel.amount)
}

##-----------------------------------------------------------------------------
## merge WES (different script) and IMPACT coverage for the plot
## ----------------------------------------------------------------------------
Feature.coverage.plot_df = data.frame(method = c('WES', 'IMPACT'),
                                      mean.coverage = c(mean(WES.coverage.summary_df$rel.covered), 
                                                        mean(IMPACT.coverage.summary_df$rel.covered)),
                                      sd.coverage = c(sd(WES.coverage.summary_df$rel.covered),
                                                      sd(IMPACT.coverage.summary_df$rel.covered)))

#' look into Marker-distribution (WES vs IMPACT) which retain for segmentation
#' create tribble
Markers.comparision = tribble(
  ~group, ~value,
  'WES', WES.breakpoints_df$breakpoints.used_facets,
  'IMPACT', IMPACT.breakpoints_df$breakpoints.used_facets
) %>%
  unnest(value)


## Visualization:
#' IMPACT coverage (Markers along the Y-chromosome):

#' Y chromosome coverage of 10 independent IMPACT-Seq samples:
breaks.max = max(IMPACT.coverage.summary$sort)
breaks.seq = seq(0, breaks.max, length.out = 6)
breaks.seq[1] = 1
y.max = max(as.integer(IMPACT.coverage.summary$Markers)) + 1

Impact.Y.coverage_plot = ggplot(IMPACT.coverage.summary, 
                    aes(x = sort, 
                        y = Markers, 
                        group = 1)) + 
  
  geom_line(size = 0.7) +
  
  annotate('segment',
           x = 50, xend = 5900 / 10,
           y = 38, yend = 38,
           size = 1.7) +
  
  annotate('text',
           x = 50, 
           y = 38,
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
  
  theme_Y(base_size = 14, base_rect_size = 0.1, base_line_size = 0.1) +
  
  theme(axis.text.y = element_blank(),
        panel.border = element_rect(size = 1.0, fill = NA),
        axis.ticks = element_blank()) +
  
  labs(title = 'IMPACT-Panel Seq. resolution across Y-chromosome', 
       subtitle = 'Peaks depict the number of markers in given bins (1bin = 1Kb) which retain for segmentation',
       x = '', y = 'Average number of markers')


#' Tile plot for protein coding genes which are directly covered by IMPACT
breaks.x = IMPACT.Feature.Y.plot_df$sort[!is.na(IMPACT.Feature.Y.plot_df$Gene)]
labels.x = c()
for(i in 1:length(breaks.x)){
  label = IMPACT.Feature.Y.plot_df$Gene[which(IMPACT.Feature.Y.plot_df$sort == breaks.x[i])]
  labels.x = c(labels.x, label)
}

#' make the geom_tile plot
IMPACT.feature.coverage.plot = ggplot(IMPACT.Feature.Y.plot_df,
                                      aes(x = sort,
                                          y = NA,
                                          fill = Gene)) +
  geom_tile(na.rm = T) +
  scale_fill_manual(values = c("USP9Y" = 'red',
                               "UTY" = 'red',
                               "KDM5D" = 'red',
                               "EIF1AY" = 'red',
                               "RPS4Y2" = 'red',
                               "SRY" = 'red',
                               "RPS4Y1" = 'red',
                               "ZFY" = 'red',
                               "PCDH11Y" = 'red',
                               "AMELY" = 'red',
                               "FAM197Y1" = 'red',
                               'NA' = 'white')) +
  
  scale_x_continuous(expand = c(0.01,0),
                     labels = labels.x, 
                     breaks = breaks.x, 
                     guide = guide_axis(n.dodge = 2)) +
  scale_y_discrete(expand = c(0, 0)) +
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
        legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = 'Genomic position', y = '')
  
Fig2A = Impact.Y.coverage_plot / IMPACT.feature.coverage.plot + 
  plot_layout(widths = unit(c(25, 25), c('cm', 'cm')),
              heights = unit(c(10, 1), c('cm', 'cm'))) 


ggsave_golden(filename = 'Figures/IMPACT-Y.coverage.pdf', plot = Fig2A, width = 16)


#' Feature coverage plot: comparing WES and IMPACT
#' Feature.coverage.plot_df
Feature.coverage.plot_df$method = factor(Feature.coverage.plot_df$method, levels = c('WES', 'IMPACT'))
Feature.coverage.comparison.plot = ggplot(Feature.coverage.plot_df, 
                                          aes(x = method, y = mean.coverage, fill = method)) +
  geom_bar(stat = 'identity', width = 0.9) +
  geom_linerange(ymin = Feature.coverage.plot_df$mean.coverage - Feature.coverage.plot_df$sd.coverage,
                ymax = Feature.coverage.plot_df$mean.coverage + Feature.coverage.plot_df$sd.coverage,
                size = 0.8) +
  scale_fill_manual(values = c('IMPACT' = '#A2DDE6',
                               'WES' = '#FBAE94'),
                    name = 'Seq-strategy') +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(0, max(Feature.coverage.plot_df$mean.coverage) + 
                                  max(Feature.coverage.plot_df$sd.coverage)),
                     breaks = c(0, 0.2, 0.4, 0.6),
                     labels = c(0, '20%', '40%', '60%')) +
  scale_x_discrete(expand = c(0.5, 0)) +
  theme_Y(base_size = 16) +
  theme(panel.border = element_blank(),
        legend.position = 'none',
        aspect.ratio = 2,
        plot.subtitle = element_text(size = 10)) +
  labs(x = 'Seq-strategy', y = 'Frequency',
       title = 'Protein coding genes covered\nacross the Y chromosome',
       subtitle = 'Data from 238 prostate samples')

ggsave_golden(filename = 'Figures/Protein.coding.genes.covered.pdf', 
              plot = Feature.coverage.comparison.plot,
              width = 12)


#' Marker-density comparison for WES and IMPACT
Marker.comparision.plot = ggplot(Markers.comparision, 
                                 aes(x = value, y = group, fill = group)) + 
  stat_halfeye(aes(thickness = stat(f * n)),
               adjust = 2,
               alpha = 0.6,
               height = 1,
               justification = -0.05) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0.1, 0)) +
  scale_fill_manual(values = c('WES' = '#A2DDE6',
                               'IMPACT' = '#FBAE94')) +
  geom_vline(xintercept = seq(0, 5000, 1000),
             linetype = 'dashed',
             color = 'grey55',
             size = 0.2) +
  
  annotate('curve',
           xend = median(WES.breakpoints_df$breakpoints.used_facets),
           x = 1450,
           yend = 2.05,
           y = 2.2,
           curvature = 0.25,
           size = 0.6, 
           arrow = arrow(length = unit(0.1, 'inches'), type = 'closed')) +
  
  annotate('text',
           x = 1460,
           y = 2.2,
           label = 'median=1,298',
           vjust = 0.5,
           hjust = 0,
           lineheight = 0.85,
           color = 'black',
           size = 3.5,
           family = "RobotoCondensed-Regular",
           fontface = 'plain') +
  annotate('curve',
           xend = median(IMPACT.breakpoints_df$breakpoints.used_facets),
           x = 500,
           yend = 1.05,
           y = 1.20,
           curvature = 0.25,
           size = 0.6, 
           arrow = arrow(length = unit(0.1, 'inches'), type = 'closed')) +
  
  annotate('text',
           x = 505,
           y = 1.20,
           label = 'median=380',
           vjust = 0.5,
           hjust = 0,
           lineheight = 0.85,
           color = 'black',
           size = 3.5,
           family = "RobotoCondensed-Regular",
           fontface = 'plain') +
  theme_Y(base_size = 18) +
  theme(legend.position = 'none',
        plot.subtitle = element_text(size = 12)) +
  labs(x = 'Number of markers', y = '', 
       title = 'WES provide significant more informative SNPs (markers) than IMPACT',
       subtitle = 'Data from 238 prostate samples')
  
  
ggsave_golden(filename = 'Figures/Breakpoint.comparision.pdf', 
              plot = Marker.comparision.plot, 
              width = 12)
  