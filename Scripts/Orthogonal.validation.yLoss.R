## PanCancer analysis of Y chromosome loss:
## Orthogonal validation of IMPACT WES and Panel sequencing results (FacetsY)
## 
## using the mLRR of Y-specific (protein coding) genes to confirm a proper loss:
## the rational: at positions where sequencing read depth in N and T are equal, we would roughly expect a ratio of 1
## say: Y-chromosome: Position 26,234,951 we have 10 reads in the Normal and 10 read in the Tumor sample --> log2(10/10) == 0
## gross deviations from 0 [- << ] indicate deep deletions
## 
## the mLRR is broadly used in literature to call Y chromosome losses
## 
## start: 18/03/2021
## revision: 29/03/2021


rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')


## Libraries:
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(cowplot))


## Input:
WES.prostate.cnlr_df = read.csv('../Data_out/Cnlr.WES.prostate.out.txt', sep = '\t')
WES.prostate.Y.calls_df = read.csv('../Data_out/WES.Prostate.Y_calls_txt', sep = '\t')
WES.Y.loss.calls_df = read.csv('../Data_out/Prostate_Y_loss.txt', sep = '\t')

#' ENSEMBL gene model: hg19
ensembl = useEnsembl(biomart = 'ensembl', 
                     dataset = 'hsapiens_gene_ensembl', 
                     GRCh = 37)
Y_genes = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 
                             'chromosome_name', 'start_position', 
                             'end_position'), 
                filters = 'chromosome_name', 
                values = "Y", 
                mart = ensembl)



## Functions:
#' fetch LRR of individual genes and calculate mean (different cohorts; Y_loss OR Y_intact)
LRR.among.genes = function(data, group = NULL){
  
  if(is.null(group)){
    stop('Error: Provide a grouping names')
  }
  if(!exists('GOI.mLRR.prostate')){
    stop('Error: GOI.mLRR.prostate must be defined in environment')
  }
  
  Y.intact.mLRR = data.frame()
  id = unique(data$ID)
  
  for(i in unique(id)){
    data.selected = data[which(data$ID == i),, drop = F]
    
    for(j in 1:nrow(GOI.mLRR.prostate)){
      cnlr.gene = data.selected[which(data.selected$maploc >= GOI.mLRR.prostate[j, 'Start'] & 
                                        data.selected$maploc <= GOI.mLRR.prostate[j, 'End']), ]
      
      if(nrow(cnlr.gene) <= 1){
        m.cnlr = 0
        SNPs = 0
      } 
      else if(nrow(cnlr.gene) > 1 & nrow(cnlr.gene) < 10){
        m.cnlr = mean(cnlr.gene$cnlr, na.rm = T)
        SNPs = nrow(cnlr.gene)
      }
      else if(nrow(cnlr.gene) >= 10){
        m.cnlr = median(cnlr.gene$cnlr, na.rm = T)
        SNPs = nrow(cnlr.gene)
      }
      
      gene = GOI.mLRR.prostate[j, 'Gene']
      
      data.merge = data.frame(ID = i,
                              gene = gene,
                              mLRR = m.cnlr,
                              group = group,
                              SNPs = SNPs)
      
      Y.intact.mLRR = rbind(Y.intact.mLRR, data.merge)
    }
  }
  return(Y.intact.mLRR)
}


#' orthogonal gene summary
gene.summary = function(data, genes, group){
  data.out_df = data.frame()
  for(i in unique(genes)){
    gene.mean = mean(data$mLRR[which(data$gene == i)], na.rm = T)
    gene.sd = sd(data$mLRR[which(data$gene == i)], na.rm = T)
    
    group = group
    gene = i
    
    out = data.frame(group = group,
                     gene = gene,
                     mean = gene.mean,
                     sd = gene.sd)
    
    data.out_df = rbind(data.out_df, out)
  }
  return(data.out_df)
}



## Processing:
Y.intact.samples = unique(WES.Y.loss.calls_df$WES.ID[which(WES.Y.loss.calls_df$WES.call == 'intact_Y_chrom')])
Y.loss.samples = unique(WES.Y.loss.calls_df$WES.ID[which(WES.Y.loss.calls_df$WES.call == 'Y_chrom_loss')])
colnames(Y_genes) = c('ID', 'Gene', 'Chromosome', 'Start', 'End')

#' Genes of interest; to look into mLRR (due to direct coverage of sequencing)
GOI.mLRR.prostate = c('KDM5D',
                      'RPS4Y1',
                      'NAP1L1P2',
                      'ZFY',
                      'ZFY-AS1',
                      'LINC00278',
                      'PCDH11Y',
                      'EIF4A1P2',
                      'RPL26P37',
                      'TBL1Y',
                      'USP9Y',
                      'UTY',
                      'PSMA6P1',
                      'EIF1AY',
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

GOI.mLRR.prostate = Y_genes[which(Y_genes$Gene %in% GOI.mLRR.prostate),, drop = F]


## Analysis:
#' Illustrative example: comparisons between Median cnlr values in INTACT and LOSS Y-chromosome group
#' here: loss
loss_extreme_cnlr = data.frame()
for(i in unique(Y.loss.samples)){
  median.sample = summary(WES.prostate.cnlr_df[which(WES.prostate.cnlr_df$ID == i), 'cnlr'])
  median.sample = median.sample[['Median']]
  out = cbind(median.sample, i)
  loss_extreme_cnlr = rbind(loss_extreme_cnlr, out)
}
loss_extreme_cnlr$median.sample = as.numeric(as.character(loss_extreme_cnlr$median.sample))
Y_loss.example.cnlr = WES.prostate.cnlr_df[which(WES.prostate.cnlr_df$ID == loss_extreme_cnlr$i[which.min(loss_extreme_cnlr$median.sample)]), ]

#' here: intact Y-chromsomes:
intact_extreme_cnlr = data.frame()
for(i in unique(Y.intact.samples)){
  median.sample = summary(WES.prostate.cnlr_df[which(WES.prostate.cnlr_df$ID == i), 'cnlr'])
  median.sample = median.sample[['Median']]
  out = cbind(median.sample, i)
  intact_extreme_cnlr = rbind(intact_extreme_cnlr, out)
}
intact_extreme_cnlr$median.sample = as.numeric(as.character(intact_extreme_cnlr$median.sample))
Y_intact.example.cnlr = WES.prostate.cnlr_df[which(WES.prostate.cnlr_df$ID == intact_extreme_cnlr$i[which.max(intact_extreme_cnlr$median.sample)]), ]

#' create tribble, and prepare for plotting:
Cnlr.comparision = tribble(
  ~group, ~value,
  'Y-chromosome loss', Y_loss.example.cnlr$cnlr,
  'Y-chromosome intatct', Y_intact.example.cnlr$cnlr
) %>%
  unnest(value)


#' look into LRR of individual genes (orthogonal validation)
#' genes from table above
#' Y intact samples
WES.cnlr.intact_df = WES.prostate.cnlr_df[which(WES.prostate.cnlr_df$ID %in% Y.intact.samples), ]
Intact.processed.gene.out = LRR.among.genes(data = WES.cnlr.intact_df, group = 'intact')

#' Y loss samples
WES.cnlr.loss_df = WES.prostate.cnlr_df[which(WES.prostate.cnlr_df$ID %in% Y.loss.samples), ]
Loss.processed.gene.out = LRR.among.genes(data = WES.cnlr.loss_df, group = 'loss')


#' make summary for gene and prepare for plotting
gene.loss = Loss.processed.gene.out$gene[which(Loss.processed.gene.out$mLRR != 0)]
gene.loss = unique(gene.loss)

gene.summary.loss = gene.summary(data = Loss.processed.gene.out, genes = gene.loss, group = 'loss')
gene.summary.intact = gene.summary(data = Intact.processed.gene.out, genes = gene.loss, group = 'intact')
Combined.gene.cnlr_df = rbind(gene.summary.intact, gene.summary.loss)

#' make a t.test to directly compare the individual genes within the two groups:
statistics.gene = data.frame()
for(i in unique(gene.loss)){
  gene.pvalue = t.test(Intact.processed.gene.out$mLRR[which(Intact.processed.gene.out$gene == i)], 
                       Loss.processed.gene.out$mLRR[which(Loss.processed.gene.out$gene == i)])$p.value
  out = cbind(gene.pvalue,
              i)
  statistics.gene = rbind(statistics.gene, out)
}

colnames(statistics.gene) = c('pvalue', 'gene')
statistics.gene$group = NA


## Visualization:
#' orthogonal validation of Y chromsome loss
#' select four illustrative genes for plotting:
Combined.gene.cnlr_df = Combined.gene.cnlr_df[Combined.gene.cnlr_df$gene %in% c('RPS4Y1', 'KDM5D', 'DDX3Y', 'EIF1AY'), ]
Combined.gene.cnlr_df$gene = factor(Combined.gene.cnlr_df$gene, 
                                    levels = c('KDM5D', 'DDX3Y', 'EIF1AY', 'RPS4Y1'))
statistics.gene = statistics.gene[statistics.gene$gene %in% c('KDM5D', 'DDX3Y', 'EIF1AY', 'RPS4Y1'), ]
statistics.gene$pvalue = as.numeric(as.character(statistics.gene$pvalue))

#' make the plot:
gene.cnlr.comparision = ggplot(Combined.gene.cnlr_df, 
                               aes(x = gene, y = mean, 
                                   color = group)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd), 
                  position = position_dodge(0.1)) +
  geom_hline(yintercept = 0, 
             linetype = 'dashed',
             color = 'grey55',
             size = 0.2) +
  geom_text(data = statistics.gene, aes(x = gene, 
                                        y = 0.45, 
                                        label = paste0('p = ', round(pvalue, 6))),
            color = 'black', size = 4.5) +
  
  scale_x_discrete(expand = c(0.1, 0)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(-3, 1),
                     breaks = c(-2, -1, 0)) +
  scale_color_manual(values = c('intact' = '#F04E33',
                                'loss' = '#006CA2'),
                     name = 'Group',
                     guide = guide_legend(direction = 'horizontal',
                                          title.theme = element_text(size = 16, face = 'bold', family = 'RobotoCondensed-Regular'),
                                          label.theme = element_text(size = 14, family = 'RobotoCondensed-Regular'))) +
  theme(text = element_text(family = 'RobotoCondensed-Regular', size = 13),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.15, 0.96)) +
  panel_border(color = 'black') +
  labs(x = '', y = 'CnLR [Normal/Tumor]')

ggsave_golden(filename = 'Figures/CnLR.among.genes.pdf', plot = gene.cnlr.comparision, width = 10)


#' Illustrative example: compare groups:
CnLR.toy.plot = ggplot(Cnlr.comparision, aes(x = value, y = group, fill = group)) + 
  stat_halfeye(aes(thickness = stat(f * n)),
               adjust = 2,
               alpha = 0.6,
               height = 1,
               justification = -0.05) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(-6, -4, -2, 0, 2, 4)) +
  scale_y_discrete(expand = c(0.1, 0)) +
  scale_fill_manual(values = c('Y-chromosome loss' = '#A2DDE6',
                               'Y-chromosome intatct' = '#FBAE94')) +
  geom_vline(xintercept = seq(-6, 2, 2),
             linetype = 'dashed',
             color = 'grey55',
             size = 0.2) +
  annotate('curve',
           xend = median(Y_loss.example.cnlr$cnlr),
           x = -4,
           yend = 2.05,
           y = 2.2,
           curvature = 0.25,
           size = 0.6, 
           arrow = arrow(length = unit(0.1, 'inches'), type = 'closed')) +
  
  annotate('text',
           x = -3.95,
           y = 2.2,
           label = 'Median\n(n=488)',
           vjust = 0.5,
           hjust = 0,
           lineheight = 0.85,
           color = 'black',
           size = 3.5,
           family = "RobotoCondensed-Regular",
           fontface = 'bold.italic') +
  annotate('curve',
           xend = median(Y_intact.example.cnlr$cnlr),
           x = 1.2,
           yend = 1.05,
           y = 1.20,
           curvature = 0.25,
           size = 0.6, 
           arrow = arrow(length = unit(0.1, 'inches'), type = 'closed')) +
  annotate('text',
           x = 1.25,
           y = 1.20,
           label = 'Median\n(n=55)',
           vjust = 0.5,
           hjust = 0,
           lineheight = 0.85,
           color = 'black',
           size = 3.5,
           family = "RobotoCondensed-Regular",
           fontface = 'bold.italic') +
  annotate('text',
           x = -4,
           y = 2.9,
           color = '#006CA2',
           label = 'Y-chromosome loss',
           vjust = 0.5,
           hjust = 0,
           lineheight = 0.85,
           color = 'black',
           size = 6.5,
           family = "RobotoCondensed-Regular",
           fontface = 'bold.italic') +
  annotate('text',
           x = 1.0,
           y = 1.9,
           color = '#F04E33',
           label = 'Y-chromosome intact',
           vjust = 0,
           hjust = 0.0,
           lineheight = 0.85,
           color = 'black',
           size = 6.5,
           family = "RobotoCondensed-Regular",
           fontface = 'bold.italic') +
  theme(text = element_text(family = 'RobotoCondensed-Regular'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none',
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  panel_border(color = 'black') +
  labs(x = 'CnLR [Normal/Tumor]',
       title = 'Copy Number Log Ratio at individual bases across the Y-chromosome',
       subtitle = 'The median CnLR between the two groups can be used as surrogate for FacetsY calls')

ggsave_golden(filename = 'Figures/ToyExample_CnLR.pdf', plot = CnLR.toy.plot, width = 10)
           
