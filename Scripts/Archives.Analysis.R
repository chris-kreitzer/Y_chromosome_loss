## PanCancer analysis of Y chromosome loss:
## Archives, Utilities and deprecated scripts
##


# Pileup = vroom::vroom('~/Documents/MSKCC/04_Y_chromo_loss/Data/FacetsPileupFile', delim = '\t', skip = 153)
# colnames(Pileup)[1] = 'Chrom'
# Pileup.Y_df = Pileup[which(Pileup$Chrom == 'chrY'), ]
# write.table(Pileup.Y_df, file = '~/Documents/MSKCC/04_Y_chromo_loss/Data/Y_chromosome_Facets_Pileup.txt', sep = '\t', row.names = F)

# MSK.WES_df = read.csv('~/Documents/MSKCC/04_Y_chromo_loss/s_C_740W28_M001_d__s_C_740W28_N001_d.snp_pileup.gz', sep = ',')
# MSK.WES_Y_df = MSK.WES_df[which(MSK.WES_df$Chromosome == 'Y'), ]
# write.table(MSK.WES_Y_df, file = '~/Documents/MSKCC/04_Y_chromo_loss/Data/Y_chromosome_IMPACT_WES.example.txt', sep = '\t', row.names = F)

# MSK.Panel_df = read.csv('~/Desktop/mnt/ATMcountdata/countsMerged____P-0013697-T01-IM5_P-0013697-N01-IM5.dat.gz', sep = ',')
# MSK.Panel_Y_df = MSK.Panel_df[which(MSK.Panel_df$Chromosome == 'Y'), ]
# write.table(MSK.Panel_Y_df, file = 'Data/Y_chromosome_IMPACT_Panel.example.txt', sep = '\t', row.names = F)


#' Seq. coverage across the Y-chromosome: Run that on Juno
#' run the same analysis on n = 15 randomly selected IMPACT Panel sequenced samples
Panel.sequencing.examples # select random samples

Y_coverage = data.frame()
for(i in sample(Panel.sequencing.examples, 2, replace = T)){
  print(i)
  input = vroom::vroom(i)
  input = input[which(input$Chromosome == 'Y'), ]
  input.modi = Facets.preprocess(data = input)
  data.merge = merge(Y.chromo.MasterFile_df, 
                            input.modi[, c('Position', 'Normal_Coverage', 'Tumor_Coverage')], 
                            by.x = 'loc', 
                            by.y = 'Position',
                            all.x = T)
  
  data.merge[is.na(data.merge)] = 0
  message('Analyse bins')
  data.merge.out = data.frame(bin = unique(data.merge$bin),
                              start = unlist(lapply(unique(data.merge$bin), function(x) min.function(data = data.merge, bin = x))),
                              end = unlist(lapply(unique(data.merge$bin), function(x) max.function(data = data.merge, bin = x))),
                              SNPs.covered = unlist(lapply(unique(data.merge$bin), function(x) coverage.count.function(data = data.merge, bin = x))))
  
  Y_coverage = rbind(Y_coverage, data.merge.out)
}
  
write.table(Y_coverage, file = 'Data_out/Ave.coverageY_PanelSeq.txt', sep = '\t', row.names = F)

#' select only one feature per bin
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

## Specific annotations for Y chromosomes
## different tracks
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





## obsolete scripts:

#' Y chromosome structure:
#' Y chromosome, split into n = bin.numer (input above)
#' full.y.chromo = ggplot() + 
#'   scale_x_continuous(limits = c(0, bin.number), expand = c(0,0)) +
#'   scale_y_continuous(limits = c(-1, 1.5), expand = c(0, 0)) +
#'   theme(panel.border = element_blank(),
#'         panel.background = element_blank(),
#'         panel.grid = element_blank(),
#'         axis.text = element_blank(),
#'         axis.ticks = element_blank(),
#'         axis.title = element_blank(),
#'         plot.margin = unit(c(0,0,0,0), 'cm'))
#' 
#' ## add PAR1 region: 
#' #' PAR1 start = 10,001 == bin:1
#' #' PAR1 end = 2649520 == Y.chromo.Masterfile_df$bin[which(Y.chromo.MasterFile_df$loc == 2649520)]
#' PAR1.start = 1
#' PAR1.end = unique(Y.chromo.MasterFile_df$bin[which(Y.chromo.MasterFile_df$loc == 2649520)])
#' PAR1.end = as.numeric(extract_numeric(PAR1.end))
#' 
#' full.y.chromo = full.y.chromo + statebins:::geom_rrect(mapping = aes(xmin = PAR1.start,
#'                                                                      xmax = PAR1.end,
#'                                                                      ymin = -1,
#'                                                                      ymax = 1),
#'                                                        color = NA,
#'                                                        fill = '#527ba4',
#'                                                        radius = grid:::unit(10, 'pt')) +
#'   annotate('text', x = (PAR1.start + PAR1.end) / 2, 
#'            y = 1.0, 
#'            label = 'PAR1',
#'            hjust = -0.5,
#'            vjust = 0.5,
#'            size = 5,
#'            angle = 90)
#' 
#' ## add MSY region:
#' MSY.start = PAR1.end
#' MSY.end = unique(Y.chromo.MasterFile_df$bin[which(Y.chromo.MasterFile_df$loc == 59034049)])
#' MSY.end = bin.number - 1
#' 
#' full.y.chromo = full.y.chromo + 
#'   statebins:::geom_rrect(mapping = aes(xmin = MSY.start,
#'                                        xmax = MSY.end,
#'                                        ymin = -1,
#'                                        ymax = 1),
#'                          color = NA,
#'                          fill = '#e0a347',
#'                          radius = grid:::unit(10, 'pt')) +
#'   
#'   annotate('text', x = (MSY.start + MSY.end) / 2, 
#'            y = 1.0, 
#'            label = 'MSY',
#'            hjust = 0.5,
#'            vjust = -1,
#'            size = 5)
#' 
#' ## add PAR2 region:
#' PAR2.start = MSY.end
#' PAR2.end = bin.number
#' 
#' full.y.chromo = full.y.chromo + statebins:::geom_rrect(mapping = aes(xmin = PAR2.start,
#'                                                                      xmax = PAR2.end,
#'                                                                      ymin = -1,
#'                                                                      ymax = 1),
#'                                                        color = NA,
#'                                                        fill = '#527ba4',
#'                                                        radius = grid:::unit(2, 'pt')) +
#'   annotate('text', x = (PAR2.start + PAR2.end) / 2, 
#'            y = 1.0, 
#'            label = 'PAR2',
#'            hjust = -0.5,
#'            vjust = 0,
#'            size = 5,
#'            angle = 90)

