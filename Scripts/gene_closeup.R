##----------------+
## gene closeup plot
## for CSF
##----------------+
library(cowplot)
gene_closeup = function(data, gene, diplogr = NULL){
  snps = data$snps
  chrom = genes_hg19$chrom[which(genes_hg19$gene == gene)]
  start = genes_hg19$start[which(genes_hg19$gene == gene)] 
  end = genes_hg19$end[which(genes_hg19$gene == gene)]
  
  arm = ifelse(start <= hg19$centromere[which(hg19$chrom == chrom)], 'p', 'q')
  
  snps = snps[which(snps$chrom == chrom), ]
  
  if(arm == 'p'){
    snps = snps[which(snps$maploc <= hg19$centromere[which(hg19$chrom == chrom)]), ]
    snps$order = seq(1, nrow(snps), 1)
    snps$gene = ifelse(snps$maploc >= start & snps$maploc <= end, 'color', 'no_color')
  } else {
    snps = snps[which(snps$maploc >= hg19$centromere[which(hg19$chrom == chrom)]), ]
    snps$order = seq(1, nrow(snps), 1)
    snps$gene = ifelse(snps$maploc >= start & snps$maploc <= end, 'color', 'no_color')
  }
  
  plot = ggplot(snps, aes(x = order, y = cnlr, color = gene)) +
    geom_point(position = position_dodge(width = 0.1), size = 1.5) +
    geom_segment(aes(x = min(snps$order[which(snps$gene == 'color')]),
                     xend = max(snps$order[which(snps$gene == 'color')]),
                     y = mean(snps$cnlr[which(snps$gene == 'color')]),
                     yend = mean(snps$cnlr[which(snps$gene == 'color')])),
                 color = 'red', size = 1.5) +
    scale_color_manual(values = c('color' = '#487EFB', 
                                  'no_color' = 'grey75'),
                       name = 'Probes', labels = c(gene, paste0('chrom ', chrom, arm))) +
    scale_y_continuous(expand = c(0.01, 1), limits = c(-4.5, 4.5)) +
    geom_hline(yintercept = ifelse(is.null(diplogr), 0, diplogr), linetype = 'dashed') +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text = element_text(size = 12, colour = 'black')) +
    panel_border(size = 2, color = 'black') +
    labs(x = paste0('Chromosome ', chrom, arm), y = 'CnLR', title = gene)
  
  return(list(snps = snps, 
              plot = plot))
  
}

