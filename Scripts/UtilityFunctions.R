## UTILITY functions for CSF analysis





## multi_snp-pileup read for FACETS analysis
## 
## start: 08/31/2022
## chris-kreitzer


multi_readSnpMatrix = function(filename, 
                               tumor_sample = 2,
                               skip = 0L, 
                               err.thresh = 10, 
                               del.thresh = 10){
  
  pileup = read.csv(filename,
                    stringsAsFactors = FALSE)
  pileup[, 1] = as.character(pileup[,1])
  pileup[, 2] = as.numeric(as.integer(pileup[, 2]))
  pileup[, 3] = as.character(pileup[, 3])
  pileup[, 4] = as.character(pileup[, 4])
  cols = seq(5, length(pileup), 1)
  pileup[, cols] = apply(pileup[, cols], 2, function(x) as.numeric(as.integer(x)))
  
  if(tumor_sample - 1 > length(cols) - 4 / 4){
    stop(paste('wrong tumor number (column) selected'), call. = F)
  }
  
  #' select the right tumor sample:
  pileup_tumor = cbind(pileup[,c(1:8)], pileup[, grep(pattern = tumor_sample, x = colnames(pileup))])
  
  if(tumor_sample != 2){
    colnames(pileup_tumor)[9:12] = gsub(pattern = tumor_sample, replacement = '2', x = colnames(pileup_tumor)[9:12])
  }
  
  # remove chr if present in Chrom
  if(grepl("chr", pileup_tumor$Chromosome[1])){
    pileup_tumor$Chromosome = gsub("chr", "", pileup_tumor$Chromosome)
  }
  
  message(paste0((length(cols) - 4) / 4, ' tumor samples in pileup file; ready to analyze'))
  
  if(nrow(pileup_tumor) == 0){
    stop(paste(filename, 'does not exist or cannot be read properly.'), call. = F)
  }

  # remove loci where errors and deletions exceeded thresholds
  ii = which(pileup_tumor$File1E <= err.thresh & 
               pileup_tumor$File1D <= del.thresh & 
               pileup_tumor$File2E <= err.thresh & 
               pileup_tumor$File2D <= del.thresh)
  rcmat = pileup_tumor[ii, 1:2]
  rcmat$NOR.DP <- pileup_tumor$File1R[ii] + pileup_tumor$File1A[ii]
  rcmat$NOR.RD <- pileup_tumor$File1R[ii]
  rcmat$TUM.DP <- pileup_tumor$File2R[ii] + pileup_tumor$File2A[ii]
  rcmat$TUM.RD <- pileup_tumor$File2R[ii]
  
  return(rcmat)
}


##-----------------------------------------------------------------------------
## Customized CnLR plot: based on FacetsSuite
##
## start: 04/22/2022
## 
## chris-kreitzer


source('~/Documents/GitHub/MSKCC/Scripts/hg19.R')


# cnlr_plot(facets_data = facets_data, genome = hg19)

#' function
cnlr_plot = function(facets_data,
                     colors = c('#5c9e75', '#272425'),
                     highlight_gene = NULL) {
  
  if('hg19' %in% ls(name = globalenv())){
    print('hg19 gene model is used')
  } else {
    stop('hg19 not loaded in current environment')
  }
  
  genome = hg19
  
  snps = facets_data$snps
  segs = facets_data$segs
  dipLogR = facets_data$dipLogR
  
  snps = get_cum_chr_maploc(snps, genome)
  mid = snps$mid[names(snps$mid) %in% snps$snps$chrom]
  centromeres = snps$centromeres
  snps = snps$snps
  
  snps$cnlr_median = rep(segs$cnlr.median, segs$num.mark)
  
  starts = cumsum(c(1, segs$num.mark))[seq_along(segs$num.mark)]
  ends = cumsum(c(segs$num.mark))
  my_starts = snps[starts, c('chr_maploc', 'cnlr_median')]
  my_ends = snps[ends, c('chr_maploc', 'cnlr_median')]
  
  ymin = floor(min(segs$cnlr.median, na.rm = T))
  ymax = ceiling(max(segs$cnlr.median, na.rm = T))
  if (ymin > -5) ymin = -5
  if (ymax < 5) ymax = 5
  
  
  pt_cols = colors[c(snps$chrom %% 2) + 1]
  
  # plot
  cnlr = ggplot(snps) +
    geom_point(aes(y = cnlr, x = chr_maploc), pch = 19, col = pt_cols, size = .4) +
    scale_x_continuous(breaks = mid, labels = names(mid), expand = c(0, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(ymin, ymax), expand = c(0,0)) +
    geom_hline(yintercept = dipLogR, color = 'sandybrown', size = .8) +
    geom_segment(data = segs, aes(x = my_starts$chr_maploc, xend = my_ends$chr_maploc,
                                  y = my_starts$cnlr_median, yend = my_ends$cnlr_median),
                 col = 'red3', size = 1, lineend = 'butt') +
    labs(x = NULL, y = 'Copy number log ratio') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, size = 10, color = 'black'),
          axis.text.y = element_blank(),
          text = element_text(size = 12),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_line(colour = 'grey', size = 0.2),
          panel.grid.major = element_blank(),
          panel.border = element_rect(fill = NA, size = 1.5))
  
  
  # if highligthing gene
  if (!is.null(highlight_gene)) {
    if (is.character(highlight_gene)) {
      highlight_gene = get_gene_position(highlight_gene)
    }
    snps$gene = FALSE
    snps$gene[which(snps$chrom %in% highlight_gene$chrom &
                      snps$chr_maploc >= highlight_gene$start &
                      snps$chr_maploc <= highlight_gene$end)] = TRUE
    cnlr = cnlr +
      geom_vline(xintercept = highlight_gene$mid, color = 'palevioletred1') +
      geom_point(data = filter(snps, gene == TRUE), aes(y = cnlr, x = chr_maploc), color = '#525252', size = .4) 
  }
  return(cnlr)
  
}



get_cum_chr_maploc = function(snps, genome = hg19) {
  
  genome = genome
  
  cum_chrom_lengths = cumsum(as.numeric(genome$size))
  mid = cum_chrom_lengths - (genome$size / 2)
  names(mid) = seq_len(nrow(genome))
  centromeres = genome$centromere + c(0, cum_chrom_lengths[-length(cum_chrom_lengths)])
  
  snps$chr_maploc = snps$maploc + c(0, cum_chrom_lengths)[snps$chrom]
  
  list(snps = snps, mid = mid, centromeres = centromeres)
}


#' extract gene position
get_gene_position = function(genes) {
  
  if('genes_hg19' %in% ls(name = globalenv())){
    message('genes_hg19 are used')
  } else {
    stop('No gene_models are loaded in global environment (eg. genes_hg19')
  }
  
  gene_loci = genes_hg19 
  genome = hg19
  cum_chrom_lengths = cumsum(as.numeric(genome$size))
  
  if(!all(genes %in% gene_loci$gene)) {
    missing_genes = setdiff(genes, gene_loci$gene)
    stop(paste('Gene(s)', paste(missing_genes, collapse = ', '), 'are not in gene annotation.'), call. = FALSE)
  } else {
    filter(gene_loci, gene %in% genes) %>%
      mutate(chrom = as.numeric(chrom),
             mid = start + (end-start) / 2,
             start = start + cum_chrom_lengths[chrom-1],
             end = end + cum_chrom_lengths[chrom-1],
             mid = mid + cum_chrom_lengths[chrom-1]) 
  }
}


##-----------------
## Save Facets-Output
##-----------------
create_facets_output = function(facets_output, 
                                directory, 
                                sample_id, 
                                counts_file, 
                                sel_run_type, 
                                run_details) {
  
  library(patchwork)
  library(dplyr)
  
  output_prefix = paste0(directory, '/', sample_id)
  dir.create(path = output_prefix)
  
  #' create cncf.txt
  facets_output$segs %>%
    mutate(ID = sample_id,
           loc.start = start,
           loc.end = end) %>%
    select(ID, chrom, loc.start, loc.end, seg, num.mark, nhet, cnlr.median, mafR, 
           segclust, cnlr.median.clust, mafR.clust, cf, tcn, lcn, 
           cf.em, tcn.em, lcn.em) %>%
    write.table(file = paste0(output_prefix, '/', basename(output_prefix), '.cncf.txt'), 
                quote = F, 
                row.names = F, 
                sep = '\t')
  
  #' create .Rdata
  out =
    list(
      jointseg = facets_output$snps,
      out = facets_output$segs %>% 
        select(chrom, seg, num.mark, nhet, cnlr.median, mafR, 
               segclust, cnlr.median.clust, mafR.clust, cf, tcn, lcn),
      nX = 23,
      chromlevels = c(1:22, "X"),
      dipLogR = facets_output$dipLogR,
      alBalLogR = facets_output$alBalLogR,
      IGV = NULL
    )
  fit = 
    list(
      loglik = facets_output$loglik,
      purity = facets_output$purity,
      ploidy = facets_output$ploidy,
      dipLogR = facets_output$dipLogR,
      seglen = -1,
      cncf = facets_output$segs,
      emflags = facets_output$em_flags
    )
  
  save(out, fit, file = paste0(output_prefix, '/', basename(output_prefix), ".Rdata"), compress = T)
  
  #' create .CNCF.png
  cnlr = facetsSuite::cnlr_plot(facets_output, return_object = T)
  valor = facetsSuite::valor_plot(facets_output, return_object = T)
  icn = facetsSuite::icn_plot(facets_output, return_object = T)
  cf = facetsSuite::cf_plot(facets_output, return_object = T)
  
  aggregate_fit = cnlr/valor/icn/cf
  ggsave(filename = paste0(output_prefix, '/', basename(output_prefix),  ".pdf"), plot = aggregate_fit)
}