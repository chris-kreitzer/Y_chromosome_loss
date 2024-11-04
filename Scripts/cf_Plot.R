cf_plot = function(facets_data,
                   method = c('em', 'cncf'),
                   plotX = FALSE,
                   genome = c('hg19', 'hg18', 'hg38'),
                   return_object = FALSE) {
  
  genome = match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE)
  method = match.arg(method, c('em', 'cncf'), several.ok = FALSE)
  
  snps = as.data.frame(facets_data$snps)
  segs = as.data.frame(facets_data$segs)
  
  if (!plotX) {
    snps = subset(snps, chrom < 23)
    segs = subset(segs, chrom < 23)
  }
  
  snps = facetsSuite:::get_cum_chr_maploc(snps, genome)
  mid = snps$mid[names(snps$mid) %in% snps$snps$chrom]
  centromeres = snps$centromeres
  snps = snps$snps
  
  if (method == 'em') {
    cols = c(grDevices::colorRampPalette(c('white', 'steelblue'))(10), 'papayawhip')[round(10 * segs$cf.em + 0.501)]
    my_ylab = 'CF (EM)'
  } else if (method == 'cncf') {
    my_ylab = 'CF (CNCF)'
    cols = c(grDevices::colorRampPalette(c('white', 'steelblue'))(10), 'papayawhip')[round(10 * segs$cf + 0.501)]
  }
  
  starts = cumsum(c(1, segs$num.mark))[seq_along(segs$num.mark)]
  ends = cumsum(c(segs$num.mark))
  
  my_starts = snps[starts, 'chr_maploc']
  my_ends = snps[ends, 'chr_maploc']
  
  cf = ggplot(segs) +
    geom_rect(aes(xmin = my_starts, xmax = my_ends, ymax = 1, ymin = 0),
              fill = cols, col = 'white', size = 0) +
    scale_x_continuous(breaks = mid, labels = names(mid), expand = c(.01, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = my_ylab) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, size = 8, color = 'black'),
          axis.text.y = element_text(angle = 0, size = 8, color = 'white'),
          axis.ticks.y = element_line(color = 'white'),
          text = element_text(size = 10),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 1, .5, 0), 'lines'))
  
  if (return_object == TRUE) {
    cf 
  } else {
    suppressMessages(print(cf))
  }
}
