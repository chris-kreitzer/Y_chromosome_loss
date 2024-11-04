##----------------+
## CDKN2A exonic structure
## mapping maploc and cnlr
##----------------+


## Exonic structure of CDKN2A
# https://grch37.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000147889;mr=9:21969729-21969792;r=9:21967751-21995300;t=ENST00000304494
# Transcript ID: ENST00000304494.5

# Exon1: 21,975,097	21,974,677
# Exon2: 21,971,207	21,970,901
# Exon3: 21,968,241	21,967,752


exon_mapping = function(data){
  cdkn2a = data[which(data$maploc >= 21967752 & data$maploc <= 21995300), ]
  cdkn2a$seq = NA
  k = 0
  for(i in 1:nrow(cdkn2a)){
    cdkn2a$seq[i] = k
    k = k + 0.1
  }
  cdkn2a$seq = cdkn2a$seq * -1
  midpoint = median(cdkn2a$cnlr)
  
  
  ##-- plot
  title_size = 14
  axis_label_size = 12
  text_size = 12
  
  exon = ggplot(cdkn2a) +
    geom_rect(aes(xmin = 21975097, xmax = 21974677, ymin = 0.5, ymax = 1.5)) +
    geom_rect(aes(xmin = 21971207, xmax = 21970901, ymin = 0.5, ymax = 1.5)) +
    geom_rect(aes(xmin = 21968241, xmax = 21967752, ymin = 0.5, ymax = 1.5)) +
    geom_rect(aes(xmin = 21994138, xmax = 21994489, ymin = 0.5, ymax = 1.5)) +
    geom_point(aes(x = maploc, y = seq, color = cnlr), size = 3.5) +
    scale_colour_gradient2(
      low = "darkblue",
      mid = "white",
      high = "red",
      midpoint = midpoint,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour", name = 'CnLR') +
    
    theme_bw() +
    geom_vline(xintercept = 21975097, color = 'grey85', linetype = 'dashed') +
    geom_vline(xintercept = 21974677, color = 'grey85', linetype = 'dashed') +
    geom_vline(xintercept = 21971207, color = 'grey85', linetype = 'dashed') +
    geom_vline(xintercept = 21970901, color = 'grey85', linetype = 'dashed') +
    geom_vline(xintercept = 21968241, color = 'grey85', linetype = 'dashed') +
    geom_vline(xintercept = 21967752, color = 'grey85', linetype = 'dashed') +
    geom_vline(xintercept = 21994138, color = 'grey85', linetype = 'dashed') +
    geom_vline(xintercept = 21994489, color = 'grey85', linetype = 'dashed') +
    
    annotate('text',
      x = (xmin = 21975097 + 21974677) / 2,
      y = 1.8,
      size = 4.5,
      label = 'Exon1') +
    
    annotate('text',
             x = (xmin = 21994138 + 21994489) / 2,
             y = 1.8,
             size = 4.5,
             label = 'Exon1 alpha\nCDKN2B') +
    
  
  annotate('text',
             x = (21971207 + 21970901) / 2,
             y = 1.8,
             size = 4.5,
             label = 'Exon2') +
    
    annotate('text',
             x = (21968241 + 21967752) / 2,
             y = 1.8,
             size = 4.5,
             label = 'Exon3') +
    
    labs(x = 'genomic coordinates',
         y = element_blank(),
         title = 'CDKN2A exonic structure\nENST00000304494.5') +
    
    theme(text = element_text(size = text_size),
          plot.title.position = 'plot',
          panel.grid = element_blank(),
          axis.title.x = element_text(size = axis_label_size),
          axis.text.x = element_text(size = axis_label_size),
          axis.ticks.x = element_line(size = 0.25),
          axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.x = element_blank()
    )
  
  return(exon)
  
}


#' out

