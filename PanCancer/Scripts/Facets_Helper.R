##----------------+
## FacetsY helper functions
##----------------+
##
## start: 12/15/2022
## chris-kreitzer

library(data.table)
library(dplyr)


##----------------+
## Read the count-matrices
##----------------+
readsnpmatrix = function(input_file,
                           err.thresh = 10,
                           del.thresh = 10) {
  
  read_counts = data.table::fread(cmd = paste('gunzip -c', input_file), key = c('Chromosome', 'Position'))
  
  if (nrow(read_counts) == 0) { # necessary since fread command doesn't throw errors for certain cases
    stop(paste(input_file, 'does not exist or cannot be read properly.'), call. = F)
  }
  
  read_counts = read_counts[File1E <= err.thresh & File2E <= err.thresh &
                              File1D <= del.thresh & File2D <= del.thresh &
                              !Chromosome %in% c('MT', 'chrM')]
  
  read_counts[, `:=`(
    NOR.DP = File1R + File1A,
    TUM.DP = File2R + File2A,
    NOR.RD = File1R,
    TUM.RD = File2R,
    Chromosome = gsub('chr', '', Chromosome)
  )][, ('Chromosome') := factor(get('Chromosome'), levels = c(1:22, 'X', 'Y'))]
  
  read_counts[order(Chromosome, Position)][, list(Chromosome, Position, NOR.DP, NOR.RD, TUM.DP, TUM.RD)]
}


##----------------+
## Format output to IGV
##----------------+
#' function for data conversion to IGV 
format_igv_seg = function(facets_output, sample_id, normalize = TRUE) {
  
  if (!all(c('snps', 'segs', 'dipLogR') %in% names(facets_output))) {
    stop(paste0('Input is missing segs, snps or dipLogR ojbect.'), call. = FALSE)
  }
  
  seg = group_by(facets_output$snps, chrom, seg) %>% 
    summarize(loc.start = min(maploc),
              loc.end = max(maploc)) %>% 
    ungroup() %>% 
    left_join(., select(facets_output$segs, chrom, seg, num.mark, seg.mean = cnlr.median),
              by = c('chrom', 'seg')) %>% 
    mutate(ID = sample_id) %>% 
    select(ID, chrom, loc.start, loc.end, num.mark, seg.mean)
  
  if (normalize) { seg = mutate(seg, seg.mean = seg.mean - facets_output$dipLogR) }
  data.frame(seg)
}

#' out