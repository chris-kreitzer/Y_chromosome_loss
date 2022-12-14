##----------------+
## Assess fraction of seq. 
## reads uniquely mapping to 
## chromosome Y features
##----------------+
## Furthermore check the 
## sequence quality at
## selected features (n=429)
## Inspired by Ed Reznik 
## (mutation-phasing) paper
##----------------+
## start: 06/21/2022
## revision: 06/28/2022
## revision: 07/13/2022
## revision: 12/14/2022
## 
## chris-kreitzer


library(Rsamtools)
library(data.table)


Masterfile = read.csv('~/Master/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')
Masterfile = Masterfile[!is.na(Masterfile$Normal_bam), ]
bed = data.table::fread('~/Master/Y_features.bed')
colnames(bed) = c('chrom', 'start', 'end', 'feature')

key_cols = c('chrom', 'start', 'end')

alignment_Y = function(file){
  try({
    name = substr(x = basename(file), start = 17, stop = 51)
    print(name)
    countmatrix = data.table::fread(file, sep = ',')
    countmatrix = countmatrix[which(countmatrix$Chromosome == 'Y'), ]
    countY = data.frame(chrom = countmatrix$Chromosome,
                        start = countmatrix$Position,
                        end = countmatrix$Position,
                        N_depth = countmatrix$File1R + countmatrix$File1A,
                        T_depth = countmatrix$File2A + countmatrix$File2R)
    countY = data.table::as.data.table(countY)
    data.table::setkeyv(bed, key_cols)
    data.table::setkeyv(countY, key_cols)
    
    overlap = data.table::foverlaps(countY, bed,
                                    nomatch = 0,
                                    by.x = key_cols,
                                    by.y = key_cols)
    overlap$id = name
    overlap
  })
}

coverage = lapply(unique(Masterfile$counts_file), function(x) alignment_Y(file = x))
coverage = data.table::rbindlist(coverage)
write.table(coverage, file = '~/Master/coverage.txt', sep = '\t', row.names = F, quote = F)

message('Part2: Starting')
##----------------+
## part 2:
##----------------+
alignment_summary = function(id){
  try({
    all_out = data.frame()
    name = id
    id_data = coverage[which(coverage$id == name), ]
    
    for(i in unique(id_data$feature)){
      gene = i
      N_cov = ifelse(length(id_data$chrom[which(id_data$feature == i)]) < 20,
                     median(id_data$N_depth[which(id_data$feature == i)], na.rm = T), 
                     mean(id_data$N_depth[which(id_data$feature == i)], na.rm = T))
      T_cov = ifelse(length(id_data$chrom[which(id_data$feature == i)]) < 20,
                     median(id_data$T_depth[which(id_data$feature == i)], na.rm = T), 
                     mean(id_data$T_depth[which(id_data$feature == i)], na.rm = T))
      
      out = data.frame(gene = gene,
                       N_cov = N_cov,
                       T_cov = T_cov)
      all_out = rbind(all_out, out)
    }
    all_out$id = name
    all_out
  })
}

coverage_summary = lapply(unique(coverage$id), function(x) alignment_summary(id = x))
coverage_summary = data.table::rbindlist(coverage_summary)
write.table(coverage_summary, file = '~/Master/coverage_summary.txt', sep = '\t', row.names = F, quote = F)

message('Part3: starting')
##----------------+
## part3: BAM quality
##----------------+
GOI = bed[, c('chrom', 'start', 'end', 'feature')]

what = c('mapq', 'mrnm', 'isize')
flag = scanBamFlag(isUnmappedQuery = FALSE,
                   isNotPassingQualityControls = FALSE,
                   isSecondaryAlignment = FALSE,
                   isDuplicate = FALSE)

bam_quality_check = function(file){
  print(file)
  bam = file
  bai = paste0(substr(x = file, start = 1, stop = nchar(file) - 1), 'i')
  bamfile = BamFile(file = bam, index = bai)
  name = basename(file)
  Quality_out = data.frame()
  for(i in 1:nrow(GOI)){
    gene = GOI$feature[i]
    df = data.frame(chr = 'Y',
                    start = GOI$start[i],
                    end = GOI$end[i])
    which = makeGRangesFromDataFrame(df)
    x = scanBam(bamfile, 
                param = ScanBamParam(which = which,
                                     what = what,
                                     flag = flag))
    mapq = x[[1]]$mapq
    mate = x[[1]]$mrnm
    mate_Y = table(mate)[names(table(mate)) == 'Y'][[1]]
    mate_nonY = sum(table(mate)) - mate_Y
    out = data.frame(name = name,
                     gene = gene,
                     mapq = mean(mapq, na.rm = T),
                     mate_Y = mate_Y,
                     mate_nonY = mate_nonY)
    
    Quality_out = rbind(Quality_out, out)
    
  }
  rm(bam, bai, bamfile, mapq, mate, out)
  return(Quality_out)
}

MAPQ_BWA = lapply(unique(Masterfile$Normal_bam), function(x) bam_quality_check(file = x))
MAPQ_BWA = data.table::rbindlist(MAPQ_BWA)

write.table(MAPQ_BWA, file = '~/Master/MAPQ_BWA.txt', sep = '\t', row.names = F, quote = F)

#' out