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
Masterfile = sample(Masterfile$Normal_bam, size = 5000, replace = F)
bed = data.table::fread('~/Master/Y_features.bed')
colnames(bed) = c('chrom', 'start', 'end', 'feature')

#key_cols = c('chrom', 'start', 'end')

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

#coverage = lapply(unique(Masterfile$counts_file), function(x) alignment_Y(file = x))
#coverage = Filter(function(x) length(x) > 1, coverage)
#coverage = data.table::rbindlist(coverage)
#write.table(coverage, file = '~/Master/coverage.txt', sep = '\t', row.names = F, quote = F)

#message('Part2: Starting')
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

#coverage_summary = lapply(unique(coverage$id), function(x) alignment_summary(id = x))
#coverage_summary = data.table::rbindlist(coverage_summary)
#write.table(coverage_summary, file = '~/Master/coverage_summary.txt', sep = '\t', row.names = F, quote = F)

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
  try({
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
    })
  rm(bam, bai, bamfile, mapq, mate, out)
  return(Quality_out)
}

print(head(Masterfile))
MAPQ_BWA = lapply(unique(Masterfile), function(x) bam_quality_check(file = x))
MAPQ_BWA = Filter(function(x) length(x) > 1, xx)
MAPQ_BWA = data.table::rbindlist(MAPQ_BWA)

write.table(MAPQ_BWA, file = '~/Master/MAPQ_BWA.txt', sep = '\t', row.names = F, quote = F)






#' #' copy .bams to /tmp and run function
#' #' make sure the Masterfile contains
#' #' normal_bam[] column
#' alignment_Y = function(files){
#'   try({
#'     dir.create(path = '~/tmp')
#'     system(command = paste0('cp ', files$tumor_bam[files], '*', ' ~/tmp'))
#'     bam = list.files(path = '~/tmp/', pattern = '.bam', full.names = T)
#'     bai = list.files(path = '~/tmp/', pattern = '.bai', full.names = T)
#'     
#'     #' assess .bam files
#'     bamfile = BamFile(file = bam, 
#'                       index = bai)
#'     what = scanBamWhat()
#'     
#'     ## which param: regions for this composite_mutation as GRangesList
#'     Y_interval = data.frame(chr = 'Y',
#'                             start = Regions$Start_hg19,
#'                             end = Regions$End_hg19,
#'                             name = Regions$Gene_name)
#'     Y_interval$end = as.integer(as.character(Y_interval$end))
#'     Y_interval = Y_interval[!is.na(Y_interval$end), ]
#'     which = makeGRangesFromDataFrame(Y_interval)
#'     flag = scanBamFlag(isUnmappedQuery = FALSE,
#'                        isNotPassingQualityControls = FALSE,
#'                        isSecondaryAlignment = FALSE,
#'                        isDuplicate = NA)
#'     
#'     what = c("qname","seq","pos","cigar","qwidth","rname","mrnm","mpos","mate_status")
#'     param_unfiltered = ScanBamParam(which = which, 
#'                                     what = what)
#'     param_filtered = ScanBamParam(which = which, 
#'                                   mapqFilter = 60,
#'                                   flag = flag, 
#'                                   what = what)
#'     count_unfiltered = countBam(bamfile, param = param_unfiltered)
#'     count_unfiltered$tag = 'unfiltered'
#'     count_filtered = countBam(bamfile, param = param_filtered)
#'     count_filtered$tag = 'filtered'
#'     
#'     count_reads = rbind(count_unfiltered, count_filtered)
#'     count_reads
#'     
#'     #' remove bam files in temporary directory
#'     system(command = paste0('rm -r',' ~/tmp'))
#'   })
#' }
#' 
#' out = lapply(unique(Masterfile$ID), function(x) alignment_Y(files = x))
#' 
#' 
#' ##---------------TEST-Normal
#' bamfile = BamFile(file = '~/Desktop/KT218731-N.bam', 
#'                   index = '~/Desktop/KT218731-N.bai')
#' what = scanBamWhat()
#' 
#' ## which param: regions for this composite_mutation as GRangesList
#' Y_interval = data.frame(chr = 'Y',
#'                         start = Regions$Start_hg19,
#'                         end = Regions$End_hg19,
#'                         name = Regions$Gene_hgnc_symbol)
#' Y_interval$end = as.integer(as.character(Y_interval$end))
#' Y_interval = Y_interval[!is.na(Y_interval$end), ]
#' which = makeGRangesFromDataFrame(Y_interval)
#' flag = scanBamFlag(isUnmappedQuery = FALSE,
#'                    isNotPassingQualityControls = FALSE,
#'                    isSecondaryAlignment = FALSE,
#'                    isDuplicate = NA)
#' what = c("qname","seq","pos","cigar","qwidth","rname","mrnm","mpos","mate_status")
#' param_unfiltered = ScanBamParam(which = which, 
#'                                 what = what)
#' param_filtered = ScanBamParam(which = which, 
#'                               mapqFilter = 60,
#'                               flag = flag, 
#'                               what = what)
#' count_unfiltered = countBam(bamfile, param = param_unfiltered)
#' count_unfiltered$tag = 'unfiltered'
#' count_filtered = countBam(bamfile, param = param_filtered)
#' count_filtered$tag = 'filtered'
#' 
#' 
#' 
#' 
#' 
#' #' out
#' 
#' 
#' ##----------------+
#' ## Assessing read count 
#' ## information per base
#' ##----------------+
#' #' library(megadepth)
#' #' 
#' #' #' convert Tumor sample to BigWig file
#' #' megadepth::bam_to_bigwig(bam_file = '~/Documents/MSKCC/10_MasterThesis/Data/P-0000140-T01-IM3/EM634900-T.bam', 
#' #'                          prefix = '~/Documents/MSKCC/10_MasterThesis/Data/P-0000140-T01-IM3/EM', 
#' #'                          min_unique_qual = 60, 
#' #'                          double_count = F)
#' #' 
#' #' #' make sure that chromosome names are equal in BED and BigWig file
#' #' BED = data.frame(chr = c(7, 17, 'X', 10),
#' #'                  start = c(55086710, 37844347, 47004620, 89623382),
#' #'                  end = c(55279321, 37884911, 47046212, 89731687),
#' #'                  name = c('EGFR', 'ERBB2', 'RBM10', 'PTEN'))
#' #' 
#' #' #' make sure there are no quotes or colnames in the output .bed file
#' #' #' check the bed file structure with
#' #' #' @example bed = rtracklayer::import('4Genes.bed')
#' #' write.table(BED, file = 'genes.bed', sep = '\t', quote = F, row.names = F, col.names = F)
#' #' 
#' #' #' Compute the coverage of selected regions
#' #' bw_cov = get_coverage(bigwig_file = 'DE840153-T.all.bw', op = 'mean', annotation = 'genes.bed')
#' #' 
#' #' 
#' #' ## bamCoverage() example
#' #' bamCoverage -b YB324274-N.bam -o Normal_PTEN.bw -of 'bigwig' -bs 1 -r 10:89623382:89731687 --minMappingQuality 30
#' #' 
#' #' library()
#' #' bamCoverage -b '~/Documents/MSKCC/10_MasterThesis/Data/P-0000140-T01-IM3/EM634900-T.bam' -o test.bw -of 'bigwig' -bs 1 -r Y:89623382:89731687 --minMappingQuality 30
