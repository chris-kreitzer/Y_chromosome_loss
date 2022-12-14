##----------------+
## Assess fraction of seq. 
## reads uniquely mapping to 
## chromosome Y features
##----------------+
##
## Inspired by Ed Reznik (mutation-phasing) script/paper
## 
## start: 06/21/2022
## revision: 06/28/2022
## revision: 07/13/2022
## revision: 12/14/2022
## 
## chris-kreitzer


Regions = as.data.frame(readxl::read_excel('~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/Y_chromosome_annotation_hg19.xlsx'))
Masterfile = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/00_CohortData/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')


#' copy .bams to /tmp and run function
#' make sure the Masterfile contains
#' normal_bam[] column
alignment_Y = function(files){
  try({
    dir.create(path = '~/tmp')
    system(command = paste0('cp ', files$tumor_bam[files], '*', ' ~/tmp'))
    bam = list.files(path = '~/tmp/', pattern = '.bam', full.names = T)
    bai = list.files(path = '~/tmp/', pattern = '.bai', full.names = T)
    
    #' assess .bam files
    bamfile = BamFile(file = bam, 
                      index = bai)
    what = scanBamWhat()
    
    ## which param: regions for this composite_mutation as GRangesList
    Y_interval = data.frame(chr = 'Y',
                            start = Regions$Start_hg19,
                            end = Regions$End_hg19,
                            name = Regions$Gene_name)
    Y_interval$end = as.integer(as.character(Y_interval$end))
    Y_interval = Y_interval[!is.na(Y_interval$end), ]
    which = makeGRangesFromDataFrame(Y_interval)
    flag = scanBamFlag(isUnmappedQuery = FALSE,
                       isNotPassingQualityControls = FALSE,
                       isSecondaryAlignment = FALSE,
                       isDuplicate = NA)
    
    what = c("qname","seq","pos","cigar","qwidth","rname","mrnm","mpos","mate_status")
    param_unfiltered = ScanBamParam(which = which, 
                                    what = what)
    param_filtered = ScanBamParam(which = which, 
                                  mapqFilter = 60,
                                  flag = flag, 
                                  what = what)
    count_unfiltered = countBam(bamfile, param = param_unfiltered)
    count_unfiltered$tag = 'unfiltered'
    count_filtered = countBam(bamfile, param = param_filtered)
    count_filtered$tag = 'filtered'
    
    count_reads = rbind(count_unfiltered, count_filtered)
    count_reads
    
    #' remove bam files in temporary directory
    system(command = paste0('rm -r',' ~/tmp'))
  })
}

out = lapply(unique(Masterfile$ID), function(x) alignment_Y(files = x))


##---------------TEST-Normal
bamfile = BamFile(file = '~/Desktop/KT218731-N.bam', 
                  index = '~/Desktop/KT218731-N.bai')
what = scanBamWhat()

## which param: regions for this composite_mutation as GRangesList
Y_interval = data.frame(chr = 'Y',
                        start = Regions$Start_hg19,
                        end = Regions$End_hg19,
                        name = Regions$Gene_hgnc_symbol)
Y_interval$end = as.integer(as.character(Y_interval$end))
Y_interval = Y_interval[!is.na(Y_interval$end), ]
which = makeGRangesFromDataFrame(Y_interval)
flag = scanBamFlag(isUnmappedQuery = FALSE,
                   isNotPassingQualityControls = FALSE,
                   isSecondaryAlignment = FALSE,
                   isDuplicate = NA)
what = c("qname","seq","pos","cigar","qwidth","rname","mrnm","mpos","mate_status")
param_unfiltered = ScanBamParam(which = which, 
                                what = what)
param_filtered = ScanBamParam(which = which, 
                              mapqFilter = 60,
                              flag = flag, 
                              what = what)
count_unfiltered = countBam(bamfile, param = param_unfiltered)
count_unfiltered$tag = 'unfiltered'
count_filtered = countBam(bamfile, param = param_filtered)
count_filtered$tag = 'filtered'





#' out


##----------------+
## Assessing read count 
## information per base
##----------------+
#' library(megadepth)
#' 
#' #' convert Tumor sample to BigWig file
#' megadepth::bam_to_bigwig(bam_file = '~/Documents/MSKCC/10_MasterThesis/Data/P-0000140-T01-IM3/EM634900-T.bam', 
#'                          prefix = '~/Documents/MSKCC/10_MasterThesis/Data/P-0000140-T01-IM3/EM', 
#'                          min_unique_qual = 60, 
#'                          double_count = F)
#' 
#' #' make sure that chromosome names are equal in BED and BigWig file
#' BED = data.frame(chr = c(7, 17, 'X', 10),
#'                  start = c(55086710, 37844347, 47004620, 89623382),
#'                  end = c(55279321, 37884911, 47046212, 89731687),
#'                  name = c('EGFR', 'ERBB2', 'RBM10', 'PTEN'))
#' 
#' #' make sure there are no quotes or colnames in the output .bed file
#' #' check the bed file structure with
#' #' @example bed = rtracklayer::import('4Genes.bed')
#' write.table(BED, file = 'genes.bed', sep = '\t', quote = F, row.names = F, col.names = F)
#' 
#' #' Compute the coverage of selected regions
#' bw_cov = get_coverage(bigwig_file = 'DE840153-T.all.bw', op = 'mean', annotation = 'genes.bed')
#' 
#' 
#' ## bamCoverage() example
#' bamCoverage -b YB324274-N.bam -o Normal_PTEN.bw -of 'bigwig' -bs 1 -r 10:89623382:89731687 --minMappingQuality 30
#' 
#' library()
#' bamCoverage -b '~/Documents/MSKCC/10_MasterThesis/Data/P-0000140-T01-IM3/EM634900-T.bam' -o test.bw -of 'bigwig' -bs 1 -r Y:89623382:89731687 --minMappingQuality 30
