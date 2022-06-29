###############################################################################
#' Assess the amount of READS that uniquely map to specific genes on the Y-chromosome
#' Which genes to consider for Y-chromosome loss study
#' Inspired by Ed Reznik (mutation-phasing) script/paper
## 
## start: 06/21/2022
## revision: 06/28/2022
## chris-kreitzer


#' /juno/res/dmpcollab/dmprequest/12-245/key.txt
#' /juno/res/dmpcollab/dmpshare/share/irb12_245/J/J/JJ958363-T

Regions = read.csv('~/Master/Data/proteinCodingY.txt', sep = '\t', skip = 1)
files = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/ProstateCohort_cbio.tsv', sep = '\t')
files = files[, c('Patient.ID', 'Sample.ID')]
files$array = strsplit(files$Sample.ID, split = '-')[[1]][4]
files$tumor_bam = NA
files$normal_bam = NA
for(i in 1:nrow(files)){
  try({
    location = as.character(system(command = paste0('grep ', files$Patient.ID[i], ' /juno/res/dmpcollab/dmprequest/12-245/key.txt'), intern = T))
    location_sample = location[grep(files$array[i], location)]
    #' tumor specific
    location_tumor = grep('-T,.*$', location_sample, value = T)
    location_tumor = strsplit(location_tumor, split = ',')[[1]][2]
    location_normal = grep('-N,.*$', location_sample, value = T)
    location_normal = strsplit(location_normal, split = ',')[[1]][2]
    files$tumor_bam[i] = paste0('/juno/res/dmpcollab/dmpshare/share/irb12_245/', substr(location_tumor, start = 0, stop = 1),
                                '/', substr(location_tumor, start = 2, stop = 2), '/', location_tumor)
    files$normal_bam[i] = paste0('/juno/res/dmpcollab/dmpshare/share/irb12_245/', substr(location_normal, start = 0, stop = 1),
                                 '/', substr(location_normal, start = 2, stop = 2), '/', location_normal)
  })
}

write.table(x = files, file = '~/Master/Data/paths.txt', sep = '\t')


#' copy the files to temporary location and run RSamtools
all_out = data.frame()
for(i in 1:nrow(files)){
  try({
    dir.create(path = '~/tmp')
    system(command = paste0('cp ', files$tumor_bam[1], '*', ' ~/tmp'))
    bam = list.files(path = '~/tmp/', pattern = '.bam', full.names = T)
    bai = list.files(path = '~/tmp/', pattern = '.bai', full.names = T)
    
    #' assess .bam files
    bamfile = BamFile(file = bam, 
                      index = bai)
    what = scanBamWhat()
    
    ## which param: regions for this composite_mutation as GRangesList
    Y_interval = data.frame(chr = 'Y',
                            start = Regions$Start,
                            end = Regions$End,
                            name = Regions$Gene_name)
    
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
    all_out = rbind(all_out, count_reads)
    
    #' remove bam files in temporary directory
    system(command = paste0('rm -r ',' ~/tmp'))
  })
}

write.table(all_out, file = '~/Master/Data/AlignmentStatsY.txt', sep = '\t', quote = F)



##-----------------------------------------------------------------------------
## Assessing Read Count Information per base
library(megadepth)

#' convert Tumor sample to BigWig file
megadepth::bam_to_bigwig(bam_file = '~/Documents/MSKCC/10_MasterThesis/Data/P-0000140-T01-IM3/EM634900-T.bam', 
                         prefix = '~/Documents/MSKCC/10_MasterThesis/Data/P-0000140-T01-IM3/EM', 
                         min_unique_qual = 60, 
                         double_count = F)

#' make sure that chromosome names are equal in BED and BigWig file
BED = data.frame(chr = c(7, 17, 'X', 10),
                 start = c(55086710, 37844347, 47004620, 89623382),
                 end = c(55279321, 37884911, 47046212, 89731687),
                 name = c('EGFR', 'ERBB2', 'RBM10', 'PTEN'))

#' make sure there are no quotes or colnames in the output .bed file
#' check the bed file structure with
#' @example bed = rtracklayer::import('4Genes.bed')
write.table(BED, file = 'genes.bed', sep = '\t', quote = F, row.names = F, col.names = F)

#' Compute the coverage of selected regions
bw_cov = get_coverage(bigwig_file = 'DE840153-T.all.bw', op = 'mean', annotation = 'genes.bed')


## bamCoverage() example
## works online on command line
bamCoverage -b YB324274-N.bam -o Normal_PTEN.bw -of 'bigwig' -bs 1 -r 10:89623382:89731687 --minMappingQuality 30

#' out
