## Unique sequence coverage of the Y-chromosome
## 
## Assess the sequence read quantity of low and high mapping quality
## Further assess the coverage at bp resolution; specifically for the Y-chromosome
## Only high-quality, non-duplicated reads with MAPQ > 45 are considered
## 
## start: 07/05/2022
## chris-kreitzer


#clean()
#gc()
#.rs.restartR()
#setwd('~/Documents/MSKCC/10_MasterThesis/')


library(Rsamtools)
system(command = paste0('module load samtools/1.9'))


#' /juno/res/dmpcollab/dmprequest/12-245/key.txt
#' /juno/res/dmpcollab/dmpshare/share/irb12_245/J/J/JJ958363-T


## Data setup;
Regions = read.csv('~/Master/Data/Y_features_all_hg38.txt', sep = '\t')
Regions = as.data.frame(Regions)
files = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/ProstateCohort_cbio.tsv', sep = '\t')
files = files[, c('Patient.ID', 'Sample.ID')]
files$array = NA
for(i in 1:nrow(files)){
  files$array[i] = strsplit(files$Sample.ID[i], split = '-')[[1]][4]
}

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

files = files[!is.na(files$tumor_bam), ]
write.table(x = files, file = '~/Master/Data/paths.txt', sep = '\t')



##-----------------
## ScanBam: filter; unfilter
all_out = data.frame()
for(i in 1:nrow(files)){
  try({
    dir.create(path = '~/tmp')
    system(command = paste0('cp ', files$tumor_bam[i], '*', ' ~/tmp'))
    bam = as.character(list.files(path = '~/tmp', pattern = paste0(basename(files$tumor_bam[i]), '.bam'), full.names = T))
    bai = as.character(list.files(path = '~/tmp', pattern = paste0(basename(files$tumor_bam[i]), '.bai'), full.names = T))
    
    #' assess .bam files
    bamfile = BamFile(file = bam, 
                      index = bai)
    what = scanBamWhat()
    
    ## which param: regions for this composite_mutation as GRangesList
    Y_interval = data.frame(chr = 'Y',
                            start = Regions$Start_hg19,
                            end = Regions$End_hg19,
                            name = Regions$hgnc_symbol)
    which = makeGRangesFromDataFrame(Y_interval)
    flag = scanBamFlag(isUnmappedQuery = FALSE,
                       isNotPassingQualityControls = FALSE,
                       isSecondaryAlignment = FALSE,
                       isDuplicate = NA)
    
    what = c("qname", "seq", 'mapq', "pos", "cigar","qwidth", "rname", "mrnm", "mpos", "mate_status")
    
    ## raw seq. read count
    param_unfiltered = ScanBamParam(which = which, 
                                    what = what)
    count_unfiltered = countBam(bamfile, param = param_unfiltered)
    count_unfiltered$tag = 'unfiltered'
    count_unfiltered$SampleID = files$Sample.ID[i]
    count_unfiltered$type = 'Tumor'
    
    ## filtered seq. read count
    param_filtered = ScanBamParam(which = which, 
                                  mapqFilter = 60,
                                  flag = flag, 
                                  what = what)
    count_filtered = countBam(bamfile, param = param_filtered)
    count_filtered$tag = 'filtered'
    count_filtered$SampleID = files$Sample.ID[i]
    count_filtered$type = 'Tumor'
    
    count_reads = rbind(count_unfiltered, count_filtered)
    print(head(count_reads))
    all_out = rbind(all_out, count_reads)
    
    #' remove bam files in temporary directory
    system(command = paste0('rm -r',' ~/tmp'))
  })
}

write.table(all_out, file = '~/Master/Data/seq.reads_filter.txt', sep = '\t', row.names = F, quote = F)


##-----------------
## filterBam; readCount; samtools 
files = read.csv('~/Master/Data/paths.txt', sep = '\t')
files = files[1:2, ]

countsOut = data.frame()
normal_countsOut = data.frame()
for(i in 1:nrow(files)){
  try({
    dir.create(path = '~/tmp')
    system(command = paste0('cp ', files$tumor_bam[i], '*', ' ~/tmp'))
    bam = as.character(list.files(path = '~/tmp', pattern = paste0(basename(files$tumor_bam[i]), '.bam'), full.names = T))
    bai = as.character(list.files(path = '~/tmp', pattern = paste0(basename(files$tumor_bam[i]), '.bai'), full.names = T))
    
    #' assess .bam files
    bamfile = BamFile(file = bam, 
                      index = bai)
    what = scanBamWhat()
    
    #' filter-rules:
    filter = FilterRules(list(chrom = function(x) x$rname == 'Y',
                              qual = function(x) x$mapq >= 45,
                              partner = function(x) x$mrnm == 'Y'))
    Y_interval = data.frame(chr = 'Y',
                            start = Regions$Start_hg19,
                            end = Regions$End_hg19,
                            name = Regions$hgnc_symbol)
    which = makeGRangesFromDataFrame(Y_interval)
    flag = scanBamFlag(isUnmappedQuery = FALSE,
                       isNotPassingQualityControls = FALSE,
                       isSecondaryAlignment = FALSE,
                       isDuplicate = NA)
    param_filtered = ScanBamParam(which = which, 
                                  mapqFilter = 60,
                                  flag = flag, 
                                  what = what)
    filterBam(file = bam, 
              destination = paste0('~/tmp/', basename(files$tumor_bam)[i], '_filtered.bam'), 
              index = bai, 
              par = param_filtered, 
              filter = filter)
    
    counts = system(command = paste0('module load samtools/1.9','\n', 'samtools depth -Q 45', ' ~/tmp/', basename(files$tumor_bam[i]), '_filtered.bam'), intern = T)
    counts = strsplit(counts, split = '\t')
    counts = data.frame(do.call(rbind, counts))
    colnames(counts) = c('Chromosome', 'Position', 'depth')
    counts$SampleID = files$Sample.ID[i]
    counts$bamID = files$tumor_bam[i]
    counts = counts[which(counts$depth >= 25), ]
    countsOut = rbind(countsOut, counts)
    
    ## concentrate on Normal samples in here
    system(command = paste0('cp ', files$normal_bam[i], '*', ' ~/tmp'))
    bam = as.character(list.files(path = '~/tmp', pattern = paste0(basename(files$normal_bam[i]), '.bam'), full.names = T))
    bai = as.character(list.files(path = '~/tmp', pattern = paste0(basename(files$normal_bam[i]), '.bai'), full.names = T))
    
    #' assess .bam files
    bamfile = BamFile(file = bam, 
                      index = bai)
    what = scanBamWhat()
    filter = FilterRules(list(chrom = function(x) x$rname == 'Y',
                              qual = function(x) x$mapq >= 45,
                              partner = function(x) x$mrnm == 'Y'))
    Y_interval = data.frame(chr = 'Y',
                            start = Regions$Start_hg19,
                            end = Regions$End_hg19,
                            name = Regions$hgnc_symbol)
    which = makeGRangesFromDataFrame(Y_interval)
    flag = scanBamFlag(isUnmappedQuery = FALSE,
                       isNotPassingQualityControls = FALSE,
                       isSecondaryAlignment = FALSE,
                       isDuplicate = NA)
    param_filtered = ScanBamParam(which = which, 
                                  mapqFilter = 60,
                                  flag = flag, 
                                  what = what)
    filterBam(file = bam, 
              destination = paste0('~/tmp/', basename(files$normal_bam)[i], '_filtered.bam'), 
              index = bai, 
              par = param_filtered, 
              filter = filter)
    
    normal_counts = system(command = paste0('module load samtools/1.9','\n', 'samtools depth -Q 45', ' ~/tmp/', basename(files$normal_bam[i]), '_filtered.bam'), intern = T)
    normal_counts = strsplit(normal_counts, split = '\t')
    normal_counts = data.frame(do.call(rbind, normal_counts))
    colnames(normal_counts) = c('Chromosome', 'Position', 'depth')
    normal_counts$SampleID = files$Sample.ID[i]
    normal_counts$bamID = files$normal_bam[i]
    normal_counts = normal_counts[which(normal_counts$depth >= 25), ]
    normal_countsOut = rbind(normal_countsOut, normal_counts)
    
    #' remove tmp directory
    system(command = paste0('rm -r',' ~/tmp'))
  })
}

write.table(countsOut, file = '~/Master/Data/Tumor_countsOut.txt', sep = '\t', row.names = F)
write.table(normal_countsOut, file = '~/Master/Data/Normal_countsOut.txt', sep = '\t', row.names = F)


# .unlist <- function (x){
#   #do.call(c)  #'coerces factor to integer, which is undesired
#   x1 <- x[[1L]]
#   if (is.factor(x1)) {
#     structure(unlist(x), class = "factor", levels = levels(x1))
#     } else {
#       do.call(c, x)
#     }
#   }
# 

