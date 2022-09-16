##-----------------
## Check sequence read quality
## of chosen features; necesseary
## to avoid biased read-depth data
## 
## 09/16/2022
## chris-kreitzer
##-----------------

clean()
gc()
.rs.restartR()

library(Rsamtools)

## data setup
genes_analysis = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/Gene_Coverage_Final.txt', sep = '\t')
genes_analysis = as.character(genes_analysis$gene)
annotations = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/Y_chromosome_annotation_hg19.txt', sep = '\t')
GOI = annotations[which(annotations$Gene_hgnc_symbol %in% genes_analysis), c('Gene_hgnc_symbol', 'Chromosome', 'Start_hg19', 'End_hg19')]

##-----------------
## BAM parameters: Rsamtools
##-----------------
what = c('mapq', 'mrnm', 'isize')
flag = scanBamFlag(isUnmappedQuery = FALSE,
                   isNotPassingQualityControls = FALSE,
                   isSecondaryAlignment = FALSE,
                   isDuplicate = FALSE)

bam = '~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/P-0004539-T01-IM5/ZU367405-N.bam'
bai = '~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/P-0004539-T01-IM5/ZU367405-N.bai'
bamfile = BamFile(file = bam, 
                  index = bai)

file = c('~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/P-0004539-T01-IM5/ZU367405-N',
         '~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/P-0004539-T01-IM5/ZU367405-T')

bam_quality_check = function(file){
  bam = paste0(file, '.bam')
  bai = paste0(file, '.bai')
  bamfile = BamFile(file = bam, index = bai)
  name = basename(file)
  Quality_out = data.frame()
  for(i in 1:nrow(GOI)){
    gene = GOI$Gene_hgnc_symbol[i]
    df = data.frame(chr = 'Y',
                    start = GOI$Start_hg19[i],
                    end = GOI$End_hg19[i])
    which = makeGRangesFromDataFrame(df)
    x = scanBam(bamfile, param = ScanBamParam(which = which,
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

lapply(file, function(x) bam_quality_check(file = x))


# all_out = data.frame()
# for(gene in 1:nrow(Genes_covered_annoation)){
#   gene_name = Genes_covered_annoation$Gene_hgnc_symbol[gene]
#   #u = c()
#   #print(Genes_covered_annoation$Gene_hgnc_symbol[gene])
#   for(i in seq(Genes_covered_annoation$Start_hg19[gene] + 300, 
#                Genes_covered_annoation$End_hg19[gene] - 300, length.out = 10)){
#     print(i)
#     gene = Genes_covered_annoation$Gene_hgnc_symbol[gene]
#     df = data.frame(chr = 'Y',
#                     start = i,
#                     end = i)
#     which = makeGRangesFromDataFrame(df)
#     x = scanBam(bamfile, param = ScanBamParam(which = which,
#                                               what = what,
#                                               flag = flag))
#     mapq = x[[1]]$mapq
#     
#     if(length(mapq) != 0){
#       values = mapq
#     } else {
#       values = NA
#     }
#     
#     out = data.frame(gene = gene_name,
#                      value = values)
#     all_out = rbind(all_out, out)
#   }
# }





##-----------------
bam = '~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/P-0004539-T01-IM5/ZU367405-N.bam'
bai = '~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/P-0004539-T01-IM5/ZU367405-N.bai'
bamfile = BamFile(file = bam, 
                  index = bai)
what = c('mapq', 'mrnm', 'isize')
flag = scanBamFlag(isUnmappedQuery = FALSE,
                   isNotPassingQualityControls = FALSE,
                   isSecondaryAlignment = FALSE,
                   isDuplicate = FALSE)

df = data.frame(chromosome = 'Y',
                start = 2654895,
                end = 2655723)
which = makeGRangesFromDataFrame(df)

x = scanBam(file = bamfile, param = ScanBamParam(flag = flag, what = what, which = which))
mapq = x$`Y:13568223-13574392`$mapq
mate = x$`Y:13568223-13574392`$mrnm
ru = table(mate)
table(mate)[names(table(mate)) == 'Y'][[1]]
ru[names(ru) == 'Y'][[1]]


plot(density(mapq))
summary(mapq)
