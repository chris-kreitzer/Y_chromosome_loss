##----------------+
## MADSEQ masterfile
##----------------+

source('~/MADSEQ/MADSEQ_helper.R')
Masterfile = read.csv(file = '~/MADSEQ/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')


#' Path to file
normal_bam = Masterfile$Normal_bam[2]
target = ifelse(Masterfile$GENE_PANEL[2] == 'IMPACT341',
                '~/MADSEQ/IMPACT341_b37_baits.bed', '~/MADSEQ/IMPACT410_b37_baits.bed')
genome_assembly = 'hg19'


#' get coverage
target_gr = getCoverage(bam = normal_bam, 
                        target_bed = target, 
                        genome_assembly = 'hg19')

target_gr = addchr(target_gr)
gc = calculateGC(range = target_gr)

if(length(gc) == length(target_gr)){
  mcols(target_gr)$GC = gc
}

## filter out regions around gaps
target_gr_deGAP = removeGap(target_gr,genome_assembly) 

## filter out HLA regions
target_gr_deHLA = removeHLA(target_gr_deGAP,genome_assembly)

## filter out other highly polymorphic regions
target_gr_deAQP = removeAQP(target_gr_deHLA,genome_assembly)

## filter out regions overlap with repeats
#target_gr_final = removeRE(target_gr_deAQP,genome_assembly)
target_gr_final = target_gr_deAQP
target_gr_final

normalizeCoverage(object = target_gr_final, control = NULL, writeToFile = TRUE, destination = '.', plot = FALSE)
