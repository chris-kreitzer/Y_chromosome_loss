##----------------+
## Get autosomal coverage
## and chromosome Y coverage
## calculate the ratio;
## hints for mosaic loss of Y
##----------------+

## start: 12/04/2022
## chris-kreitzer



##----------------+
## MADSEQ masterfile
##----------------+

source('~/MADSEQ/MADSEQ_helper.R')
Masterfile = read.csv(file = '~/MADSEQ/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')
Masterfile = Masterfile[!is.na(Masterfile$Normal_bam), ]
Masterfile = Masterfile[1:2, ]

#' Path to file
mLOY = function(bams){
  try({
    name = Masterfile$SAMPLE_ID[bams]
    normal_bam = Masterfile$Normal_bam[bams]
    target = ifelse(Masterfile$GENE_PANEL[bams] == 'IMPACT341',
                    '~/MADSEQ/IMPACT341_b37_baits.bed',
                    ifelse(Masterfile$GENE_PANEL[bams] == 'IMPACT410',
                           '~/MADSEQ/IMPACT410_b37_baits.bed',
                           ifelse(Masterfile$GENE_PANEL[bams] == 'IMPACT468',
                                  '~/MADSEQ/IMPACT468_b37_baits.bed', '~/MADSEQ/IMPACT505_b37_baits.bed')))
    print(bams)
    genome_assembly = 'hg19'
    
    target_gr = getCoverage(bam = normal_bam, 
                            target_bed = target, 
                            genome_assembly = genome_assembly)
    
    target_gr = addchr(target_gr)
    gc = calculateGC(range = target_gr)
    
    if(length(gc) == length(target_gr)){
      mcols(target_gr)$GC = gc
    }
    
    target_gr_deGAP = removeGap(target_gr,genome_assembly) 
    target_gr_deHLA = removeHLA(target_gr_deGAP,genome_assembly)
    target_gr_deAQP = removeAQP(target_gr_deHLA,genome_assembly)
    target_gr_final = target_gr_deAQP
    
    normalizeCoverage(object = target_gr_final, 
                      control = NULL, 
                      writeToFile = TRUE, 
                      destination = './Normalized', 
                      sampleID = name,
                      plot = FALSE)
  })
}

sapply(1:nrow(Masterfile), function(x) mLOY(bams = x))

#' out