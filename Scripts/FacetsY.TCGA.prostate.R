## FacetsY on WES-TCGA Prostate Samples:
## we have n=333 samples from Cell 2015 paper
## Confirm proper parameter selection for FacetsY (see below) 
## create normal cncf file and IGV view 
## 
## look at specific losses (gene-level); Affymetrix6.0 (gold standard)


# rm(list = ls())
# .rs.restartR()
# setwd('~/Documents/MSKCC/04_Y_chromo_loss/')

## Libraries:
library(facetsY)
library(pctGCdata)
library(dplyr)
library(data.table)

## Input:
# TCGA.PRAD.annotation_df = read.csv('Data/TCGA-PRAD-clinical.cbio.tsv', sep = '\t') #' Cell 2015 cohort
# TCGA.PRAD.paths = read.csv('Data/PRAD_summary.out', sep = '\t')


## Processing:

#' create file for TCGA-PRAD path to run FacetsY on those
# TCGA.PRAD.paths = TCGA.PRAD.paths[, c('Sample', 'path')]
# TCGA.path = strsplit(TCGA.PRAD.paths$path, split = '/')
# TCGA.path = sapply(TCGA.path, function(x) x[5])
# TCGA.PRAD.paths$short.path = TCGA.path
# 
# TCGA.PRAD.samplepath = data.frame()
# for(i in unique(TCGA.PRAD.annotation_df$Sample.ID)){
#   if(any(grepl(pattern = i, x = TCGA.PRAD.paths$short.path))){
#     count.path = grep(pattern = i, x = TCGA.PRAD.paths$short.path, value = T)
#     sampleID = i
#     out = data.frame(TCGA.ID = sampleID,
#                      TCGA.pileup = count.path)
#   } else next
#   TCGA.PRAD.samplepath = rbind(TCGA.PRAD.samplepath, out)
# }
# 
# TCGA.PRAD.samplepath = TCGA.PRAD.samplepath[!duplicated(TCGA.PRAD.samplepath$TCGA.ID), ]
# TCGA.PRAD.samplepath$TCGA.pileup = paste0('/ifs/res/taylorlab/tcga_facets/snp-pileup/', TCGA.PRAD.samplepath$TCGA.pileup)
# TCGA.PRAD.samplepath$TCGA.pileup = gsub(pattern = '(.*)_\\w+', '\\1', TCGA.PRAD.samplepath$TCGA.pileup)
# TCGA.PRAD.samplepath$TCGA.pileup = substr(TCGA.PRAD.samplepath$TCGA.pileup, start = 1, stop = 143)
# write.table(TCGA.PRAD.samplepath[, c('TCGA.pileup')], file = 'Data_out/TCGA.PRAD.samplepath.txt', 
#             col.names = F, row.names = F, sep = '\t', quote = F)


## Analysis:
# load FACETS countsfiles for Prostate WES samples
WES.prostate.samples = read.csv('/home/kreitzec/WES_Prostate/TCGA.PRAD.samplepath.txt', header = F)
WES.prostate.samples = as.character(WES.prostate.samples$V1)
print(WES.prostate.samples)

# loop through every file and run FacetsY
# rbind output into comprehensive cncf file (downstream analysis on local PC)

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


# FacetsY parameters:
cval.preprocess = 25
cval.postprocess = 150
snp.nbhd = 250

Prostate.out_df = data.frame()
flags.out = c()
Cnlr.prostate.out_df = data.frame()
WES.prostate.IGV_out = data.frame()

for(i in unique(WES.prostate.samples)){
  try({
    print(i)
    data.in = facetsY::readSnpMatrix(i)
    data.pre = facetsY::preProcSample(data.in, 
                                      cval = cval.preprocess, 
                                      gbuild = 'hg19', 
                                      snp.nbhd = snp.nbhd)
    data.process = facetsY::procSample(data.pre, 
                                       cval = cval.postprocess)
    data.out = facetsY::emcncf(data.process)
    
    ID = substr(i, start = 71, stop = 99)
    purity = data.out$purity
    purity = ifelse(is.na(purity), 0, purity)
    
    parameter.selected = paste0('cval=', cval.preprocess, '/', cval.postprocess, '::snp.nbhd=', snp.nbhd)
    data.return = data.out$cncf
    data.return$ID = ID
    data.return$parameter = parameter.selected
    data.return$purity = purity
    data.return$ploidy = data.out$ploidy
    
    # catch flags, if there are some
    
    facets.flags = data.out$emflags
    if(!is.null(facets.flags)){
      flags.out = c(flags.out, paste0(facets.flags, ';', ID))
    }
    
    Prostate.out_df = rbind(Prostate.out_df, data.return)
    
    # fetch cnlr values for chrY; orthogonal validation
    data.cnlr = data.process$jointseg
    data.cnlr = data.cnlr[which(data.cnlr$chrom == 24), ]
    data.cnlr$ID = ID
    Cnlr.prostate.out_df = rbind(Cnlr.prostate.out_df, data.cnlr)
    
    #' run FacetsY to IGV converstion
    #' create list
    FacetsY_output = list(snps = data.process$jointseg,
                          segs = data.out$cncf,
                          dipLogR = data.out$dipLogR)
    
    IGV_out = format_igv_seg(facets_output = FacetsY_output,
                             sample_id = ID)
    
    WES.prostate.IGV_out = rbind(WES.prostate.IGV_out, IGV_out)
    
  })
}

write.table(Prostate.out_df, file = '/home/kreitzec/WES_Prostate/WES.TCGA.Prostate.processed.out.txt', row.names = F, sep = '\t')
write.table(flags.out, file = '/home/kreitzec/WES_Prostate/WES.TCGA_errorflags.txt', row.names = F)
write.table(Cnlr.prostate.out_df, file = '/home/kreitzec/WES_Prostate/Cnlr.WES.TCGA.prostate.out.txt', row.names = F, sep = '\t')
write.table(WES.prostate.IGV_out, file = '/home/kreitzec/WES_Prostate/WES.TCGA.prostate.IGV.txt', row.names = F, sep = '\t')

## qualitative calling of WES loss (50% rule) and adherent analysis
## Functions:
#' calling Y_chromosome loss from Facets output
message('evaluating binary Y loss')
Y_loss_call = function(data, sample.id){
  data.sub = data[which(data$ID == sample.id), ]
  data.sub = data.sub[which(data.sub$chrom == 24), ]
  data.sub$length = data.sub$end - data.sub$start
  Y_ratio = sum(data.sub$length[which(data.sub$tcn.em == 0)], na.rm = T) / sum(data.sub$length, na.rm = T)
  
  # return values:
  purity = unique(data.sub$purity)
  length = unique(sum(data.sub$length, na.rm = T))
  ploidy = unique(data.sub$ploidy)
  min.segment = unique(min(data.sub$start))
  max.segment = unique(max(data.sub$end))
  
  Y_call = ifelse(Y_ratio >= 0.5, 'Y_chrom_loss', 'intact_Y_chrom')
  Y_call = unique(Y_call)
  
  sample.id = unique(sample.id)
  
  return(cbind(Y_call, 
               sample.id,
               length,
               ploidy,
               purity,
               min.segment = min.segment,
               max.segment = max.segment))
}

## apply to data, analyzed above
WES.Prostate.Y.out_df = lapply(unique(Prostate.out_df$ID),
                               function(x) Y_loss_call(data = Prostate.out_df, sample.id = x))

WES.Prostate.Y.out_df = data.frame(do.call('rbind', WES.Prostate.Y.out_df))
write.table(WES.Prostate.Y.out_df, 
            file = '/home/kreitzec/WES_Prostate/WES.TCGA.Prostate.binary.txt', 
            row.names = F,
            sep = '\t')

## out