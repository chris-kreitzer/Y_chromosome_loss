## FacetsY on WES  Prostate cancer samples; n = 238
## This script is intened to be run on the juno cluster

## run purity run (cval = 700) first; then get dipLogR and run with cval 150
## exclude base coverage > 30 Mb (see coverage plot) (PAR2 region not relevant)
## hinders Segmentation
## run FacetsQC and flag samples with QC == FALSE

## 04/14/2021

# install local FacetsY and pctGCdata
require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
source('/juno/home/kreitzec/WES_Prostate/FacetsQC.R')
library(facetsSuite)
library(dplyr)
library(data.table)

# load FACETS countsfiles for Prostate WES samples
WES.prostate.samples = read.csv('/juno/home/kreitzec/WES_Prostate/WES.prostate.samplepath.txt', header = F)
WES.prostate.samples = as.character(WES.prostate.samples$V1)

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


cval.preprocess = 25
cval.purity = 700
cval.postprocess = 150
snp.nbhd = 250

Prostate.out_df = data.frame()
flags.out = c()
Cnlr.prostate.out_df = data.frame()
WES.prostate.IGV_out = data.frame()
WES.prostate.arm_change_out = data.frame()
WES.prostate.fit.out = data.frame()

for(i in unique(WES.prostate.samples)){
  try({
    print(i)
    data.in = facetsY::readSnpMatrix(i)
    data.in = rbind(data.in[which(data.in$Chromosome != 'Y'), ], 
                    data.in[which(data.in$Chromosome == 'Y' & data.in$Position < 30000000), ])
    
    data.pre = facetsY::preProcSample(data.in, 
                                      cval = cval.preprocess, 
                                      gbuild = 'hg19', 
                                      snp.nbhd = snp.nbhd)
    
    #' first run FacetsY on wider cval (purity) to determine dipLogR
    data.process = facetsY::procSample(data.pre, 
                                       cval = cval.purity)
    data.out = facetsY::emcncf(data.process)
    
    #' run fine-tuning segmentation
    dipLogR.purity = data.out$dipLogR
    data.process_out = facetsY::procSample(data.pre,
                                           cval = cval.postprocess,
                                           dipLogR = dipLogR.purity)
    data.out = facetsY::emcncf(data.process_out)
    
    ID = substr(i, start = 70, stop = 98)
    purity = data.out$purity
    purity = ifelse(is.na(purity), 0, purity)
    
    parameter.selected = paste0('cval=', cval.preprocess, '/', cval.purity, '/', cval.postprocess, '::snp.nbhd=', snp.nbhd)
    data.return = data.out$cncf
    data.return$ID = ID
    data.return$parameter = parameter.selected
    data.return$purity = purity
    data.return$ploidy = data.out$ploidy
    
    #' run quality metric on WES samples
    FacetsQuality = facetsSuite::run_facets(data.in,
                                            genome = 'hg19',
                                            cval = cval.postprocess, 
                                            dipLogR = dipLogR.purity, 
                                            snp_nbhd = snp.nbhd)
    
    fit.out = facets_fit_qc(facets_output = FacetsQuality)
    fit.out_df = data.frame(QC = fit.out$facets_qc,
                            ID = ID)
    
    WES.prostate.fit.out = rbind(WES.prostate.fit.out, fit.out_df)
    
    # catch flags, if there are some
    
    facets.flags = data.out$emflags
    if(!is.null(facets.flags)){
      flags.out = c(flags.out, paste0(facets.flags, ';', ID))
    }
    
    Prostate.out_df = rbind(Prostate.out_df, data.return)
    
    # fetch cnlr values for chrY; orthogonal validation
    data.cnlr = data.process_out$jointseg
    data.cnlr = data.cnlr[which(data.cnlr$chrom == 24), ]
    data.cnlr$ID = ID
    Cnlr.prostate.out_df = rbind(Cnlr.prostate.out_df, data.cnlr)
    
    #' run FacetsY to IGV converstion
    #' create list
    FacetsY_output = list(snps = data.process_out$jointseg,
                          segs = data.out$cncf,
                          dipLogR = data.out$dipLogR)
    
    IGV_out = format_igv_seg(facets_output = FacetsY_output,
                             sample_id = ID)
    
    WES.prostate.IGV_out = rbind(WES.prostate.IGV_out, IGV_out)
    
    #' run facetsSuite and create arm-level-change
    arm_change = facetsSuite::arm_level_changes(segs = FacetsY_output$segs,
                                   ploidy = data.out$ploidy,
                                   genome = 'hg19',
                                   algorithm = 'em')
    
    arm_change.out = arm_change$full_output
    arm_change.out$WGD = arm_change$genome_doubled
    arm_change.out$fcna = arm_change$fraction_cna
    arm_change.out$AS = arm_change$aneuploidy_score
    arm_change.out$ID = ID
    
    WES.prostate.arm_change_out = rbind(WES.prostate.arm_change_out, arm_change.out)
    
    
  })
}

write.table(Prostate.out_df, file = '/home/kreitzec/WES_Prostate/WES.Prostate.processed.out.14.4.txt', row.names = F, sep = '\t')
write.table(flags.out, file = '/home/kreitzec/WES_Prostate/WES_errorflags.14.4.txt', row.names = F)
write.table(Cnlr.prostate.out_df, file = '/home/kreitzec/WES_Prostate/Cnlr.WES.prostate.out.14.4.txt', row.names = F, sep = '\t')
write.table(WES.prostate.IGV_out, file = '/home/kreitzec/WES_Prostate/WES.prostate.IGV.14.4.txt', row.names = F, sep = '\t')
write.table(WES.prostate.arm_change_out, file = '/home/kreitzec/WES_Prostate/WES.prostate.arm_change.14.4.txt', row.names = F, sep = '\t')
write.table(WES.prostate.fit.out, file = '/home/kreitzec/WES_Prostate/WES.prostate.QC.14.4.txt', row.names = F, sep = '\t')

#' qualitatively call y chromosome loss
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
            file = '/juno/home/kreitzec/WES_Prostate/WES.Prostate.binary.14.4.txt', 
            row.names = F,
            sep = '\t')

## out