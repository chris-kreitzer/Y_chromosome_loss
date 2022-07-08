## FacetsY on IMPACT  samples;
## This script is intended to be run on the juno
##
## start: 04/16/2021
## revision: 08/25/2021
## revision: 07/08/2022
## chris kreitzer


## install local FacetsY and pctGCdata
require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
source('/juno/home/kreitzec/Y_chromosome_loss/Loss/FacetsQC.R')
library(facetsSuite)
library(dplyr)
library(data.table)


# load IMPACT samplepaths
IMPACT.samples = readRDS('/juno/home/kreitzec/Y_chromosome_loss/Data/cohort_data.rds')
IMPACT.samples = as.character(IMPACT.samples$IMPACT.cohort$counts_file)

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



cval.purity = 300
cval.postprocess = 100


CopyNumber_out = data.frame()
flags.out = c()
Cnlr_out = data.frame()
IMPACT_IGV_out = data.frame()
IMPACT_arm_change_out = data.frame()
QC_metrics = data.frame()

for(i in unique(IMPACT.samples)){
  try({
    print(i)
    data.in = facetsY::readSnpMatrix(i, err.thresh = 10, del.thresh = 10)
    
    #' split allosomes and autosomes
    Y_chromosome = data.in[which(data.in$Chromosome == 'Y' & 
                                   data.in$Position >= 2654550 & 
                                   data.in$Position <= 28000000), ]
    
    #' exclude PCDH11Y and centromeric region
    Y_chromosome = Y_chromosome[with(Y_chromosome, !((Position %in% 4922131:5612269))), ]
    Y_chromosome = Y_chromosome[with(Y_chromosome, !((Position %in% 10500000:14000000))), ]
    
    
    #--------------
    # Segmentation on Y:
    set.seed(100)
    Y_chromo_segmentation = facetsY::preProcSample(rcmat = Y_chromosome,
                                                   ndepth = 20,
                                                   snp.nbhd = 50,
                                                   gbuild = 'hg19')
    
    pmat_Y = Y_chromo_segmentation$pmat
    joint_Y = Y_chromo_segmentation$jointseg
    seg_Y = Y_chromo_segmentation$seg.tree
    
    
    ##-------------
    ## autosomes:
    autosomes = data.in[which(data.in$Chromosome %in% c(seq(1, 22, 1), 'X')), ]
    set.seed(100)
    auto.segmentation = facetsY::preProcSample(rcmat = autosomes,
                                               ndepth = 35,
                                               snp.nbhd = 150, 
                                               gbuild = 'hg19')
    
    pmat_auto = auto.segmentation$pmat
    joint_auto = auto.segmentation$jointseg
    seg_auto = auto.segmentation$seg.tree
    
    
    ##-------------
    ## Merge and run together
    pmat = rbind(pmat_Y, pmat_auto)
    joint = rbind(joint_Y, joint_auto)
    seg = c(seg_auto, seg_Y)
    attr(x = seg, which = 'cval') = 25
    
    facets_pre = list(pmat = pmat,
                      gbuild = 'hg19',
                      nX = 23,
                      seg.tree = seg,
                      jointseg = joint,
                      hscl = auto.segmentation$hscl,
                      chromlevels = seq(1,24, 1))
    

    #' first run FacetsY on wider cval (purity) to determine dipLogR
    data.process = facetsY::procSample(x = facets_pre, 
                                       cval = cval.purity)
    data.out = facetsY::emcncf(data.process)
    
    #' run fine-tuning segmentation
    dipLogR.purity = data.out$dipLogR
    data.process_out = facetsY::procSample(facets_pre,
                                           cval = cval.postprocess,
                                           dipLogR = dipLogR.purity)
    data.out = facetsY::emcncf(data.process_out)
    
    ID = substr(i, start = 59, stop = 75)
    purity = data.out$purity
    purity = ifelse(is.na(purity), 0, purity)
    
    parameter.selected = paste0('cval_purity=', cval.purity, '/', cval.postprocess, '::snp.nbhd=150|50')
    data.return = data.out$cncf
    data.return$ID = ID
    data.return$parameter = parameter.selected
    data.return$purity = purity
    data.return$ploidy = data.out$ploidy
    data.return$XY_ratio = median(data.process_out$jointseg$cnlr[which(data.process_out$jointseg$chrom == 23)], na.rm = T) / 
      median(data.process_out$jointseg$cnlr[which(data.process_out$jointseg$chrom == 24)], na.rm = T)
    
    #' run quality metric
    FacetsQuality = facetsSuite::run_facets(data.in,
                                            genome = 'hg19',
                                            cval = cval.postprocess, 
                                            dipLogR = dipLogR.purity, 
                                            snp_nbhd = 150, 
                                            seed = 100)
    
    Quality = facetsSuite::check_fit(FacetsQuality)
    fit.out = facets_fit_qc(facets_output = FacetsQuality)
    fit.out_df = data.frame(QC = fit.out$facets_qc,
                            wgd = Quality$wgd,
                            ID = ID)
    
    QC_metrics = rbind(QC_metrics, fit.out_df)
    
    # catch flags, if there are some
    facets.flags = data.out$emflags
    if(!is.null(facets.flags)){
      flags.out = c(flags.out, paste0(facets.flags, ';', ID))
    }
    
    CopyNumber_out = rbind(CopyNumber_out, data.return)
    
    # fetch cnlr values for chrY; orthogonal validation
    data.cnlr = data.process_out$jointseg
    data.cnlr = data.cnlr[which(data.cnlr$chrom == 24), ]
    data.cnlr$ID = ID
    Cnlr_out = rbind(Cnlr_out, data.cnlr)
    
    #' run FacetsY to IGV converstion
    #' create list
    FacetsY_output = list(snps = data.process_out$jointseg,
                          segs = data.out$cncf,
                          dipLogR = data.out$dipLogR)
    
    IGV_out = format_igv_seg(facets_output = FacetsY_output,
                             sample_id = ID)
    
    IMPACT_IGV_out = rbind(IMPACT_IGV_out, IGV_out)
    
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
    
    IMPACT_arm_change_out = rbind(IMPACT_arm_change_out , arm_change.out)
    
    
  })
}

write.table(x = CopyNumber_out, file = '/juno/home/kreitzec/Y_chromosome_loss/Loss/IMPACT_copynumber_out.txt', row.names = F, sep = '\t')
write.table(flags.out, file = '/juno/home/kreitzec/Y_chromosome_loss/Loss/IMPACT_errorflags', row.names = F)
write.table(Cnlr_out, file = '/juno/home/kreitzec/Y_chromosome_loss/Loss/Cnlr_out.txt', row.names = F, sep = '\t')
write.table(IMPACT_IGV_out, file = '/juno/home/kreitzec/Y_chromosome_loss/Loss/IMPACT_IGV_out.seg', row.names = F, sep = '\t', quote = F)
write.table(IMPACT_arm_change_out, file = '/juno/home/kreitzec/Y_chromosome_loss/Loss/IMPACT_arm_change_out.txt', row.names = F, sep = '\t')
write.table(QC_metrics, file = '/juno/home/kreitzec/Y_chromosome_loss/Loss/QC_metrics.txt', row.names = F, sep = '\t')


#' qualitatively call y chromosome loss
## qualitative calling of WES loss (50% rule) and adherent analysis
## Functions:
#' calling Y_chromosome loss from Facets output
message('evaluating binary Y loss')

Y_loss_call = function(data, sample.id){
  data.sub = data[which(data$ID == sample.id), ]
  data.sub = data.sub[which(data.sub$chrom == 24), ]
  if(nrow(data.sub) == 0){
    return(data.frame(cbind(Y_call = NA, 
                            sample.id,
                            length = NA,
                            ploidy = NA,
                            purity = NA,
                            min.segment = NA,
                            max.segment = NA)))
  } else {
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
    
    return(data.frame(cbind(Y_call, 
                            sample.id,
                            length,
                            ploidy,
                            purity,
                            min.segment = min.segment,
                            max.segment = max.segment)))
  }
}
## apply to data, analyzed above
IMPACT.binaryLoss_out = lapply(unique(CopyNumber_out$ID),
                               function(x) Y_loss_call(data = CopyNumber_out, sample.id = x))

IMPACT.binaryLoss_out = data.frame(do.call('rbind', IMPACT.binaryLoss_out))
write.table(IMPACT.binaryLoss_out, 
            file = '/juno/home/kreitzec/Y_chromosome_loss/Loss/IMPACT.binaryLoss_out.txt', 
            row.names = F,
            sep = '\t')


## out