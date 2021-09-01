## FacetsY on TCGA data; n = 4,875 male samples
## equal procedure as for IMPACT-WES samples (same parameters)
##
## 08/31/2021
## chris-kreitzer



###############################################################################
# install local FacetsY and pctGCdata
require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
source('/juno/home/kreitzec/WES_Prostate/FacetsQC.R')
library(facetsSuite)
library(dplyr)
library(data.table)

# load FACETS countsfiles for Prostate WES samples
TCGA.samples = read.csv('/juno/home/kreitzec/Y_chromosome_loss/Data/TCGA_paths.txt', sep = '\t')
TCGA.samples = as.character(TCGA.samples$path)


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
cval.postprocess = 200
snp.nbhd = 250

TCGA_copyNumber_out = data.frame()
flags.out = c()
TCGA_cnLR_out = data.frame()
TCGA_IGV_out = data.frame()
TCGA_arm_alteration_out = data.frame()
TCGA_QC_out = data.frame()

for(i in unique(TCGA.samples)){
  try({
    print(i)
    data.in = facetsY::readSnpMatrix(i, err.thresh = 10, del.thresh = 10)
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
    
    ID = substr(i, start = 39, stop = 50)
    purity = data.out$purity
    purity = ifelse(is.na(purity), 0, purity)
    
    parameter.selected = paste0('cval=', cval.preprocess, '/', cval.purity, '/', cval.postprocess, '::snp.nbhd=', snp.nbhd)
    data.return = data.out$cncf
    data.return$ID = ID
    data.return$parameter = parameter.selected
    data.return$purity = purity
    data.return$ploidy = data.out$ploidy
    
    #' run quality metric on TCGA samples
    FacetsQuality = facetsSuite::run_facets(data.in,
                                            genome = 'hg19',
                                            cval = cval.postprocess, 
                                            dipLogR = dipLogR.purity, 
                                            snp_nbhd = snp.nbhd)
    
    Quality = facetsSuite::check_fit(FacetsQuality)
    fit.out = facets_fit_qc(facets_output = FacetsQuality)
    fit.out_df = data.frame(QC = fit.out$facets_qc,
                            wgd = Quality$wgd,
                            ID = ID)
    
    TCGA_QC_out = rbind(TCGA_QC_out, fit.out_df)
    
    # catch flags, if there are some
    
    facets.flags = data.out$emflags
    if(!is.null(facets.flags)){
      flags.out = c(flags.out, paste0(facets.flags, ';', ID))
    }
    
    TCGA_copyNumber_out = rbind(TCGA_copyNumber_out, data.return)
    
    # fetch cnlr values for chrY; orthogonal validation
    data.cnlr = data.process_out$jointseg
    data.cnlr = data.cnlr[which(data.cnlr$chrom == 24), ]
    data.cnlr$ID = ID
    TCGA_cnLR_out = rbind(TCGA_cnLR_out, data.cnlr)
    
    #' run FacetsY to IGV converstion
    #' create list
    FacetsY_output = list(snps = data.process_out$jointseg,
                          segs = data.out$cncf,
                          dipLogR = data.out$dipLogR)
    
    IGV_out = format_igv_seg(facets_output = FacetsY_output,
                             sample_id = ID)
    
    TCGA_IGV_out = rbind(TCGA_IGV_out, IGV_out)
    
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
    
    TCGA_arm_alteration_out = rbind(TCGA_arm_alteration_out, arm_change.out)
    
    
  })
}

write.table(TCGA_copyNumber_out, file = '/home/kreitzec/Y_chromosome_loss/Loss/TCGA_copyNumber_out.txt', row.names = F, sep = '\t', quote = F)
write.table(flags.out, file = '/home/kreitzec/Y_chromosome_loss/Loss/TCGA_errorflags.txt', row.names = F)
write.table(TCGA_cnLR_out, file = '/home/kreitzec/Y_chromosome_loss/Loss/TCGA_cnLR_out.txt', row.names = F, sep = '\t', quote = F)
write.table(TCGA_IGV_out, file = '/home/kreitzec/Y_chromosome_loss/Loss/TCGA_IGV_out.txt', row.names = F, sep = '\t', quote = F)
write.table(TCGA_arm_alteration_out, file = '/home/kreitzec/Y_chromosome_loss/Loss/TCGA_arm_alteration_out.txt', row.names = F, sep = '\t', quote = F)
write.table(TCGA_QC_out, file = '/home/kreitzec/Y_chromosome_loss/Loss/TCGA_QC_out.txt', row.names = F, sep = '\t', quote = F)

#' qualitatively call y chromosome loss
## qualitative calling of TCGA loss (50% rule) and adherent analysis
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
TCGA.binaryLoss_out = lapply(unique(TCGA_copyNumber_out$ID),
                            function(x) Y_loss_call(data = TCGA_copyNumber_out, sample.id = x))

TCGA.binaryLoss_out = data.frame(do.call('rbind', TCGA.binaryLoss_out))
write.table(TCGA.binaryLoss_out, 
            file = '/juno/home/kreitzec/Y_chromosome_loss/Loss/TCGA.binaryLoss_out.txt', 
            row.names = F,
            sep = '\t',
            quote = F)

## out