## FacetsY on IMPACT  samples;
## This script is intended to be run on the juno
##
## start: 04/16/2021
## revision: 08/25/2021
## revision: 07/08/2022
## revision: 07/11/2022
## revision: 07/12/2022
## revision: 07/13/2022
## revision: 12/20/2022
## chris kreitzer


## install local FacetsY and pctGCdata
require('pctGCdata', '~/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '~/R/x86_64-pc-linux-gnu-library/4.0/')
source('~/Master/Scripts/FacetsQC.R')
source('~/Master/Facets_Helper.R')
library(facetsSuite)
library(dplyr)
library(data.table)


# IMPACT masterfile
IMPACT.samples = read.csv('~/Master/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')
IMPACT.samples = as.character(IMPACT.samples$counts_file)
copyNumber = facetsSuite:::copy_number_states
cn_state_loss = copyNumber$call[which(copyNumber$numeric_call %in% c(-2, -1))]

cval_purity = 100
cval_hisens = 50

##----------------+
## Functions autom.
##----------------+

IMPACT_Y_loss = function(files){
  try({
    CopyNumberCalls = data.frame()
    Cnlr_Y = data.frame()
    Cnlr_X = data.frame()
    IMPACT_IGV_out = data.frame()
    Arm_changes = data.frame()
    
    
    ID = substr(files, start = 59, stop = 75)
    print(ID)
    data.in = readsnpmatrix(files)
    data.in = as.data.frame(data.in)
    
    if(nrow(data.in) == 0) next
    
    #' exclude repetitive stretches, pseudogenes, and centromere
    Y_chromosome = data.in[which(data.in$Chromosome == 'Y' & 
                                   data.in$Position >= 2654800 & 
                                   data.in$Position <= 28000000), ]
    #' exclude PCDH11Y
    Y_chromosome = Y_chromosome[with(Y_chromosome, !((Position %in% 4922131:5612269))), ]
    #' exclude centromere (+CHEK2P1)
    Y_chromosome = Y_chromosome[with(Y_chromosome, !((Position %in% 10500000:14000000))), ]
    #' exclude FAM197Y1
    Y_chromosome = Y_chromosome[with(Y_chromosome, !((Position %in% 9382525:9384756))), ]
    
    ##-------------
    ## gather read-depth data;
    ##-------------
    countmatrix = rbind(data.in[which(data.in$Chromosome %in% c(seq(1, 22, 1), 'X')), ], Y_chromosome)
    set.seed(100)
    segmentation = facetsY::preProcSample(rcmat = countmatrix,
                                          ndepth = 35,
                                          cval = 25, 
                                          het.thresh = 0.25, 
                                          snp.nbhd = 250, 
                                          hetscale = T, 
                                          unmatched = FALSE, 
                                          ndepthmax = 1000,
                                          gbuild = 'hg19')
    
    ##-------------
    ## PURITY-run:
    ##-------------
    data.process = facetsY::procSample(x = segmentation,
                                       cval = cval_purity, 
                                       min.nhet = 15)
    
    data.out = facetsY::emcncf(data.process)
    
    ##-------------
    ## Hisens-run:
    ##-------------
    dipLogR.purity = data.out$dipLogR
    data.process_out = facetsY::procSample(x = segmentation,
                                           cval = cval_hisens,
                                           dipLogR = dipLogR.purity)
    
    data.out = facetsY::emcncf(data.process_out)
    
    
    ##-------------
    ## Gather information
    ##-------------
    data.out$cncf = cbind(data.out$cncf, 
                          cf = data.process_out$out$cf, 
                          tcn = data.process_out$out$tcn, 
                          lcn = data.process_out$out$lcn)
    data.out$cncf$lcn[data.out$cncf$tcn == 1] = 0
    data.out$cncf$lcn.em[data.out$cncf$tcn.em == 1] = 0
    facets_out = list(snps = data.process_out$jointseg,
                      segs = data.out$cncf,
                      purity = as.numeric(data.out$purity),
                      ploidy = as.numeric(data.out$ploidy),
                      dipLogR = data.process_out$dipLogR,
                      alBalLogR = data.process_out$alBalLogR,
                      flags = data.process_out$flags,
                      em_flags = data.out$emflags,
                      loglik = data.out$loglik)
    
    QC_metrics = facets_fit_qc(facets_output = facets_out)
    QC = QC_metrics$facets_qc
    WGD = QC_metrics$wgd
    purity = QC_metrics$purity
    ploidy = QC_metrics$ploidy
    fga = QC_metrics$fga
    
    data_return = facets_out$segs
    data_return$id = ID
    data_return$QC = QC
    data_return$wgd = WGD
    data_return$purity = purity
    data_return$ploidy = ploidy
    data_return$fga
    
    CopyNumberCalls = rbind(CopyNumberCalls, data_return)
    CopyNumberCalls
    
    
    ##------------+
    ## CnLR values for
    ## chromosome X and Y
    ##------------+
    data_cnlr = data.process_out$jointseg
    data.cnlr_Y = data_cnlr[which(data_cnlr$chrom == 24), ]
    data.cnlr_X = data_cnlr[which(data_cnlr$chrom == 23), ]
    data.cnlr_Y$ID = ID
    data.cnlr_X$ID = ID
    Cnlr_Y = rbind(Cnlr_Y, data.cnlr_Y)
    Cnlr_Y
    Cnlr_X = rbind(Cnlr_X, data.cnlr_X)
    Cnlr_X
    
    
    ##------------+
    ## convert FacetsY to
    ## IGV-like format
    ##------------+
    IGV_out = format_igv_seg(facets_output = facets_out, 
                             sample_id = ID, 
                             normalize = T)
    
    IMPACT_IGV_out = rbind(IMPACT_IGV_out, IGV_out)
    IMPACT_IGV_out
    
    
    ##------------+ 
    ## Arm-Level changes
    ##------------+
    arm_change = facetsSuite::arm_level_changes(segs = facets_out$segs,
                                                ploidy = facets_out$ploidy,
                                                genome = 'hg19',
                                                algorithm = 'em')
    arm_change_df = arm_change$full_output
    aneuploidy = filter(arm_change_df, cn_state != 'DIPLOID')
    loss = filter(arm_change_df, cn_state %in% cn_state_loss)
    out = data.frame(id = ID,
                     purity = purity,
                     ploidy = ploidy,
                     genome_doubled = arm_change$genome_doubled,
                     fraction_cna = arm_change$fraction_cna,
                     weighted_fraction_cna = arm_change$weighted_fraction_cna,
                     AS_score = nrow(aneuploidy),
                     losses_n = nrow(loss))
    
    Arm_changes = rbind(Arm_changes, out)
    Arm_changes
    
    ##------------+
    ## OUTPUT
    ##------------+
    output = list(CopyNumberCalls = CopyNumberCalls,
                  Cnlr_Y = Cnlr_Y,
                  Cnlr_X = Cnlr_X,
                  IMPACT_IGV_out = IMPACT_IGV_out,
                  Arm_changes = Arm_changes)
    output
  })
}

xx = lapply(unique(IMPACT_samples), function(x) IMPACT_Y_loss(files = x))
xx = Filter(function(x) length(x) > 1, xx)
CopyNumberStates = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$CopyNumberCalls)))
Cnlr_Y = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$Cnlr_Y)))
Cnlr_X = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$Cnlr_X)))
IMPACT_IGV_out = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$IMPACT_IGV_out)))
Arm_changes = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$Arm_changes)))

##----------------+
## Write Output
##----------------+
write.table(x = CopyNumberStates, file = '~/Master/CopyNumberStates.txt', row.names = F, sep = '\t', quote = F)
write.table(x = Cnlr_Y, file = '~/Master/Cnlr_Y.txt', row.names = F, sep = '\t', quote = F)
write.table(x = Cnlr_X, file = '~/Master/Cnlr_X.txt', row.names = F, sep = '\t', quote = T)
write.table(x = IMPACT_IGV_out, file = '~/Master/IGV_out.seg', row.names = F, sep = '\t', quote = F)
write.table(x = Arm_changes, file = '~/Master/IMPACT_arm_change_out.txt', row.names = F, sep = '\t', quote = F)


#' #' qualitatively call y chromosome loss
#' ## qualitative calling of WES loss (50% rule) and adherent analysis
#' ## Functions:
#' #' calling Y_chromosome loss from Facets output
#' message('evaluating binary Y loss')
#' 
#' Y_loss_call = function(data, sample.id){
#'   data.sub = data[which(data$ID == sample.id), ]
#'   data.sub = data.sub[which(data.sub$chrom == 24), ]
#'   if(nrow(data.sub) == 0){
#'     return(data.frame(cbind(Y_call = NA, 
#'                             sample.id,
#'                             length = NA,
#'                             ploidy = NA,
#'                             purity = NA,
#'                             min.segment = NA,
#'                             max.segment = NA)))
#'   } else {
#'     data.sub$length = data.sub$end - data.sub$start
#'     Y_ratio = sum(data.sub$length[which(data.sub$tcn.em == 0)], na.rm = T) / sum(data.sub$length, na.rm = T)
#'     
#'     # return values:
#'     purity = unique(data.sub$purity)
#'     length = unique(sum(data.sub$length, na.rm = T))
#'     ploidy = unique(data.sub$ploidy)
#'     min.segment = unique(min(data.sub$start))
#'     max.segment = unique(max(data.sub$end))
#'     
#'     Y_call = ifelse(Y_ratio >= 0.5, 'Y_chrom_loss', 'intact_Y_chrom')
#'     Y_call = unique(Y_call)
#'     
#'     sample.id = unique(sample.id)
#'     
#'     return(data.frame(cbind(Y_call, 
#'                             sample.id,
#'                             length,
#'                             ploidy,
#'                             purity,
#'                             min.segment = min.segment,
#'                             max.segment = max.segment)))
#'   }
#' }
#' 
#' ## apply to data, analyzed above
#' IMPACT.binaryLoss_out = lapply(unique(CopyNumber_out$ID),
#'                                function(x) Y_loss_call(data = CopyNumber_out, sample.id = x))
#' 
#' IMPACT.binaryLoss_out = data.frame(do.call('rbind', IMPACT.binaryLoss_out))
#' write.table(IMPACT.binaryLoss_out, 
#'             file = '/juno/home/kreitzec/Master/Data/Loss/IMPACT.binaryLoss_out.txt', 
#'             row.names = F,
#'             sep = '\t')
#' 
#' 
#' ## out
