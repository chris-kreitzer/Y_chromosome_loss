##----------------+
## FacetsY on TCGA samples
## Juno-scirpt
##----------------+

## run purity run (cval = 500) first; then get dipLogR and run with cval 200
## exclude base coverage > 30 Mb (see coverage plot) (PAR2 region not relevant)
## run FacetsQC and flag samples with QC == FALSE

## start: 03/07/2023
##
## chris Kreitzer


# install local FacetsY and pctGCdata
## install local FacetsY and pctGCdata
require('pctGCdata', '~/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '~/R/x86_64-pc-linux-gnu-library/4.0/')
source('~/Master/Scripts/FacetsQC.R')
source('~/Master/Facets_Helper.R')
library(facetsSuite)
library(dplyr)
library(data.table)



# samplePaths = read.csv('Data/07_TCGA/TCGA_paths.txt', sep = '\t')
# prePrint = read.csv('Data/07_TCGA/LOY_PanCancerTCGA.txt', sep = '\t')
# prePrint_samples = prePrint$case_tcga_sample_id[which(prePrint$Flag == "")]
# 
# 
# TCGA_paths = data.frame()
# for(i in unique(prePrint_samples)){
#   print(i)
#   id = i
#   path = samplePaths$path[grep(pattern = i, x = samplePaths$path)]
#   cohort = prePrint$Cohort[which(prePrint$case_tcga_sample_id == i)]
#   age = prePrint$Age[which(prePrint$case_tcga_sample_id == i)]
#   control = prePrint$control_tcga_sample_id[which(prePrint$case_tcga_sample_id == i)]
#   path = ifelse(length(path) > 1, path[grep(pattern = control, x = path)], path)
#   
#   out = data_frame(id = id,
#                    path = path,
#                    control = control,
#                    cohort = cohort,
#                    age = age)
#   TCGA_paths = rbind(TCGA_paths, out)
# }
# 
# TCGA_paths = TCGA_paths[!is.na(TCGA_paths$path), ]
# write.table(x = TCGA_paths, file = 'Data/07_TCGA/TCGA_paths.txt', sep = '\t', row.names = F, quote = F)



# load FACETS countsfiles for Prostate WES samples
TCGA.samples = read.csv('~/TCGA_paths.txt', sep = '\t')
TCGA.samples = as.character(TCGA.samples$path)
head(TCGA.samples)

TCGA.samples = TCGA.samples[1:2]

copyNumber = facetsSuite:::copy_number_states
cn_state_loss = copyNumber$call[which(copyNumber$numeric_call %in% c(-2, -1))]

cval_purity = 500
cval_hisens = 200


##----------------+
## Functions autom.
##----------------+

WES_Y_loss = function(files){
  try({
    CopyNumberCalls = data.frame()
    Cnlr_Y = data.frame()
    Cnlr_X = data.frame()
    IMPACT_IGV_out = data.frame()
    Arm_changes = data.frame()
    
    
    ID = substr(x = basename(files), start = 1, stop = 36)
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

xx = lapply(unique(TCGA.samples), function(x) WES_Y_loss(files = x))
xx = Filter(function(x) length(x) > 1, xx)
CopyNumberStates = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$CopyNumberCalls)))
Cnlr_Y = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$Cnlr_Y)))
Cnlr_X = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$Cnlr_X)))
WES_IGV_out = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$IMPACT_IGV_out)))
Arm_changes = setDF(rbindlist(lapply(1:length(xx), function(x) xx[[x]]$Arm_changes)))

##----------------+
## Write Output
##----------------+
write.table(x = CopyNumberStates, file = '~/Master/Data/TCGA/WES_CopyNumberStates.txt', row.names = F, sep = '\t', quote = F)
write.table(x = Cnlr_Y, file = '~/Master/Data/TCGA/WES_Cnlr_Y.txt', row.names = F, sep = '\t', quote = F)
write.table(x = Cnlr_X, file = '~/Master/Data/TCGA/WES_Cnlr_X.txt', row.names = F, sep = '\t', quote = T)
write.table(x = WES_IGV_out, file = '~/Master/Data/TCGA/WES_IGV_out.seg', row.names = F, sep = '\t', quote = F)
write.table(x = Arm_changes, file = '~/Master/Data/TCGA/WES_arm_change_out.txt', row.names = F, sep = '\t', quote = F)


#' out