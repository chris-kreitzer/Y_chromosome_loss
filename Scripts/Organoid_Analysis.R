## Investigate 5 IMPACT samples; of which organoids where made (Harisha);
## Here, we will investigate the IMPACT samples to call Y-chromosome loss
## 
## later we will also check whether we can apply FacetsY on those sequenced
## organoids to call the same status as we confirmed in the lab
## 
## Harisha confirmed a Y-chromosome loss in Organoids with FISH
## 
## Start: 04/23/2021


rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/Harisha_Subanalysis/')
source('../Y_chromosome_loss/FacetsQC.R')
source('../Y_chromosome_loss/Plotting_theme.R')
set.seed(395)


## Libraries
library(facetsY)
library(facetsSuite)
library(pctGCdata)


## Input:
Sample.path_df = read.csv('../05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
Samples_to_investigate = c('P-0003033-T02-IM5',
                           'P-0001445-T01-IM3',
                           'P-0005953-T01-IM5',
                           'P-0001899-T01-IM3',
                           'P-0004597-T01-IM5')

Path = c()
for(i in unique(Samples_to_investigate)){
  gre = grep(pattern = i, x = Sample.path_df$counts_file, value = T)
  Path = c(Path, gre)
}


## Example run:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## P-0003033-T02-IM5 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Metastasis Liver:
# Fraction Genome altered (cBIO): 0.3121
# Coverage: 596
# Purity: 50

# 3 Mutations (TP53, NOTCH1, MDM4)
# 
# 1p lost; 1q gain
# 3q gain
# 8p lost; 8q gain
# 13q lost
# 14q gain

P1_FSuite = facetsSuite::read_snp_matrix(input_file = 'Harisha_Subanalysis/countsMerged____P-0003033-T02-IM5_P-0003033-N01-IM5.dat.gz')
P1_FSuite = facetsSuite::run_facets(read_counts = P1_FSuite,
                                    cval = 100, 
                                    genome = 'hg19')
P1_seg = facetsSuite::format_igv_seg(facets_output = P1_FSuite, sample_id = 'P-0003033')
write.table(P1_seg, file = 'Harisha_Subanalysis/P-0003033.seg', sep = '\t', row.names = F, quote = F)

# the copy number aberration profile between cBIO and Facets-segmented profile looks pretty much the same;
# All of the above mentioned chromosome-arm changes are found within Facets-Segmentation
P1_arm.change = facetsSuite::arm_level_changes(segs = P1_FSuite$segs,
                                               ploidy = P1_FSuite$ploidy, 
                                               genome = 'hg19', 
                                               algorithm = 'cncf')
P1_arm.change$full_output
P1_arm.change$genome_doubled
P1_arm.change$fraction_cna # slightly overestimated - compared to cBIO (0.3121)

# the arm-level alteration calling method is not suitable in looking into chromosome-arm
# changes, because it captures just segments which are >80% of respective arms.
# This means, if there is 75% of chromosome 1q altered (clearly) it will not appear within 
# the facetsSuite::arm_level_change() output:


P1_gene.change = facetsSuite::gene_level_changes(facets_output = P1_FSuite,
                                                 genome = 'hg19')
P1_gene.change

## if we work on selected genes; this output might provide useful insights
## looks at every segment (called by Facets) and map respective genes onto the segments
## the median cnlr (of given segment) as well as copy numbers are provided within the output;
## also suggested cn_state will be found
## the number of heterozygous SNP on respective gene as well as the number of markers on respective genes is found;

# the facetsSuite::closeup_plot() turns out to be super useful, when someone wants to inspect certain genes
# e.g. the copy number of PTEN on chromosome 10 (it highlights the individual cnlr at this particular gene)
# use this function to run on chrY (if possible)
x = facetsSuite::closeup_plot(facets_data = P1_FSuite, highlight_gene = 'TP53')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## try running the same with FacetsY:
P1_FY = facetsY::readSnpMatrix(filename = 'Harisha_Subanalysis/countsMerged____P-0003033-T02-IM5_P-0003033-N01-IM5.dat.gz')
P1_Y = facetsY::preProcSample(P1_FY, gbuild = 'hg19')
P1_Y2 = facetsY::procSample(P1_Y, cval = 100)
P1_Y3 = facetsY::emcncf(P1_Y2)

P1_Y3$cncf$lcn.em = ifelse(is.na(P1_Y3$cncf$lcn.em), 0, P1_Y3$cncf$lcn.em)


P1_manual = list(P1_Y3$cncf,
                 P1_Y3$purity,
                 P1_Y3$ploidy,
                 P1_Y2$jointseg,
                 P1_Y3$dipLogR)
names(P1_manual) = c('segs', 'purity', 'ploidy', 'snps', 'dipLogR')

## Arm_manual (facetsY)
P1_arm_manual = facetsSuite::arm_level_changes(segs = P1_manual$segs, ploidy = P1_manual$ploidy, genome = 'hg19')
P1_IGV_manual = facetsSuite::format_igv_seg(facets_output = P1_manual, sample_id = 'FacetsY')
P1_IGV_manual
P1_manual$segs
write.table(P1_IGV_manual, file = 'Harisha_Subanalysis/P-0003033.manual.seg', sep = '\t', row.names = F, quote = F)




##-----------------------------------------------------------------------------
P1 = facetsY::readSnpMatrix(filename = 'Harisha_Subanalysis/countsMerged____P-0005953-T01-IM5_P-0005953-N01-IM5.dat.gz')
P1_pre = facetsY::preProcSample(rcmat = P1, 
                                ndepth = 35, 
                                het.thresh = 0.25, 
                                snp.nbhd = 250, 
                                cval = 25, 
                                gbuild = 'hg19', 
                                hetscale = TRUE)
P1_post = facetsY::procSample(x = P1_pre, 
                              cval = 100,
                              min.nhet = 15)
P1_fit = facetsY::emcncf(P1_post)

## sub modifications:
P1_fit$cncf = cbind(P1_fit$cncf, cf = P1_post$out$cf, tcn = P1_post$out$tcn, lcn = P1_post$out$lcn)
P1_fit$cncf$lcn[P1_fit$cncf$tcn == 1] = 0
P1_fit$cncf$lcn.em[P1_fit$cncf$tcn.em == 1] = 0

P1_out = list(segs = P1_fit$cncf,
              snps = P1_post$jointseg,
              purity = P1_fit$purity,
              ploidy = P1_fit$ploidy,
              dipLogR = P1_fit$dipLogR)

## output:
P1_segment = P1_fit$cncf[which(P1_fit$cncf$chrom == 24), ]
P1_purity = P1_fit$purity
P1_ploidy = P1_fit$ploidy
P1_fit$dipLogR
P1_IGV_seg = facetsSuite::format_igv_seg(facets_output = P1_out, sample_id = 'P-0005953-T01-IM3')

# P-0005953-T01 = intact!

P1_Markers = P1_post$jointseg[which(P1_post$jointseg$chrom == 24), ]
P1_Markers = P1_Markers[,c('maploc', 'rCountT', 'rCountN', 'het')]
P1_Markers$ID = 'P-0005953'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P2 = facetsY::readSnpMatrix(filename = 'Harisha_Subanalysis/countsMerged____P-0001445-T01-IM3_P-0001445-N01-IM3.dat.gz')
P2_pre = facetsY::preProcSample(rcmat = P2, 
                                ndepth = 35, 
                                het.thresh = 0.25, 
                                snp.nbhd = 250, 
                                cval = 25, 
                                gbuild = 'hg19', 
                                hetscale = TRUE)
P2_post = facetsY::procSample(x = P2_pre, 
                              cval = 100,
                              min.nhet = 15)
P2_fit = facetsY::emcncf(P2_post)

## sub modifications:
P2_fit$cncf = cbind(P2_fit$cncf, cf = P2_post$out$cf, tcn = P2_post$out$tcn, lcn = P2_post$out$lcn)
P2_fit$cncf$lcn[P2_fit$cncf$tcn == 1] = 0
P2_fit$cncf$lcn.em[P2_fit$cncf$tcn.em == 1] = 0

P2_out = list(segs = P2_fit$cncf,
              snps = P2_post$jointseg,
              purity = P2_fit$purity,
              ploidy = P2_fit$ploidy,
              dipLogR = P2_fit$dipLogR)

## output:
P2_segment = P2_fit$cncf[which(P2_fit$cncf$chrom == 24), ]
P2_purity = P2_fit$purity
P2_ploidy = P2_fit$ploidy
P2_IGV_seg = facetsSuite::format_igv_seg(facets_output = P2_out, sample_id = 'P-0001445-T01-IM3')
P2_IGV_unadjusted_seg = facetsSuite::format_igv_seg(facets_output = P2_out, sample_id = 'P-0001445-unadjusted', normalize = F)

# P-0001445-T01 = intact! (most likely doubled)

P2_Markers = P2_post$jointseg[which(P2_post$jointseg$chrom == 24), ]
P2_Markers = P2_Markers[,c('maploc', 'rCountT', 'rCountN', 'het')]
P2_Markers$ID = 'P-0001445'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P3 = facetsY::readSnpMatrix(filename = 'Harisha_Subanalysis/countsMerged____P-0001899-T01-IM3_P-0001899-N01-IM3.dat.gz')
P3_pre = facetsY::preProcSample(rcmat = P3, 
                                ndepth = 35, 
                                het.thresh = 0.25, 
                                snp.nbhd = 250, 
                                cval = 25, 
                                gbuild = 'hg19', 
                                hetscale = TRUE)
P3_post = facetsY::procSample(x = P3_pre, 
                              cval = 100,
                              min.nhet = 15)
P3_fit = facetsY::emcncf(P3_post)

## sub modifications:
P3_fit$cncf = cbind(P3_fit$cncf, cf = P3_post$out$cf, tcn = P3_post$out$tcn, lcn = P3_post$out$lcn)
P3_fit$cncf$lcn[P3_fit$cncf$tcn == 1] = 0
P3_fit$cncf$lcn.em[P3_fit$cncf$tcn.em == 1] = 0

P3_out = list(segs = P3_fit$cncf,
              snps = P3_post$jointseg,
              purity = P3_fit$purity,
              ploidy = P3_fit$ploidy,
              dipLogR = P3_fit$dipLogR)

## output:
P3_segment = P3_fit$cncf[which(P3_fit$cncf$chrom == 24), ]
P3_purity = P3_fit$purity
P3_ploidy = P3_fit$ploidy
P3_IGV_seg = facetsSuite::format_igv_seg(facets_output = P3_out, sample_id = 'P-0001899-T01-IM3')
P3_IGV_seg.unadjusted = facetsSuite::format_igv_seg(facets_output = P3_out, sample_id = 'P-0001899_unadjusted', normalize = F)


# P-0001899-T01 = intact! (most likely doubled)

P3_Markers = P3_post$jointseg[which(P3_post$jointseg$chrom == 24), ]
P3_Markers = P3_Markers[,c('maploc', 'rCountT', 'rCountN', 'het')]
P3_Markers$ID = 'P-0001899'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P4 = facetsY::readSnpMatrix(filename = 'Harisha_Subanalysis/countsMerged____P-0003033-T02-IM5_P-0003033-N01-IM5.dat.gz')
P4_pre = facetsY::preProcSample(rcmat = P4, 
                                ndepth = 35, 
                                het.thresh = 0.25, 
                                snp.nbhd = 250, 
                                cval = 25, 
                                gbuild = 'hg19', 
                                hetscale = TRUE)
P4_post = facetsY::procSample(x = P4_pre, 
                              cval = 100,
                              min.nhet = 15)
P4_fit = facetsY::emcncf(P4_post)

## sub modifications:
P4_fit$cncf = cbind(P4_fit$cncf, cf = P4_post$out$cf, tcn = P4_post$out$tcn, lcn = P4_post$out$lcn)
P4_fit$cncf$lcn[P4_fit$cncf$tcn == 1] = 0
P4_fit$cncf$lcn.em[P4_fit$cncf$tcn.em == 1] = 0

P4_out = list(segs = P4_fit$cncf,
              snps = P4_post$jointseg,
              purity = P4_fit$purity,
              ploidy = P4_fit$ploidy,
              dipLogR = P4_fit$dipLogR)

QC = facetsSuite::check_fit(facets_output = P)
## output:
P4_segment = P4_fit$cncf[which(P4_fit$cncf$chrom == 24), ]
P4_purity = P4_fit$purity
P4_ploidy = P4_fit$ploidy
P4_IGV_seg = facetsSuite::format_igv_seg(facets_output = P4_out, sample_id = 'P-0003033-T01-IM5')

# P-0003033-T01 = intact! (most likely doubled)

P4_Markers = P4_post$jointseg[which(P4_post$jointseg$chrom == 24), ]
P4_Markers = P4_Markers[,c('maploc', 'rCountT', 'rCountN', 'het')]
P4_Markers$ID = 'P-0003033'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P5 = facetsY::readSnpMatrix(filename = 'Harisha_Subanalysis/countsMerged____P-0004597-T01-IM5_P-0004597-N01-IM5.dat.gz')
P5_pre = facetsY::preProcSample(rcmat = P5, 
                                ndepth = 35, 
                                het.thresh = 0.25, 
                                snp.nbhd = 250, 
                                cval = 25, 
                                gbuild = 'hg19', 
                                hetscale = TRUE)
P5_post = facetsY::procSample(x = P5_pre, 
                              cval = 100,
                              min.nhet = 15)
P5_fit = facetsY::emcncf(P5_post)

## sub modifications:
P5_fit$cncf = cbind(P5_fit$cncf, cf = P5_post$out$cf, tcn = P5_post$out$tcn, lcn = P5_post$out$lcn)
P5_fit$cncf$lcn[P5_fit$cncf$tcn == 1] = 0
P5_fit$cncf$lcn.em[P5_fit$cncf$tcn.em == 1] = 0

P5_out = list(segs = P5_fit$cncf,
              snps = P5_post$jointseg,
              purity = P5_fit$purity,
              ploidy = P5_fit$ploidy,
              dipLogR = P5_fit$dipLogR)

## output:
P5_segment = P5_fit$cncf[which(P5_fit$cncf$chrom == 24), ]
P5_purity = P5_fit$purity
P5_ploidy = P5_fit$ploidy
P5_IGV_seg = facetsSuite::format_igv_seg(facets_output = P5_out, sample_id = 'P-0004597-T01-IM5')

# P-0004597-T01 = intact! (most likely doubled)

P5_Markers = P5_post$jointseg[which(P5_post$jointseg$chrom == 24), ]
P5_Markers = P5_Markers[,c('maploc', 'rCountT', 'rCountN', 'het')]
P5_Markers$ID = 'P-0004597'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create common IGV file
Harisha_IGV = rbind(P1_IGV_seg,
                    P2_IGV_seg,
                    P3_IGV_seg,
                    P4_IGV_seg,
                    P5_IGV_seg)

Harisha_IGV_unadjusted = rbind(P2_IGV_unadjusted_seg, P3_IGV_seg.unadjusted)

write.table(Harisha_IGV, file = 'Harisha_Subanalysis/5_Samples.FacetsY.seg', row.names = F, quote = F, sep = '\t')
write.table(Harisha_IGV_unadjusted, file = 'Harisha_Subanalysis/2samples_unadjusted.seg', sep = '\t', quote = F, row.names = F)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# common marker within the 5 samples;
common_marker = rbind(P1_Markers, 
                      P2_Markers,
                      P3_Markers,
                      P4_Markers,
                      P5_Markers)
x = common_marker %>%
  group_by(maploc) %>%
  summarise(n = n(),
            n_rel = n / 5,
            ave.coverage = mean(rCountT, na.rm = T))
x$threshold = ifelse(x$n_rel >= 0.75, 1, 0)

View(x)








##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sequencing of Organoid samples;
library(pctGCdata)
source('../Y_chromosome_loss/FacetsQC.R')
source('../Y_chromosome_loss/Plotting_theme.R')

setwd('~/Documents/MSKCC/04_Y_chromo_loss/Harisha_Subanalysis/')
list.files()

ASC1 = facetsY::readSnpMatrix(filename = 'ASC1.gz')
ASC1Y = ASC1[which(ASC1$Chromosome == 'Y' & ASC1$NOR.DP >= 35), ]
a.x = log2(ASC1Y$TUM.DP / ASC1Y$NOR.DP)
mean(a.x)

ASC1_facetsSuite = facetsSuite::run_facets(read_counts = ASC1, 
                                           cval = 100, 
                                           genome = 'hg19')

facetsSuite::cnlr_plot(ASC1_facetsSuite)
ASC1_QC = facets_fit_qc(facets_output = ASC1_facetsSuite)

ASC1_pre = facetsY::preProcSample(rcmat = ASC1,
                                  ndepth = 35,
                                  het.thresh = 0.25,
                                  snp.nbhd = 250,
                                  cval = 25,
                                  gbuild = 'hg19',
                                  hetscale = TRUE)

ASC1_post = facetsY::procSample(x = ASC1_pre,
                                cval = 100,
                                min.nhet = 15)

ASC1_fit = facetsY::emcncf(ASC1_post)
View(ASC1_fit$cncf[which(ASC1_fit$cncf$chrom == 24), ])


## sub modifications:
ASC1_fit$cncf = cbind(ASC1_fit$cncf, cf = ASC1_post$out$cf, tcn = ASC1_post$out$tcn, lcn = ASC1_post$out$lcn)
ASC1_fit$cncf$lcn[ASC1_fit$cncf$tcn == 1] = 0
ASC1_fit$cncf$lcn.em[ASC1_fit$cncf$tcn.em == 1] = 0

ASC1_out = list(segs = ASC1_fit$cncf,
              snps = ASC1_post$jointseg,
              purity = ASC1_fit$purity,
              ploidy = ASC1_fit$ploidy,
              dipLogR = ASC1_fit$dipLogR)

ASC1_IGV_seg = facetsSuite::format_igv_seg(facets_output = ASC1_out, sample_id = 'ASC1')

write.table(ASC1_IGV_seg, file = 'ASC1.seg', row.names = F, quote = F, sep = '\t')

## GMM on Chromosome Y
ASC1_local = ASC1_post$jointseg[which(ASC1_post$jointseg$chrom == 24), ]


GMM.processed = mixtools::normalmixEM(x = ASC1_local$cnlr, k = 2)
GMM.lambda = GMM.processed$lambda
GMM.mu = GMM.processed$mu

# make plot
plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}

cnlr = ASC1_local$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio') +
  #     title = 'WES Cn-LogR Distribution across the Y-chromosome',
   #    subtitle = '51.1% of markers belong to \'red\' distribution with mu = -3.25') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))


##-------------------------------------------------------------------
BM61 = facetsY::readSnpMatrix(filename = 'BM61.gz')
BM61Y = BM61[which(BM61$Chromosome == 'Y' & BM61$NOR.DP > 35), ]
b.x = log2(BM61Y$TUM.DP / BM61Y$NOR.DP)
hist(b.x, nclass = 100, xlab = 'cnlr', ylab = 'density', freq = F, main = 'T/N cnlr signal distribution: BM61')
box()

BM61_facetsSuite = facetsSuite::run_facets(read_counts = BM61, 
                                           cval = 100, 
                                           genome = 'hg19')

facetsSuite::cnlr_plot(BM61_facetsSuite)
BM61_QC = facets_fit_qc(facets_output = BM61_facetsSuite)

BM61_pre = facetsY::preProcSample(rcmat = BM61,
                                  ndepth = 35,
                                  het.thresh = 0.25,
                                  snp.nbhd = 250,
                                  cval = 25,
                                  gbuild = 'hg19',
                                  hetscale = TRUE)

BM61_post = facetsY::procSample(x = BM61_pre,
                                cval = 100,
                                min.nhet = 15)

BM61_fit = facetsY::emcncf(BM61_post)
View(BM61_fit$cncf[which(BM61_fit$cncf$chrom == 24), ])

## sub modifications:
BM61_fit$cncf = cbind(BM61_fit$cncf, cf = BM61_post$out$cf, tcn = BM61_post$out$tcn, lcn = BM61_post$out$lcn)
BM61_fit$cncf$lcn[BM61_fit$cncf$tcn == 1] = 0
BM61_fit$cncf$lcn.em[BM61_fit$cncf$tcn.em == 1] = 0

BM61_out = list(segs = BM61_fit$cncf,
                snps = BM61_post$jointseg,
                purity = BM61_fit$purity,
                ploidy = BM61_fit$ploidy,
                dipLogR = BM61_fit$dipLogR)

BM61_IGV_seg = facetsSuite::format_igv_seg(facets_output = BM61_out, sample_id = 'BM61')

write.table(BM61_IGV_seg, file = 'BM61.seg', row.names = F, quote = F, sep = '\t')

## GMM on Chromosome Y
BM61_local = BM61_post$jointseg[which(BM61_post$jointseg$chrom == 24), ]


GMM.processed = mixtools::normalmixEM(x = BM61_local$cnlr, k = 2)
GMM.lambda = GMM.processed$lambda
GMM.mu = GMM.processed$mu

# make plot
plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}

cnlr = BM61_local$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 100,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio',
       title = 'BM61') +
  #    subtitle = '51.1% of markers belong to \'red\' distribution with mu = -3.25') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))
WES.mix.example



##-------------------------------------------------------------------
ST111 = facetsY::readSnpMatrix(filename = 'ST111.gz')
ST111Y = ST111[which(ST111$Chromosome == 'Y' & ST111$NOR.DP >= 35), ]
st111.x = log2(ST111Y$TUM.DP / ST111Y$NOR.DP)
median(st111.x)

ST111_facetsSuite = facetsSuite::run_facets(read_counts = ST111, 
                                           cval = 100, 
                                           genome = 'hg19')
facetsSuite::cnlr_plot(ST111_facetsSuite)
ST111_QC = facets_fit_qc(facets_output = ST111_facetsSuite)

ST111_pre = facetsY::preProcSample(rcmat = ST111,
                                  ndepth = 35,
                                  het.thresh = 0.25,
                                  snp.nbhd = 250,
                                  cval = 25,
                                  gbuild = 'hg19',
                                  hetscale = TRUE)

ST111_post = facetsY::procSample(x = ST111_pre,
                                cval = 100,
                                min.nhet = 15)

ST111_fit = facetsY::emcncf(ST111_post)
View(ST111_fit$cncf[which(ST111_fit$cncf$chrom == 24), ])

## sub modifications:
ST111_fit$cncf = cbind(ST111_fit$cncf, cf = ST111_post$out$cf, tcn = ST111_post$out$tcn, lcn = ST111_post$out$lcn)
ST111_fit$cncf$lcn[ST111_fit$cncf$tcn == 1] = 0
ST111_fit$cncf$lcn.em[ST111_fit$cncf$tcn.em == 1] = 0

ST111_out = list(segs = ST111_fit$cncf,
                snps = ST111_post$jointseg,
                purity = ST111_fit$purity,
                ploidy = ST111_fit$ploidy,
                dipLogR = ST111_fit$dipLogR)

ST111_IGV_seg = facetsSuite::format_igv_seg(facets_output = ST111_out, sample_id = 'ST111')

write.table(ST111_IGV_seg, file = 'ST111.seg', row.names = F, quote = F, sep = '\t')

## GMM on Chromosome Y
ST111_local = ST111_post$jointseg[which(ST111_post$jointseg$chrom == 24), ]


GMM.processed = mixtools::normalmixEM(x = ST111_local$cnlr, k = 2)
GMM.lambda = GMM.processed$lambda
GMM.mu = GMM.processed$mu

# make plot
plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}

cnlr = ST111_local$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio') +
  #     title = 'WES Cn-LogR Distribution across the Y-chromosome',
  #    subtitle = '51.1% of markers belong to \'red\' distribution with mu = -3.25') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))
WES.mix.example



##-------------------------------------------------------------------
ST121 = facetsY::readSnpMatrix(filename = 'ST121.gz')
ST121Y = ST121[which(ST121$Chromosome == 'Y' & ST121$NOR.DP >= 35), ]
ST121.x = log2(ST121Y$TUM.DP / ST121Y$NOR.DP)
mean(ST121.x)

ST121_facetsSuite = facetsSuite::run_facets(read_counts = ST121, 
                                            cval = 100, 
                                            genome = 'hg19')
facetsSuite::cnlr_plot(ST121_facetsSuite)
ST121_QC = facets_fit_qc(facets_output = ST121_facetsSuite)

ST121_pre = facetsY::preProcSample(rcmat = ST121,
                                   ndepth = 35,
                                   het.thresh = 0.25,
                                   snp.nbhd = 250,
                                   cval = 25,
                                   gbuild = 'hg19',
                                   hetscale = TRUE)

ST121_post = facetsY::procSample(x = ST121_pre,
                                 cval = 100,
                                 min.nhet = 15)

ST121_fit = facetsY::emcncf(ST121_post)
View(ST121_fit$cncf[which(ST121_fit$cncf$chrom == 24), ])

## sub modifications:
ST121_fit$cncf = cbind(ST121_fit$cncf, cf = ST121_post$out$cf, tcn = ST121_post$out$tcn, lcn = ST121_post$out$lcn)
ST121_fit$cncf$lcn[ST121_fit$cncf$tcn == 1] = 0
ST121_fit$cncf$lcn.em[ST121_fit$cncf$tcn.em == 1] = 0

ST121_out = list(segs = ST121_fit$cncf,
                 snps = ST121_post$jointseg,
                 purity = ST121_fit$purity,
                 ploidy = ST121_fit$ploidy,
                 dipLogR = ST121_fit$dipLogR)

ST121_IGV_seg = facetsSuite::format_igv_seg(facets_output = ST121_out, sample_id = 'ST121')

write.table(ST121_IGV_seg, file = 'ST121.seg', row.names = F, quote = F, sep = '\t')

## GMM on Chromosome Y
ST121_local = ST121_post$jointseg[which(ST121_post$jointseg$chrom == 24), ]
GMM.processed = mixtools::normalmixEM(x = ST121_local$cnlr, k = 2)
GMM.lambda = GMM.processed$lambda
GMM.mu = GMM.processed$mu

# make plot
plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}

mixmdl = mixtools::normalmixEM(st111.x, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio') +
  #     title = 'WES Cn-LogR Distribution across the Y-chromosome',
  #    subtitle = '51.1% of markers belong to \'red\' distribution with mu = -3.25') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))
WES.mix.example





##-----------------------------------------------------------------------------
##-------------------------------------------------------------------
ST171 = facetsY::readSnpMatrix(filename = 'ST171.gz')
ST171Y = ST171[which(ST171$Chromosome == 'Y' & ST171$NOR.DP >= 35), ]
ST171.x = log2(ST171Y$TUM.DP / ST171Y$NOR.DP)
median(ST171.x)

# analysis pipeline
ST171_facetsSuite = facetsSuite::run_facets(read_counts = ST171, 
                                            cval = 100, 
                                            genome = 'hg19')
facetsSuite::cnlr_plot(ST171_facetsSuite)
ST171_QC = facets_fit_qc(facets_output = ST171_facetsSuite)

ST171_pre = facetsY::preProcSample(rcmat = ST171,
                                   ndepth = 35,
                                   het.thresh = 0.25,
                                   snp.nbhd = 250,
                                   cval = 25,
                                   gbuild = 'hg19',
                                   hetscale = TRUE)

ST171_post = facetsY::procSample(x = ST171_pre,
                                 cval = 100,
                                 min.nhet = 15)

ST171_fit = facetsY::emcncf(ST171_post)
View(ST171_fit$cncf[which(ST171_fit$cncf$chrom == 24), ])

## sub modifications:
ST171_fit$cncf = cbind(ST171_fit$cncf, cf = ST171_post$out$cf, tcn = ST171_post$out$tcn, lcn = ST171_post$out$lcn)
ST171_fit$cncf$lcn[ST171_fit$cncf$tcn == 1] = 0
ST171_fit$cncf$lcn.em[ST171_fit$cncf$tcn.em == 1] = 0

ST171_out = list(segs = ST171_fit$cncf,
                 snps = ST171_post$jointseg,
                 purity = ST171_fit$purity,
                 ploidy = ST171_fit$ploidy,
                 dipLogR = ST171_fit$dipLogR)

ST171_IGV_seg = facetsSuite::format_igv_seg(facets_output = ST171_out, sample_id = 'ST171')

write.table(ST171_IGV_seg, file = 'ST171.seg', row.names = F, quote = F, sep = '\t')

## GMM on Chromosome Y
ST171_local = ST171_post$jointseg[which(ST171_post$jointseg$chrom == 24), ]
GMM.processed = mixtools::normalmixEM(x = ST171_local$cnlr, k = 2)
GMM.lambda = GMM.processed$lambda
GMM.mu = GMM.processed$mu

# make plot
plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}

cnlr = ST171_local$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio') +
  #     title = 'WES Cn-LogR Distribution across the Y-chromosome',
  #    subtitle = '51.1% of markers belong to \'red\' distribution with mu = -3.25') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))
WES.mix.example



##-------------------------------------------------------------------
ST43 = facetsY::readSnpMatrix(filename = 'ST43.gz')
ST43Y = ST43[which(ST43$Chromosome == 'Y' & ST43$NOR.DP >= 35), ]
ST43.x = log2(ST43Y$TUM.DP / ST43Y$NOR.DP)
median(ST43.x)

# analysis pipeline
ST43_facetsSuite = facetsSuite::run_facets(read_counts = ST43, 
                                            cval = 200, 
                                            genome = 'hg19')
facetsSuite::cnlr_plot(ST43_facetsSuite)
ST43_QC = facets_fit_qc(facets_output = ST43_facetsSuite)

ST43_pre = facetsY::preProcSample(rcmat = ST43,
                                   ndepth = 35,
                                   het.thresh = 0.25,
                                   snp.nbhd = 250,
                                   cval = 25,
                                   gbuild = 'hg19',
                                   hetscale = TRUE)

ST43_post = facetsY::procSample(x = ST43_pre,
                                 cval = 200,
                                 min.nhet = 15)

ST43_fit = facetsY::emcncf(ST43_post)
View(ST43_fit$cncf[which(ST43_fit$cncf$chrom == 24), ])

## sub modifications:
ST43_fit$cncf = cbind(ST43_fit$cncf, cf = ST43_post$out$cf, tcn = ST43_post$out$tcn, lcn = ST43_post$out$lcn)
ST43_fit$cncf$lcn[ST43_fit$cncf$tcn == 1] = 0
ST43_fit$cncf$lcn.em[ST43_fit$cncf$tcn.em == 1] = 0

ST43_out = list(segs = ST43_fit$cncf,
                 snps = ST43_post$jointseg,
                 purity = ST43_fit$purity,
                 ploidy = ST43_fit$ploidy,
                 dipLogR = ST43_fit$dipLogR)

ST43_IGV_seg = facetsSuite::format_igv_seg(facets_output = ST43_out, sample_id = 'ST43')

write.table(ST43_IGV_seg, file = 'ST43.seg', row.names = F, quote = F, sep = '\t')

## GMM on Chromosome Y
ST43_local = ST43_post$jointseg[which(ST43_post$jointseg$chrom == 24), ]
GMM.processed = mixtools::normalmixEM(x = ST43_local$cnlr, k = 2)
GMM.lambda = GMM.processed$lambda
GMM.mu = GMM.processed$mu

# make plot
plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}

cnlr = ST43_local$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio') +
  #     title = 'WES Cn-LogR Distribution across the Y-chromosome',
  #    subtitle = '51.1% of markers belong to \'red\' distribution with mu = -3.25') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))
WES.mix.example


##-------------------------------------------------------------------
ST60 = facetsY::readSnpMatrix(filename = 'ST60.gz')
ST60Y = ST60[which(ST60$Chromosome == 'Y' & ST60$NOR.DP >= 35), ]
ST60.x = log2(ST60Y$TUM.DP / ST60Y$NOR.DP)
median(ST60.x)

# analysis pipeline
ST60_facetsSuite = facetsSuite::run_facets(read_counts = ST60, 
                                           cval = 200, 
                                           genome = 'hg19')
facetsSuite::cnlr_plot(ST60_facetsSuite)
ST60_QC = facets_fit_qc(facets_output = ST60_facetsSuite)

ST60_pre = facetsY::preProcSample(rcmat = ST60,
                                  ndepth = 35,
                                  het.thresh = 0.25,
                                  snp.nbhd = 250,
                                  cval = 25,
                                  gbuild = 'hg19',
                                  hetscale = TRUE)

ST60_post = facetsY::procSample(x = ST60_pre,
                                cval = 200,
                                min.nhet = 15)

ST60_fit = facetsY::emcncf(ST60_post)
View(ST60_fit$cncf[which(ST60_fit$cncf$chrom == 24), ])

## sub modifications:
ST60_fit$cncf = cbind(ST60_fit$cncf, cf = ST60_post$out$cf, tcn = ST60_post$out$tcn, lcn = ST60_post$out$lcn)
ST60_fit$cncf$lcn[ST60_fit$cncf$tcn == 1] = 0
ST60_fit$cncf$lcn.em[ST60_fit$cncf$tcn.em == 1] = 0

ST60_out = list(segs = ST60_fit$cncf,
                snps = ST60_post$jointseg,
                purity = ST60_fit$purity,
                ploidy = ST60_fit$ploidy,
                dipLogR = ST60_fit$dipLogR)

ST60_IGV_seg = facetsSuite::format_igv_seg(facets_output = ST60_out, sample_id = 'ST60')

write.table(ST60_IGV_seg, file = 'ST60.seg', row.names = F, quote = F, sep = '\t')

## GMM on Chromosome Y
ST60_local = ST60_post$jointseg[which(ST60_post$jointseg$chrom == 24), ]
GMM.processed = mixtools::normalmixEM(x = ST60_local$cnlr, k = 2)
GMM.lambda = GMM.processed$lambda
GMM.mu = GMM.processed$mu

# make plot
plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}

cnlr = ST60_local$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio') +
  #     title = 'WES Cn-LogR Distribution across the Y-chromosome',
  #    subtitle = '51.1% of markers belong to \'red\' distribution with mu = -3.25') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))
WES.mix.example


##-------------------------------------------------------------------
ST88 = facetsY::readSnpMatrix(filename = 'ST88.gz')
ST88Y = ST88[which(ST88$Chromosome == 'Y' & ST88$NOR.DP >= 35), ]
ST88.x = log2(ST88Y$TUM.DP / ST88Y$NOR.DP)
median(ST88.x)

# analysis pipeline
ST88_facetsSuite = facetsSuite::run_facets(read_counts = ST88, 
                                           cval = 200, 
                                           genome = 'hg19')
facetsSuite::cnlr_plot(ST88_facetsSuite)
ST88_QC = facets_fit_qc(facets_output = ST88_facetsSuite)

ST88_pre = facetsY::preProcSample(rcmat = ST88,
                                  ndepth = 35,
                                  het.thresh = 0.25,
                                  snp.nbhd = 250,
                                  cval = 25,
                                  gbuild = 'hg19',
                                  hetscale = TRUE)

ST88_post = facetsY::procSample(x = ST88_pre,
                                cval = 200,
                                min.nhet = 15)

ST88_fit = facetsY::emcncf(ST88_post)
View(ST88_fit$cncf[which(ST88_fit$cncf$chrom == 24), ])

## sub modifications:
ST88_fit$cncf = cbind(ST88_fit$cncf, cf = ST88_post$out$cf, tcn = ST88_post$out$tcn, lcn = ST88_post$out$lcn)
ST88_fit$cncf$lcn[ST88_fit$cncf$tcn == 1] = 0
ST88_fit$cncf$lcn.em[ST88_fit$cncf$tcn.em == 1] = 0

ST88_out = list(segs = ST88_fit$cncf,
                snps = ST88_post$jointseg,
                purity = ST88_fit$purity,
                ploidy = ST88_fit$ploidy,
                dipLogR = ST88_fit$dipLogR)

ST88_IGV_seg = facetsSuite::format_igv_seg(facets_output = ST88_out, sample_id = 'ST88')

write.table(ST88_IGV_seg, file = 'ST88.seg', row.names = F, quote = F, sep = '\t')

## GMM on Chromosome Y
ST88_local = ST88_post$jointseg[which(ST88_post$jointseg$chrom == 24), ]
GMM.processed = mixtools::normalmixEM(x = ST88_local$cnlr, k = 2)
GMM.lambda = GMM.processed$lambda
GMM.mu = GMM.processed$mu

# make plot
plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}

cnlr = ST88_local$cnlr
mixmdl = mixtools::normalmixEM(cnlr, k = 2)

#' IMPACT example
WES.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(y = "Density", x = 'Copy Number Log Ratio') +
  #     title = 'WES Cn-LogR Distribution across the Y-chromosome',
  #    subtitle = '51.1% of markers belong to \'red\' distribution with mu = -3.25') +
  theme_Y() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 16))
WES.mix.example






