## PanCancer analysis of Y-chromosome loss:
## 
## Run FacetsY on TCGA-PRAD samples, and show that overall segmentation is comparable to Affymetrix 6.0 data
## I consider Affymetrix 6.0. as Golden-standard (independent validation)
## 
## Start: 04/12/2021


rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')
source('Y_chromosome_loss/genome.data.R')
source('Y_chromosome_loss/Plotting_theme.R')

## Libraries:
library(data.table)


## Input:
WES.TCGA.PRAD.segmentation_df = read.csv('Data_out/WES.TCGA.Prostate.processed.out.txt', sep = '\t')
WES.TCGA.PRAD.IGV_df = read.csv('Data_out/WES.TCGA.prostate.IGV.seg', sep = '\t')
Affy.TCGA_df = read.csv('~/Documents/MSKCC/CPNA_analysis/TCGA/RawData/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg', sep = '\t')


## Processing:
WES.TCGA.PRAD.segmentation_df$ID = substr(WES.TCGA.PRAD.segmentation_df$ID, start = 2, stop = nchar(WES.TCGA.PRAD.segmentation_df$ID))
WES.TCGA.PRAD.segmentation_df$ID = substr(WES.TCGA.PRAD.segmentation_df$ID, start = 1, stop = 16)
WES.TCGA.PRAD.IGV_df$ID = substr(WES.TCGA.PRAD.IGV_df$ID, start = 2, stop = 16)
colnames(WES.TCGA.PRAD.IGV_df) = c('Sample', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean')
WES.TCGA.PRAD.IGV_df$length = WES.TCGA.PRAD.IGV_df$End - WES.TCGA.PRAD.IGV_df$Start

#' fetch TCGA-PRAD samples from Affymetrix data.frame
Affy.TCGA_df$ID = substr(Affy.TCGA_df$Sample, start = 14, stop = 15)
Affy.TCGA_df$Sample = substr(Affy.TCGA_df$Sample, start = 1, stop = 16)

#' sample identifier: 01, 02, ... 09 {are tumor samples}
Affy.TCGA_df = Affy.TCGA_df[which(Affy.TCGA_df$ID %in% c('01', '02', '03', '05', '06')) ,]

#' select Affydata from TCGA-PRAD samples (Cell 2015)
Affy.TCGA.PRAD = data.frame()
for(i in unique(WES.TCGA.PRAD.segmentation_df$ID)){
  if(any(grepl(pattern = i, x = unique(Affy.TCGA_df$Sample)))){
    id = unique(grep(pattern = i, x = Affy.TCGA_df$Sample, value = T))
    sub = Affy.TCGA_df[which(Affy.TCGA_df$Sample == id), ]
  } else next
  Affy.TCGA.PRAD = rbind(Affy.TCGA.PRAD, sub)
}

Affy.TCGA.PRAD$ID = NULL
row.names(Affy.TCGA.PRAD) = NULL
Affy.TCGA.PRAD$length = Affy.TCGA.PRAD$End - Affy.TCGA.PRAD$Start
Affy.TCGA.PRAD = merge(Affy.TCGA.PRAD, hg19[, c('chrom', 'size', 'centromere')], 
                       by.x = 'Chromosome', by.y = 'chrom', all.x = T)

# write.table(Affy.TCGA.PRAD, 
#             file = 'Data_out/Affy.TCGA.PRAD.seg', 
#             sep = '\t', 
#             row.names = F, 
#             quote = F)

WES.TCGA.PRAD.IGV_df = merge(WES.TCGA.PRAD.IGV_df, hg19[,c('chrom', 'size', 'centromere')],
                             by.x = 'Chromosome', by.y = 'chrom', all.x = T)


## Analysis:
#' calculate arm-level alteration in WES-TCGA-PRAD cohort (run with FacetsY):
#' >80% of the arm with abs(log2R) > 0.2
#' segments spanning centromere are excluded for now
Data.for.Analysis = Affy.TCGA.PRAD
Data.summary = Data.for.Analysis %>%
  rowwise() %>%
  filter(abs(Segment_Mean) > 0.2) %>%
  mutate(
    arm = case_when(
      Start < centromere & End <= centromere ~ 'p',
      Start >= centromere ~ 'q',
      length > size * 0.8 ~ 'p;q',
      #Start < centromere & End > centromere & length > centromere * 0.8 ~ 'p',
      #Start < centromere & End > centromere & length > (size - centromere) * 0.8 ~ 'q'
    )
  ) %>%
  filter(!is.na(arm))

#' subset whole chromosome arms:
Data.summary.WCA = Data.summary[which(Data.summary$arm == 'p;q'), ]

WCA.alteration = data.frame()
for(i in unique(Data.summary.WCA$Sample)){
  sub = Data.summary.WCA[which(Data.summary.WCA$Sample == i), ]
  for(j in unique(sub$Chromosome)){
    out = data.frame(ID = i,
                     Chromosome = j,
                     Arm = c('p', 'q'),
                     Call = c(1, 1))
    WCA.alteration = rbind(WCA.alteration, out)
  }
}

#' analyse regular arms
Data.summary.arm = Data.summary[which(Data.summary$arm != 'p;q'), ]

arms.out = data.frame()
for(i in unique(Data.summary.arm$Sample)){
  data.sub = Data.summary.arm[which(Data.summary.arm$Sample == i), ]
  for(j in unique(data.sub$Chromosome)){
    for(arm in unique(data.sub$arm)){
      length.arm = sum(data.sub$length[which(data.sub$Chromosome == j & 
                                               data.sub$arm == arm)])
      length.arm.reference = Chromosome.arms.hg19.out$length[which(Chromosome.arms.hg19.out$Chromosome == j &
                                                                     Chromosome.arms.hg19.out$Arm == arm)]
      arm.fraction = length.arm / length.arm.reference
      arm.call = ifelse(arm.fraction >= 0.8, 1, 0)
      
      out = data.frame(ID = i,
                       Chromosome = j,
                       Arm = arm,
                       Call = arm.call)
      arms.out = rbind(arms.out, out)
    }
  }
}

#' merge data
Merged.alterations_df = rbind(arms.out, WCA.alteration)
Merged.alterations_df$arms.merged = paste0(Merged.alterations_df$Chromosome, Merged.alterations_df$Arm)

method = 'AffyMetrix'

Merged.summary.stats = data.frame()
for(i in unique(Merged.alterations_df$arms.merged)){
  freq = sum(Merged.alterations_df$Call[which(Merged.alterations_df$arms.merged == i)])
  freq.rel = freq / length(unique(Merged.alterations_df$ID))
  out = data.frame(arm = i,
                   freq = freq.rel * 100,
                   strategy = method)
  Merged.summary.stats = rbind(Merged.summary.stats, out)
}

#' run same logic with FacetsY
Affy.TCGA.PRAD_df = Merged.summary.stats
WES.TCGA.FACETSY_df = Merged.summary.stats

Methods.combined_df = merge(WES.TCGA.FACETSY_df,
                            Affy.TCGA.PRAD_df, 
                            by.x = 'arm',
                            by.y = 'arm',
                            all.x = T)

Method.comparison.plot = ggplot(Methods.combined_df,
       aes(x = freq.x, y = freq.y)) + geom_point() +
  geom_text(aes(label = arm), vjust = -0.5) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     breaks = c(0, 5, 10, 15, 20, 25),
                     limits = c(0, 27)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     breaks = c(0, 5, 10, 15, 20, 25),
                     limits = c(0, 27)) +
  theme_Y(base_size = 12) +
  theme(aspect.ratio = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey55', size = 0.2) +
  labs(x = 'FacetsY', y = 'Affymetrix 6.0', title = 'FacetsY decipher the overall SCNA aberration profile in TCGA-PRAD samples', 
       subtitle = 'FacetsY was run on n=333 WES TCGA-PRAD samples (default parameters) and recovers chromosome-arm\nalteration frequencies, similar to those called by Affymetrix 6.0')
       
ggsave_golden(filename = 'Figures/Method.comparison.pdf', plot = Method.comparison.plot, width = 16)  
