##----------------+
## IGV representation of 
## Papillary and Hurthle cell
## tumors; thyroid tumors; 
## vast differences in LOY rates
##----------------+
##
## start: 02/02/23
## chris-kreitzer
## 

clean()
gc()
setwd('~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/IGV-like_plot.R')


cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
IGV = read.csv('Data/04_Loss/010523/IGV_out.seg', sep = '\t')


papillary_IGV = IGV[which(IGV$ID %in% cohort$SAMPLE_ID[which(cohort$CANCER_TYPE_DETAILED_ritika == 'Papillary Thyroid Cancer (THPA)')] &
                            IGV$chrom %in% c(1:23)), ]
hurthl_IGV = IGV[which(IGV$ID %in% cohort$SAMPLE_ID[which(cohort$CANCER_TYPE_DETAILED_ritika == 'Hurthle Cell Thyroid Cancer (THHC)')] &
                         IGV$chrom %in% c(1:23)), ]

papillary_plot = plot_segmentation(segs = papillary_IGV, plotX = F, cap_log_ratios = 2, return_object = F)
papillary_plot = papillary_plot + theme(legend.position = 'none') + labs(title = 'Papillary Thyroid Cancer (THPA)')
hurthle_plot = plot_segmentation(segs = hurthl_IGV, plotX = F, cap_log_ratios = 2)
hurthle_plot = hurthle_plot + labs(title = 'Hurthle Cell Thyroid Cancer (THHC)')

thyroid_cancer = papillary_plot / hurthle_plot
ggsave_golden(filename = 'Figures_original/ThyroidCancer_Segmentation.pdf', plot = thyroid_cancer, width = 12)


##----------------+
## SCNA for MSI vs MSS
##----------------+
MSI_SCNA = plot_segmentation(segs = IMPACT_instable[which(IMPACT_instable$chrom %in% c(1:22)), ], cap_log_ratios = 3)
MSI_SCNA = MSI_SCNA + theme(legend.position = 'none') + labs(title = 'MSI (n=299)')
MSS_SCNA = plot_segmentation(segs = IMPACT_stable[which(IMPACT_stable$chrom %in% c(1:22)), ], cap_log_ratios = 3)
MSS_SCNA = MSS_SCNA + labs(title = 'MSS (n=11,876)')

MSI_SCNA = MSI_SCNA / MSS_SCNA
MSI_SCNA
ggsave_golden(filename = 'Figures_original/MSI_SCNA.pdf', plot = MSI_SCNA, width = 14)


#' out