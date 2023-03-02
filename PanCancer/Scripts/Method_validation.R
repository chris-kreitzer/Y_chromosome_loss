##----------------+
## Confirm IMPACT chromosome Y 
## calls with WES sequenced
## samples; using CnLR as surrogate
## for binary Y-chromosome calls;
## - further check FacetsY: the method
## - on KIRC TCGA samples
## - confirmed via RNA-Seq experiments
## 
## - PanCancer mRNA expression overview
##    - take the calls from the preprint and compare
##    with mRNA expression of chrom. Y genes
##----------------+ 
## 
## Start (revision): 08/30/2021
## revision: 07/26/2022
## revision: 08/21/2022
## revision: 12/16/2022
## revision: 01/30/2023
## revision: 01/31/23
## revision: 02/28/2023
## revision: 03/02/2023
##
## chris-kreitzer

clean()
gc()
.rs.restartR()
set.seed(100)
setwd('~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')


## Input
WES_cnLR = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/WES/WES_Cnlr_Y.txt', sep = '\t')
IMPACT_cnLR = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/04_Loss/010523/Cnlr_Y.txt', sep = '\t')
WES.IMPACT.ids_df = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/WES_cBio_ID.match.txt', sep = '\t')



##----------------+
## median LRR for WES 
## chromosome Y
##-----------------+
WES_cnLR$mLRR = NA
for(i in unique(WES_cnLR$ID)){
  print(i)
  if(length(WES_cnLR$ID[which(WES_cnLR$ID == i)]) >= 10){
    WES_cnLR$mLRR[which(WES_cnLR$ID == i)] = median(WES_cnLR$cnlr[which(WES_cnLR$ID == i)], na.rm = T)
  } else {
    WES_cnLR$mLRR[which(WES_cnLR$ID == i)] = 0
  }
}


##----------------+
## median LRR for IMPACT
##----------------+
IMPACT_cnLR$mLRR = NA
for(i in unique(IMPACT_cnLR$ID)){
  print(i)
  if(length(IMPACT_cnLR$ID[which(IMPACT_cnLR$ID == i)]) >= 10){
    median.LRR = median(IMPACT_cnLR$cnlr[which(IMPACT_cnLR$ID == i)])
    IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == i)] = median.LRR
  } else {
    median.LRR = 0
    IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == i)] = median.LRR
  }
}


##----------------+
## merge the matched
## validation cohort
## - add further info
## - purity, QC
##----------------+
WES_cnLR$ID = substr(x = WES_cnLR$ID, start = 3, stop = 17)
Seq_merged = WES.IMPACT.ids_df[which(WES.IMPACT.ids_df$CMO_Sample_ID %in% unique(WES_cnLR$ID)), c('DMP_Sample_ID', 'CMO_Sample_ID', 'Sex')]

Seq_merged$mLRR_IMPACT = NA
Seq_merged$mLRR_WES = NA

for(i in 1:nrow(Seq_merged)){
  print(i)
  if(Seq_merged$DMP_Sample_ID[i] %in% IMPACT_cnLR$ID){
    Seq_merged$mLRR_IMPACT[i] = unique(IMPACT_cnLR$mLRR[which(IMPACT_cnLR$ID == Seq_merged$DMP_Sample_ID[i])])
  } else {
    Seq_merged$mLRR_IMPACT[i] = NA
  }
}

for(i in 1:nrow(Seq_merged)){
  print(i)
  if(Seq_merged$CMO_Sample_ID[i] %in% WES_cnLR$ID){
    Seq_merged$mLRR_WES[i] = unique(WES_cnLR$mLRR[which(WES_cnLR$ID == Seq_merged$CMO_Sample_ID[i])])
  } else {
    Seq_merged$mLRR_WES[i] = NA
  }
}


##-------
## add. annotation
##-------
IMPACT_cna = read.csv('Data/04_Loss/010523/CopyNumberStates.txt', sep = '\t')
WES_cna = read.csv('Data/02_Method_Validation/WES/WES_CopyNumberStates.txt', sep = '\t')
WES_cna$id = substr(x = WES_cna$id, start = 3, stop = 17)

Seq_merged = merge(Seq_merged, IMPACT_cna[, c('id', 'QC', 'purity', 'ploidy')], by.x = 'DMP_Sample_ID', by.y = 'id', all.x = T)
Seq_merged = Seq_merged[!duplicated(Seq_merged), ]
colnames(Seq_merged)[6:8] = c('IMPACT_QC', 'IMPACT_purity', 'IMPACT_ploidy')

Seq_merged = merge(Seq_merged, WES_cna[,c('id', 'QC', 'purity', 'ploidy')], by.x = 'CMO_Sample_ID', by.y = 'id', all.x = T)
Seq_merged = Seq_merged[!duplicated(Seq_merged), ]
colnames(Seq_merged)[9:11] = c('WES_QC', 'WES_purity', 'WES_ploidy')


Seq_merged = Seq_merged[!with(Seq_merged, is.na(mLRR_IMPACT) | is.na(mLRR_WES)), ]
Seq_merged = Seq_merged[which(Seq_merged$IMPACT_QC == TRUE), ]
Seq_merged = Seq_merged[which(Seq_merged$WES_QC == TRUE), ]
write.table(Seq_merged, file = '~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/ValidationCohort_out_n769.txt', sep = '\t', quote = F, row.names = F)




##----------------+
## Visualization
## - first mLRR of both cohorts
## - then purity
##----------------+
Y_continious_loss = ggplot(Seq_merged, 
                           aes(x = mLRR_WES, 
                               y = mLRR_IMPACT)) +
  geom_point() +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed',
              linewidth = 0.35,
              color = 'grey35') +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(-5, 2.5)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(-5, 2.5)) +
  theme_std(base_size = 14, base_line_size = 0.5) +
  annotate('text',
           x = -3.8,
           y = 1.5,
           label = paste0('r = ', round(cor.test(Seq_merged$mLRR_IMPACT, Seq_merged$mLRR_WES)[[4]][[1]], 2), 
                          '\np < 2.2e-16'),
           hjust = 0, 
           vjust = 0,
           family = 'ArialMT',
           size = 6) +
  
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, linewidth = 2),
        panel.background = element_rect(fill = 'white')) +
  
  labs(x = 'WES-recapture [median Cnlr]',
       y = 'MSK-IMPACT [median Cnlr]',
       title = paste0('Y-chromosome; n = ', length(unique(Seq_merged$CMO_Sample_ID))))


Y_continious_loss


##----------------+
## Purity comparison
##----------------+
Purity_plot = ggplot(Seq_merged, 
                     aes(x = WES_purity, 
                         y = IMPACT_purity)) +
  geom_point() +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed',
              linewidth = 0.35,
              color = 'grey35') +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 1)) +
  theme_std(base_size = 14, base_line_size = 0.5) +
  annotate('text',
           x = 0.15,
           y = 0.85,
           label = paste0('r = ', round(cor.test(Seq_merged$IMPACT_purity, Seq_merged$WES_purity)[[4]][[1]], 2), 
                          '\np < 2.2e-16'),
           hjust = 0, 
           vjust = 0,
           family = 'ArialMT',
           size = 6) +
  
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, linewidth = 2),
        panel.background = element_rect(fill = 'white')) +
  
  labs(x = 'purity WES-recapture',
       y = 'purity MSK-IMPACT',
       title = paste0('Y-chromosome; n = ', length(unique(Seq_merged$CMO_Sample_ID))))



ggsave_golden(filename = 'Figures_original/Method_validation_mCnLR.pdf', plot = Y_continious_loss, width = 8)
ggsave_golden(filename = 'Figures_original/Method_Purity_validation.pdf', plot = Purity_plot, width = 8)


method_validation = Y_continious_loss + Purity_plot
ggsave_golden(filename = '~/Documents/MSKCC/10_MasterThesis/Figures_original/Method_validation_MSK.pdf', plot = method_validation, width = 12)



##----------------+
## SAMPLE OUTLIER:
##----------------+
# P-0019114-T01-IM6
# P-0037069-T01-IM6
# P-0047344-T01-IM6
##----------------+


##----------------+
## Method validation on 
## TCGA, KIRC cases;
## validation through RNA-seq
##----------------+
TCGA_cohort = read.csv('Data/02_Method_Validation/TCGA/TCGA_all_CNAfits_out.txt', sep = '\t')
KIRC_cna = TCGA_cohort[which(TCGA_cohort$type == 'KIRC'), ]
mRNA = read.csv('Data/02_Method_Validation/TCGA/KIRC_TCGA_mRNA.txt', sep = '\t')
mRNA$SAMPLE_ID = substr(x = mRNA$SAMPLE_ID, start = 1, stop = 12)
mRNA$STUDY_ID = NULL

TCGA_KIRC = merge(KIRC_cna[,c('bcr_patient_barcode', 'Y_call', 'ploidy', 'purity')],
                  mRNA[,c('SAMPLE_ID', 'KDM5D', 'EIF1AX', 'DDX3Y', 'UTY', 'VHL', 'KDM5C')],
                  by.x = 'bcr_patient_barcode', by.y = 'SAMPLE_ID', all.x = T)

TCGA_KIRC = TCGA_KIRC[complete.cases(TCGA_KIRC), ]
TCGA_KIRC$KDM5D = log2(TCGA_KIRC$KDM5D+1)
TCGA_KIRC$KDM5C = log2(TCGA_KIRC$KDM5C+1)
TCGA_KIRC$EIF1AX = log2(TCGA_KIRC$EIF1AX+1)
TCGA_KIRC$DDX3Y = log2(TCGA_KIRC$DDX3Y+1)
TCGA_KIRC$UTY = log2(TCGA_KIRC$UTY+1)
TCGA_KIRC$VHL = log2(TCGA_KIRC$VHL+1)


##-------
## Visualization
##-------
vhl = ggplot(TCGA_KIRC, aes(y = VHL, x = Y_call)) +
  geom_quasirandom(width = 0.25, alpha = 0.45) +
  stat_summary(fun.y = "median", geom = "crossbar", size = 0.5, color = 'black', width = 0.35) +
  stat_summary(fun.y = "median", geom = "point", size = 3, color = 'red') +
  stat_compare_means(label.x = 1.3,
                     label.y = 10.5) +
  scale_x_discrete(labels = c('chrY wt', 'LOY')) +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(6, 11)) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 10),
        panel.border = element_rect(fill = NA, linewidth = 2)) +
  labs(x = '', y = 'mRNA expression (log2)', title = 'von Hippel-Lindau tumor suppressor, VHL')
  
#' KDM5C
kdm5c = ggplot(TCGA_KIRC, aes(y = KDM5C, x = Y_call)) +
  geom_quasirandom(width = 0.25, alpha = 0.45) +
  stat_summary(fun.y = "median", geom = "crossbar", size = 0.5, color = 'black', width = 0.35) +
  stat_summary(fun.y = "median", geom = "point", size = 3, color = 'red') +
  stat_compare_means(label.x = 1.3,
                     label.y = 12.5) +
  scale_x_discrete(labels = c('chrY wt', 'LOY')) +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(8, 13)) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 10),
        panel.border = element_rect(fill = NA, linewidth = 2)) +
  labs(x = '', y = '', title = 'Lysine Demethylase 5C, KDM5C')


##-------
## chrY genes:
kdm5d = ggplot(TCGA_KIRC, aes(y = KDM5D, x = Y_call)) +
  geom_quasirandom(width = 0.25, alpha = 0.45) +
  stat_summary(fun.y = "median", geom = "crossbar", size = 0.5, color = 'black', width = 0.35) +
  stat_summary(fun.y = "median", geom = "point", size = 3, color = 'red') +
  stat_compare_means(label.x = 1.3,
                     label.y = 12.5) +
  scale_x_discrete(labels = c('chrY wt', 'LOY')) +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(6, 13)) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 10),
        panel.border = element_rect(fill = NA, linewidth = 2)) +
  labs(x = '', y = '', title = 'Lysine Demethylase 5D, KDM5D')


ddx3y = ggplot(TCGA_KIRC, aes(y = DDX3Y, x = Y_call)) +
  geom_quasirandom(width = 0.25, alpha = 0.45) +
  stat_summary(fun.y = "median", geom = "crossbar", size = 0.5, color = 'black', width = 0.35) +
  stat_summary(fun.y = "median", geom = "point", size = 3, color = 'red') +
  stat_compare_means(label.x = 1.3,
                     label.y = 12.5) +
  scale_x_discrete(labels = c('chrY wt', 'LOY')) +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(6, 13)) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 10),
        panel.border = element_rect(fill = NA, linewidth = 2)) +
  labs(x = '', y = '', title = 'DEAD-Box Helicase 3 Y-Linked, DDX3Y')


uty = ggplot(TCGA_KIRC, aes(y = UTY, x = Y_call)) +
  geom_quasirandom(width = 0.25, alpha = 0.45) +
  stat_summary(fun.y = "median", geom = "crossbar", size = 0.5, color = 'black', width = 0.35) +
  stat_summary(fun.y = "median", geom = "point", size = 3, color = 'red') +
  stat_compare_means(label.x = 1.3,
                     label.y = 10.5) +
  scale_x_discrete(labels = c('chrY wt', 'LOY')) +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(6, 11)) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 10),
        panel.border = element_rect(fill = NA, linewidth = 2)) +
  labs(x = '', y = '', title = 'Ubiquitously Transcribed Tetratricopeptide\nRepeat Containing, Y-Linked, UTY')

grid_p = vhl + kdm5c + kdm5d + ddx3y + uty + plot_layout(ncol = 5)
ggsave_golden(filename = 'Figures_original/TCGA_mRNA_confirmation.pdf', plot = grid_p, width = 19)




##----------------+
## mRNA expression of
## chromosome Y genes depending 
## on tissue and chromosome Y status
## - taking LOY calls from the prePrint
##----------------+
tcga_loy = read.csv('Data/02_Method_Validation/TCGA/LOY_PanCancerTCGA.txt', sep = '\t')
tcga_loy$case_id = NULL
tcga_loy$X...tcn.em.means.tcn.are.from.estimated.tcn.by.FACETS..tcn.tcga.means.tcn.is.corrected.with.purity.and.ploidy.from.TCGA = NULL
tcga_loy$Gender = NULL
tcga_loy$Y_status_by_expression = NULL
tcga_loy = tcga_loy[which(tcga_loy$Y_status_source == 'tcn.em'), ]
tcga_loy$Genome_doublings_tcga = NULL
tcga_loy$control_tcga_sample_id = NULL
mrna = read.csv('Data/02_Method_Validation/TCGA/mRNA_Expression_PanCancerTCGA.txt', sep = '\t')
mrna$STUDY_ID = NULL

tcga_cohort = merge(tcga_loy, mrna, by.x = 'case_tcga_sample_id', by.y = 'SAMPLE_ID', all.x = T)

tcga_summary = data.frame()
for(i in unique(tcga_cohort$Cohort)){
  data_tumor = tcga_cohort[which(tcga_cohort$Cohort == i), ]
  for(j in c('ATM', 'AMELY', 'DDX3Y', 'KDM5D', 'RPS4Y1', 'RPS4Y2', 'SRY', 'UTY', 'ZFY')){
    gene_wt = log2(data_tumor[which(data_tumor$Y_status == 'WT'), j] + 1)
    gene_wt = ifelse(is.na(gene_wt), 0 , gene_wt)
    gene_loss = log2(data_tumor[which(data_tumor$Y_status == 'LOY'), j] + 1)
    gene_loss = ifelse(is.na(gene_loss), 0 , gene_loss)
    gene_gain = log2(data_tumor[which(data_tumor$Y_status == 'Gain'), j] + 1)
    gene_gain = ifelse(is.na(gene_gain), 0 , gene_gain)
    
    wt_loss_diff = median(gene_wt, na.rm = T) - median(gene_loss, na.rm = T)
    wt_gain_diff = median(gene_wt, na.rm = T) - median(gene_gain, na.rm = T)
    
    out = data.frame(cohort = i,
                     gene = j,
                     median_wt = median(gene_wt),
                     median_loss = median(gene_loss),
                     median_gain = median(gene_gain),
                     stat = wilcox.test(gene_wt, gene_loss)$p.value,
                     difference_wt_loss = wt_loss_diff,
                     difference_wt_gain = wt_gain_diff)
    tcga_summary = rbind(tcga_summary, out)
  }
}

tcga_summary$adj = p.adjust(tcga_summary$stat, method = 'BH')
tcga_summary$plot = ifelse(tcga_summary$adj < 0.01, 'plot', 'not_plot')
tcga_summary$plot = ifelse(tcga_summary$difference_wt_loss == 0, 'not_plot', tcga_summary$plot)

tcga_mrna = ggplot(tcga_summary, aes(x = reorder(cohort, difference_wt_loss), 
                         y = difference_wt_loss, color = plot)) +
  geom_point(size = 4) +
  scale_color_manual(values = c('not_plot' = 'grey',
                                'plot' = 'red'),
                     name = '',
                     labels = c('ns', 'Q<0.01')) +
  coord_flip() +
  theme_minimal() +
  panel_border(color = 'black') +
  theme(axis.text = element_text(size = 12, color = 'black'),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12, colour = 'black')) +
  facet_grid(~gene) +
  labs(x = '', y = 'mRNA expression difference [wildtype VS LOY]')

ggsave_golden(filename = 'Figures_original/TCGA_mRNA_correlation.pdf', plot = tcga_mrna, width = 12)


#' out