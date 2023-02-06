##----------------+
## GLM model (gene-wise) and 
## Fisher's exact test
## for mutations in genes
## associated with LOY;
## general inspection of mutations
## and loy;
##----------------+
##
## start: 11/11/2022
## revision: 11/20/2022
## revision: 01/17/2023
## revision: 01/18/2023
## revision: 01/23/2023
## 
## chris-kreitzer


## TODO:
## - QC T/F samples included
## - just work with SOMATIC mutations?
##   (currently also GERMLINE variants included)
## - mutation cutoff
## - FUSIONS; old cohort (Bastien)
## - EVENT: only for complete_loss
##   (partial, and relative loss not considered)



clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')
library(patchwork)
library(data.table)
library(ggpubr)


cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
cohort_samples = unique(cohort$SAMPLE_ID)
oncoKB_mutations = read.csv('Data/05_Mutation/data_mutations_extended.oncokb.txt', sep = '\t')
oncoKB_cna = data.table::fread('Data/05_Mutation/data_CNA.txt', sep = '\t')
oncoKB_sv = read.csv('Data/05_Mutation/data_sv.oncokb.txt', sep = '\t')
data_gam = list()


##----------------+
## set-up mutation matrix
## - all mutations are INCLUDED;
##----------------+
mutation_matrix = setNames(data.frame(matrix(ncol = length(unique(oncoKB_mutations$Hugo_Symbol)), nrow = 0)), unique(oncoKB_mutations$Hugo_Symbol))
for(id in 1:length(cohort_samples)){
  print(id)
  data.mut.sub = oncoKB_mutations[which(oncoKB_mutations$Tumor_Sample_Barcode == cohort_samples[id]), ]
  if(nrow(data.mut.sub) != 0){
    for(j in unique(data.mut.sub$Hugo_Symbol)){
      if(j %in% colnames(mutation_matrix)){
        mutation_matrix[id, j] = 1
      } else {
        mutation_matrix[id, j] = 0
      }
    }
    mutation_matrix[id, 'Sample.ID'] = cohort_samples[id]
  } else {
    mutation_matrix[id, ] = 0
    mutation_matrix[id, 'Sample.ID'] = cohort_samples[id]
}

  }
row.names(mutation_matrix) = mutation_matrix$Sample.ID
row.names(mutation_matrix) = NULL
mutation_matrix[is.na(mutation_matrix)] = 0
mutation_matrix = as.data.frame.matrix(mutation_matrix)
colnames(mutation_matrix) = paste0(colnames(mutation_matrix), '_mut')
colnames(mutation_matrix)[ncol(mutation_matrix)] = 'sample'

data_gam = c(data_gam, list(mut = mutation_matrix))
saveRDS(object = data_gam, file = 'Data/05_Mutation/data_gam.rds')



##----------------+
## set-up CNA matrix
## - all CNAs are INCLUDED;
##----------------+
colnames(oncoKB_cna)[2:ncol(oncoKB_cna)] = gsub(pattern = '\\.', replacement = '-', x = colnames(oncoKB_cna)[2:ncol(oncoKB_cna)])
oncoKB_cna = t(oncoKB_cna)
colnames(oncoKB_cna) = oncoKB_cna[1, ]
oncoKB_cna = oncoKB_cna[-1, ]
oncoKB_cna = as.data.frame.matrix(oncoKB_cna)
oncoKB_cna$sample = row.names(oncoKB_cna)
oncoKB_cna = as.data.frame(apply(oncoKB_cna, 2, function(x) gsub("\\s+", "", x)))
oncoKB_cna[oncoKB_cna == '0.0'] = 0
oncoKB_cna[oncoKB_cna == '-2.0'] = -2
oncoKB_cna[oncoKB_cna == '2.0'] = 2
oncoKB_cna[oncoKB_cna == '-1.5'] = -2

##-------
## check if there
## are any more weird
## encodings;
##-------
df2 = as.vector(as.matrix(oncoKB_cna[,-ncol(oncoKB_cna)]))
unique(df2)

##-------
## alteration matrix;
## - first: Amplifications
## - second: Deletions
##-------
oncokb_amp = oncoKB_cna
oncokb_amp$sample = NULL
colnames(oncokb_amp) = paste0(colnames(oncokb_amp), '_Amplification')
oncokb_amp[oncokb_amp == '2'] = 1
oncokb_amp[oncokb_amp == '-2'] = 0

oncokb_del = oncoKB_cna
oncokb_del$sample = NULL
colnames(oncokb_del) = paste0(colnames(oncokb_del), '_Deletion')
oncokb_del[oncokb_del == '-2'] = 1
oncokb_del[oncokb_del == '2'] = 0

cna_matrix = merge(oncokb_amp, oncokb_del, by = 'row.names', all.x = T)
colnames(cna_matrix)[1] = 'sample'

row.names(cna_matrix) = cna_matrix$sample
cna_matrix$sample = NULL

cna_matrix = as.data.frame(lapply(cna_matrix, as.integer), row.names = row.names(cna_matrix))
cna_matrix$sample = row.names(cna_matrix)
row.names(cna_matrix) = NULL

data_gam = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = list(mut = data_gam$mut, cna = cna_matrix)
saveRDS(object = data_gam, file = 'Data/05_Mutation/data_gam.rds')



##----------------+
## set-up Fusion matrix
## - all Fusions are INCLUDED;
## - use Bastien's matrix for now - no update
##----------------+
gam_cna = load('Data/signedOut/1.0.genomic_data_all.Rdata')
data_fusion = data_gam$fusion
rm(data_gam, data_gam_onc, data_MAF, data_broad_genomic)
gc()

data_fusion = as.data.frame.matrix(data_fusion) 
data_fusion$sample = row.names(data_fusion)
row.names(data_fusion) = NULL
data_gam = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = list(mutations = data_gam$mut, CNAs = data_gam$cna, fusions = data_fusion)
saveRDS(object = data_gam, file = 'Data/05_Mutation/data_gam.rds')


##----------------+
## MERGE all three object;
## Mutations, CNAs, Fusions;
##----------------+
data_gam = readRDS('Data/05_Mutation/data_gam.rds')
muts = data_gam$mutations
cnas = data_gam$CNAs
fusions = data_gam$fusions
first = merge(muts, cnas, by = 'sample', all.x = T)
GAM_all = merge(first, fusions, by = 'sample', all.x = T)

data_gam = list(Mutations = data_gam$mutations, 
                CNAs = data_gam$CNAs,
                Fusions = data_gam$fusions,
                All = GAM_all)
saveRDS(object = data_gam, file = 'Data/05_Mutation/data_gam.rds')



##----------------+
## Genes to keep;
## frequently occurring mutations;
##----------------+
data_gam_raw = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = data_gam_raw$GAM_Analysis
alterations_out = data.frame()
for(i in 2:length(data_gam)){
  gene = colnames(data_gam)[i]
  total = nrow(data_gam)
  n_specific = sum(data_gam[, i], na.rm = T)
  out = data.frame(alteration = gene,
                   fraction = (n_specific / total) * 100)
  alterations_out = rbind(alterations_out, out)
}

##-------
## exclude genes:
## < 0.1% occurrence
##-------
genes_remove = alterations_out$alteration[which(alterations_out$fraction <= 0.1)]
genes_remove = genes_remove[!genes_remove %in% c('KDM5C_Deletion', 'VHL_Deletion', 'EIF1AX_Deletion',
                                                 'PHF6_mut', 'PHF6_Deletion', 'PHF6_Amplification',
                                                 'ZRSR2_Deletion')]

data_gam_analysis = data_gam[, !colnames(data_gam) %in% genes_remove]
data_gam = list(Mutations = data_gam_raw$Mutations,
                CNAs = data_gam_raw$CNAs,
                Fusions = data_gam_raw$Fusions,
                All = data_gam_raw$All,
                GAM_Analysis = data_gam_analysis)
saveRDS(object = data_gam, file = 'Data/05_Mutation/data_gam.rds')




##----------------+
## Prepare genomic/clinical 
## data table for gene-wise
## glm;
##----------------+
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
clinical_glm = cohort[,c('SAMPLE_ID', 'QC', 'ploidy', 'classification', 'purity', 'CANCER_TYPE', 'MSI_TYPE',
                         'MSI_SCORE', 'Age_Sequencing', 'IMPACT_TMB_SCORE', 'MUTATION_COUNT', 'SAMPLE_TYPE',
                         'genome_doubled', 'fraction_cna')]
clinical_glm$classification[which(clinical_glm$classification %in% c('complete_loss', 'partial_loss', 'relative_loss'))] = 1
clinical_glm$classification[which(clinical_glm$classification %in% c('wt', 'gain', 'gain_loss', 'partial_gain'))] = 0
clinical_glm = clinical_glm[!is.na(clinical_glm$classification), ]
colnames(clinical_glm) = c('sample', 'QC', 'ploidy', 'Y_call', 'purity', 'CANCER_TYPE', 'MSI_TYPE', 
                           'MSI_SCORE', 'Age', 'TMB', 'Mut_Count', 'SAMPLE_TYPE', 'WGD', 'FGA')



##----------------+
## TMB and LOY;
##----------------+
clinical_glm = clinical_glm[which(clinical_glm$CANCER_TYPE %in% ctypes_keep), ]
clinical_glm$Y_call = ifelse(clinical_glm$Y_call == '1', 'LOY', 'wt')

#' All together
TMB_LOY_All = ggplot(clinical_glm[!is.na(clinical_glm$TMB) & clinical_glm$TMB > 0, ], 
                     aes(x = Y_call, y = log10(TMB))) +
  geom_boxplot(aes(color = Y_call), position = position_dodge2(width = 0.75), 
               outlier.alpha = 0.25, outlier.stroke = 0.5) +
  scale_color_manual(values = c('wt' = '#3d4397',
                                'LOY' = '#a22231'),
                     name = '',
                     label = c('nonLOY', 'LOY')) +
  scale_y_continuous(limits = c(0 , 3)) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 2,
        axis.line = element_blank(),
        legend.position = 'none') +
  labs(x = '', y = 'Mutational burden (log10)') +
  panel_border(size = 2, color = 'black')


ggsave_golden(filename = 'Figures_original/TMB_LOY_all.pdf', plot = TMB_LOY_All, width = 6)


#' CancerType specific
TMB_LOY_plot = ggplot(clinical_glm[!is.na(clinical_glm$TMB) & clinical_glm$TMB > 0, ], 
                      aes(x = Y_call, y = log10(TMB))) +
  geom_boxplot(aes(color = Y_call), position = position_dodge2(width = 1), 
               outlier.alpha = 0.25, outlier.stroke = 0.5) +
  facet_wrap(~CANCER_TYPE, nrow = 1, scales = 'fixed') +
  stat_compare_means(comparisons = list(c('LOY', 'wt')),
                     method = 'wilcox.test', label = 'p.signif', hide.ns = T) +
  scale_color_manual(values = c('wt' = '#3d4397',
                                'LOY' = '#a22231'),
                     name = '',
                     label = c('wt', 'LOY')) +
  scale_y_continuous(limits = c(0 , 3)) +
  geom_vline(xintercept = 2.5, linetype = 'dashed', linewidth = 0.25) +
  theme_std(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(0.8), margin = margin(), angle = 90, hjust = 0),
        panel.spacing = unit(0, "pt"),
        legend.position = 'none',
        axis.text.x = element_blank()) +
  labs(x = '', y = 'Mutational burden (log10)')

TMB_LOY = TMB_LOY_All + TMB_LOY_plot        
ggsave_golden(filename = 'Figures_original/TMB_LOY.pdf', plot = TMB_LOY, width = 14)



##----------------+
## TP53 mutations;
##----------------+
data_gam = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = data_gam$GAM_Analysis
clinical_glm = merge(clinical_glm, data_gam[, c('sample', 'TP53_mut', 'TP53_Deletion')],
                     by.x = 'sample', by.y = 'sample', all.x = T)
clinical_glm$TP53_Deletion = as.numeric(as.integer(clinical_glm$TP53_Deletion))
TP53 = data.frame()
for(i in unique(clinical_glm$Y_call)){
  n = length(clinical_glm$sample[which(clinical_glm$Y_call == i)])
  n_muts = sum(clinical_glm$TP53_mut[which(clinical_glm$Y_call == i)], na.rm = T)
  n_dels = sum(clinical_glm$TP53_Deletion[which(clinical_glm$Y_call == i)], na.rm = T)
  n_alts = n_muts + n_dels
  out = data.frame(group = i,
                   n = n,
                   fraction_altered = (n_alts / n) * 100,
                   fraction_wt = ((n - n_alts) / n) * 100)
  TP53 = rbind(TP53, out)
}

TP53_muts = ggplot(TP53, aes(x = group, y = fraction_altered)) +
  geom_bar(stat = 'identity', color = 'black', fill = 'black') +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_discrete(label = c(paste0('LOY\n', '(n=', TP53$n[which(TP53$group == 'LOY')], ')'),
                             paste0('WT\n', '(n=', TP53$n[which(TP53$group == 'wt')], ')'))) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1) +
  labs(x = '', y = expression(paste('% male samples\nwith TP53 mutation')))


ggsave_golden(filename = 'Figures_original/TP53_mutations.pdf', plot = TP53_muts, width = 6)



##-------
## LOY and TP53 wild-type
##-------
samples_TP53wt_loy = clinical_glm$sample[which(clinical_glm$Y_call == 'LOY' & clinical_glm$TP53_Deletion == 0 & clinical_glm$TP53_mut == 0)]
samples_TP53mut_loy = clinical_glm$sample[which(clinical_glm$Y_call == 'LOY' & clinical_glm$TP53_Deletion == 1 | clinical_glm$TP53_mut == 1)]
cbio(samples_TP53wt_loy)
cbio(samples_TP53mut_loy)



##----------------+
## Fishers exact test
##----------------+
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
data_gam = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = data_gam$GAM_Analysis

cancer_type_list = unique(cohort$CANCER_TYPE)
cohort_samples = unique(cohort$SAMPLE_ID)

cohort$Y_call = ifelse(cohort$classification %in% c('complete_loss', 'relative_loss', 'partial_loss'), 'Y_chrom_loss', 'intact_Y_chrom')

##----------------+
## GOI: which genes to keep
##----------------+
GOI = unique(colnames(data_gam))
GOI = GOI[!GOI %in% 'sample']


##----------------+
## populate alteration
## matrix;
##----------------+
PanCancer = data.frame()
for(i in unique(GOI)){
  try({
    loss_n = table(cohort$Y_call)['Y_chrom_loss'][[1]]
    wt_n = table(cohort$Y_call)['intact_Y_chrom'][[1]]
    gene = i
    print(gene)
    cancer = 'PanCancer'
    gene_index = which(colnames(data_gam) == i)
    Muts = data_gam[, c(1, gene_index)]
    
    #' LOY with at least 1 mut
    LOY_oneMut = length(intersect(cohort$SAMPLE_ID[which(cohort$Y_call == 'Y_chrom_loss')],
                                  Muts$sample[which(Muts[,2] != 0)]))
    
    #' LOY with no Mutation
    LOY_wt = loss_n - LOY_oneMut
    
    #' WT tumor with at least 1 mutation
    WT_oneMut = length(intersect(cohort$SAMPLE_ID[which(cohort$Y_call == 'intact_Y_chrom')],
                                 Muts$sample[which(Muts[,2] != 0)]))
    
    #' WT tumor with no Mutation
    WT_noMut = wt_n - WT_oneMut
    
    out = data.frame(cohort = cancer,
                     gene = gene,
                     LOY_1Mut = LOY_oneMut,
                     LOY_wt = LOY_wt,
                     WT_oneMut = WT_oneMut,
                     WT_wt = WT_noMut)
    PanCancer = rbind(PanCancer, out)
  })
}



##----------------+
## Fisher Test on 
## PanCancer cohort
##----------------+
all_Fisher = PanCancer
all_Fisher = all_Fisher[which(all_Fisher$LOY_1Mut >= 5), ]
all_Fisher$one_sided_fisher_fisher_p.value = NA
for(i in 1:nrow(all_Fisher)){
  print(i)
  all_Fisher$one_sided_fisher_fisher_p.value[i] = fisher.test(matrix(c(all_Fisher$LOY_1Mut[i], 
                                                                       all_Fisher$WT_oneMut[i],
                                                                       all_Fisher$LOY_wt[i],
                                                                       all_Fisher$WT_wt[i]), ncol = 2), 
                                                              alternative = 'greater')$p.value
}

all_Fisher$FDR = p.adjust(all_Fisher$one_sided_fisher_fisher_p.value, method = 'fdr')
all_Fisher$log10p = -log10(all_Fisher$one_sided_fisher_fisher_p.value)
all_Fisher$FDR_pass = ifelse(all_Fisher$FDR < 0.01, 'plot', 'not')
LOY_enriched_all = all_Fisher[which(all_Fisher$gene %in% all_Fisher$gene[which(all_Fisher$LOY_1Mut > all_Fisher$WT_oneMut & all_Fisher$FDR_pass == 'plot')]), ]
View(LOY_enriched_all)
write.table(x = all_Fisher, file = 'Data/05_Mutation/011823/PanCancer_Fisher_enrichment.txt', sep = '\t', row.names = F, quote = F)


##----------------+
## PanCancer Fisher-test
## Visualization;
##----------------+
all_Fisher$FDR_pass = ifelse(all_Fisher$log10p > 10, 'plot', 'not')
pos = position_jitter(width = 0.3, seed = 2)

PanCancer_plot = ggplot(all_Fisher, 
                        aes(x = cohort, 
                            y = log10p, 
                            color = FDR_pass)) +
  geom_jitter(position = pos) +
  scale_color_manual(values = c('not' = 'grey75',
                                'plot' = 'red'),
                     name = '') +
  geom_text_repel(aes(label = ifelse(FDR_pass == 'plot', gene, '')), 
                  size = 2, max.overlaps = 400, position = pos) +
  scale_y_continuous(position = 'right', 
                     expand = c(0.01, 0)) +
  coord_flip() +
  theme_std(base_size = 14) +
  theme(axis.line = element_blank(),
        legend.position = 'none',
        axis.text.y = element_text(angle = 90, size = 16, hjust = 0.5),
        axis.ticks.y = element_blank()) +
  labs(x = '', y = '-log10(p-value)') +
  panel_border(size = 2, color = 'black')

ggsave_golden(filename = 'Figures_original/FISHER_PanCancer.pdf', plot = PanCancer_plot, width = 11)



##----------------+
## Fishers exact test;
## Per Cancer-Type basis;
## Fisher Test subsequent;
## +Visualization
##----------------+
all_cancer_gene = data.frame()
for(i in unique(cohort$CANCER_TYPE)){
  try({
    cancer = i
    print(i)
    loss_n = table(cohort$Y_call[which(cohort$CANCER_TYPE == i)])[['Y_chrom_loss']]
    wt_n = sum(table(cohort$Y_call[which(cohort$CANCER_TYPE == i)])) - loss_n
    
    for(j in unique(GOI)){
      gene = j
      print(gene)
      gene_index = which(colnames(data_gam) == j)
      Muts = data_gam[, c(1, gene_index)]
      
      #' LOY with at least 1 mut
      LOY_oneMut = length(intersect(cohort$SAMPLE_ID[which(cohort$Y_call == 'Y_chrom_loss' & cohort$CANCER_TYPE == i)],
                                    Muts$sample[which(Muts[,2] != 0)]))
      
      #' LOY with no Mutation
      LOY_wt = loss_n - LOY_oneMut
      
      #' WT tumor with at least 1 mutation
      WT_oneMut = length(intersect(cohort$SAMPLE_ID[which(cohort$Y_call == 'intact_Y_chrom' & cohort$CANCER_TYPE == i)],
                                   Muts$sample[which(Muts[,2] != 0)]))
      
      #' WT tumor with no Mutation
      WT_noMut = wt_n - WT_oneMut
      
      out = data.frame(cohort = cancer,
                       gene = gene,
                       LOY_1Mut = LOY_oneMut,
                       LOY_wt = LOY_wt,
                       WT_oneMut = WT_oneMut,
                       WT_wt = WT_noMut)
      all_cancer_gene = rbind(all_cancer_gene, out)
    }
  })
}

Fisher_gene = all_cancer_gene
Fisher_gene$one_sided_fisher_fisher_p.value = NA
for(i in 1:nrow(Fisher_gene)){
  print(i)
  Fisher_gene$one_sided_fisher_fisher_p.value[i] = fisher.test(matrix(c(Fisher_gene$LOY_1Mut[i], 
                                                                        Fisher_gene$WT_oneMut[i],
                                                                        Fisher_gene$LOY_wt[i],
                                                                        Fisher_gene$WT_wt[i]), ncol = 2), 
                                                               alternative = 'greater')$p.value
}

Fisher_gene$FDR = p.adjust(Fisher_gene$one_sided_fisher_fisher_p.value, method = 'fdr')
Fisher_gene$log10p = -log10(Fisher_gene$one_sided_fisher_fisher_p.value)
Fisher_gene$FDR_pass = ifelse(Fisher_gene$FDR < 0.05, 'plot', 'not')


##----------------+
## Plot the observations
##----------------+
Fisher_gene = Fisher_gene[which(Fisher_gene$cohort %in% ctypes_keep), ]
cancerTypes = ggplot(Fisher_gene, aes(x = cohort, 
                                      y = log10p, 
                                      color = FDR_pass)) +
  geom_jitter(position = pos) +
  scale_color_manual(values = c('not' = 'grey75',
                                'plot' = 'red'),
                     name = '') +
  geom_text_repel(aes(label = ifelse(FDR_pass == 'plot', gene, '')), 
                  color = 'black', 
                  size = 3, 
                  max.overlaps = 50,
                  position = pos) +
  coord_flip() +
  geom_vline(xintercept = seq(1.5, 55, 1), linetype = 'dashed', color = 'grey35', size = 0.4) +
  scale_y_continuous(position = 'right', 
                     expand = c(0.01, 0)) +
  theme_std(base_size = 14) +
  theme(panel.border = element_rect(fill = NA, linewidth = 2, color = 'black'),
        legend.position = 'none') +
  labs(x = '', y = '-log10(p-value)')


cancerTypes

##----------------+
## combine the plots
##----------------+
geneLevel_all = PanCancer_plot / cancerTypes + plot_layout(heights = c(0.15, 1))
ggsave_golden(filename = 'Figures_original/FISHER_geneLevel_all.pdf', plot = geneLevel_all, width = 18)





##----------------+
## run gene-wise glm
##----------------+
# Filter for gene then for cancer type
# Run logistic regression if:
## 1. gene has at least 1 mutant and at least 1 wild type
## 2. gene has alteration frequency of >= 0.03
clean()
gc()
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
clinical_glm = cohort[,c('SAMPLE_ID', 'QC', 'ploidy', 'classification', 'purity', 'CANCER_TYPE', 'MSI_TYPE',
                         'MSI_SCORE', 'Age_Sequencing', 'IMPACT_TMB_SCORE', 'MUTATION_COUNT', 'SAMPLE_TYPE',
                         'genome_doubled', 'fraction_cna')]
clinical_glm$classification[which(clinical_glm$classification %in% c('relative_loss', 'partial_loss','complete_loss'))] = 1
clinical_glm$classification[which(clinical_glm$classification %in% c('wt', 'gain', 'gain_loss', 'partial_gain'))] = 0
clinical_glm = clinical_glm[!is.na(clinical_glm$classification), ]
colnames(clinical_glm) = c('sample', 'QC', 'ploidy', 'Y_call', 'purity', 'CANCER_TYPE', 'MSI_TYPE', 
                           'MSI_SCORE', 'Age', 'TMB', 'Mut_Count', 'SAMPLE_TYPE', 'WGD', 'FGA')
data_gam = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = data_gam$GAM_Analysis


##----------------+
## Merging
##----------------+
clinical_glm = merge(clinical_glm, data_gam, by.x = 'sample', by.y = 'sample', all.x = T)
clinical_glm$QC = NULL
clinical_glm$TMB = NULL
clinical_glm$Mut_Count = NULL
clinical_glm$SAMPLE_TYPE = NULL
clinical_glm$MSI_SCORE = NULL
clinical_glm$Age = as.numeric(as.character(clinical_glm$Age))
clinical_glm$Y_call = as.integer(as.character(clinical_glm$Y_call))


##-------
## TP53 status fixed variable;
##-------
clinical_glm$TP53_Status = ifelse(clinical_glm$TP53_Deletion == 1, 1,
                                  ifelse(clinical_glm$TP53_fusion == 1, 1,
                                         ifelse(clinical_glm$TP53_mut == 1, 1,
                                                ifelse(clinical_glm$TP53_intragenic == 1, 1, 0))))

clinical_glm$TP53_Deletion = NULL
clinical_glm$TP53_mut = NULL
clinical_glm$TP53_fusion = NULL
clinical_glm$TP53_intragenic = NULL


##----------------+
## Gene-wise GLM;
##----------------+
# Filter for gene then for cancer type
# Run logistic regression if:
## 1. gene has at least 1 mutant and at least 1 wild type
## 2. gene has alteration frequency of >= 0.03
logistic_reg_fun = function(data_frame, gene, cancer_type){
  try({
    data_frame = data_frame[is.na(data_frame[, gene]) == F &
                              data_frame[, "CANCER_TYPE"] == cancer_type, ]
    if (length(which(data_frame[, gene] == 1)) == 0 | 
        length(which(data_frame[, gene] == 0)) == 0){
      log_results_df = data.frame(variable = "Not Tested",
                                  gene = gene, 
                                  cancer_type = cancer_type, 
                                  comments = "No Mutations in this gene") 
    } else if (length(which(data_frame[, gene] == 1)) / nrow(data_frame) < 0.005) {
      log_results_df = data.frame(variable = "Not Tested",
                                  gene = gene, 
                                  cancer_type = cancer_type, 
                                  comments = "Mutation frequency <3%") 
    } else { 
      formula = as.formula(paste0(gene, "~ Y_call + ploidy + FGA + purity + TP53_Status + MSI_TYPE"))
      log_results = glm(formula, data = data_frame, family = binomial)
      log_results_df = as.data.frame(summary(log_results)$coefficients)
      log_results_df$variable = row.names(log_results_df)
      log_results_df$gene = gene
      log_results_df$cancer_type = cancer_type
      colnames(log_results_df)[1:4] = c("estimate", "std_err", "z_value", "p_value")
      conf_df = as.data.frame(confint.default(log_results))
      conf_df$variable = row.names(conf_df)
      log_results_df = left_join(log_results_df, conf_df, by = "variable")
      log_results_df = log_results_df[c(5,7,6,1:4,8,9)]
      print(gene)
      print(cancer_type)
      return(log_results_df)
    }
  })
}

##----------------+
## Set gene list and 
## cancer type list
##----------------+
gene_list = colnames(clinical_glm)[10:ncol(clinical_glm)]
gene_list = gene_list[!gene_list %in% c('TP53_Status')]
# gene_list = gene_list[!gene_list %in% c('TP53', "HLA-A_Deletion", "HLA-B_Deletion", "NKX2-1_Amplification", "NKX3-1_Deletion")]
cancer_type_list = ctypes_keep

# Expand table for gene list and cancer list (get every combo)
gene_cancer_df = expand.grid(gene_list, cancer_type_list)
gene_cancer_df =  gene_cancer_df %>%
  mutate_if(is.factor, as.character)


# Run logistic regression
log_results_df = mapply(logistic_reg_fun,
                         gene = gene_cancer_df$Var1,
                         cancer_type = gene_cancer_df$Var2,
                         MoreArgs = list(data_frame = clinical_glm),
                         SIMPLIFY = FALSE)

for(i in 1:length(log_results_df)){
  if(class(log_results_df[[i]]) != 'data.frame'){
    print(i)
  } else next
}

log_results_df = Filter(function(x) length(x) > 1, log_results_df)
log_results_df = do.call("rbind.fill", log_results_df)
log_results_df$p_adj = p.adjust(log_results_df$p_value, method = "fdr")
log_results_df = log_results_df[!is.na(log_results_df$p_adj), ]
log_results_df = log_results_df[which(log_results_df$variable == 'Y_call'), ]


write.table(x = log_results_df, file = 'Data/05_Mutation/011823/GLM_gene_level_full_out.txt', sep = '\t', quote = F, row.names = F)


##----------------+
## Visualization;
##----------------+
gene_glm = log_results_df[which(log_results_df$p_adj <= 0.15 & 
                                  log_results_df$estimate > -15 & 
                                  log_results_df$estimate < 15), ]
length_genes = length(unique(gene_glm$gene))
pos = position_jitter(width = 0.2, seed = 2)

gene_glm_plot = ggplot(gene_glm, 
                       aes(x = reorder(gene, estimate),
                           y = estimate, 
                           label = cancer_type)) +
  geom_pointrange(aes(ymin = estimate - std_err,
                      ymax = estimate + std_err),
                  size = 0.75,
                  position = pos,
                  fatten = 4) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey35', size = 0.25) +
  geom_vline(xintercept = seq(1.5, length_genes, 1), linetype = 'dashed', color = 'grey35', size = 0.4) +
  geom_text_repel(aes(label = cancer_type), color = 'black', size = 4, position = pos) +
  scale_y_continuous(expand = c(0.1, 0),
                     limits = c(-2, 5),
                     sec.axis = dup_axis()) +
  coord_flip() +
  theme_std(base_size = 14) +
  theme(legend.position = 'none',
        axis.line = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = '', y = 'log ODDS') +
  panel_border(size = 2, color = 'black')

gene_glm_plot

ggsave_golden(filename = 'Figures_original/GeneLevel_CancerTypes_out.pdf', plot = gene_glm_plot, width = 12)


#' out