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
library(patchwork)
library(data.table)

cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
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
data_gam = data_gam_raw$All
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
clinical_glm = cohort[,c('SAMPLE_ID', 'QC', 'ploidy', 'classification', 'purity', 'CANCER_TYPE', 'MSI_TYPE',
                         'MSI_SCORE', 'Age_Sequencing', 'IMPACT_TMB_SCORE', 'MUTATION_COUNT', 'SAMPLE_TYPE',
                         'genome_doubled', 'fraction_cna')]
clinical_glm$classification[which(clinical_glm$classification == 'complete_loss')] = 1




clinical_glm$classification[which(clinical_glm$classification %in% c('loss', 'relative_loss'))] = 1
clinical_glm$classification[which(clinical_glm$classification %in% c('gain', 'wt', 'gain_loss'))] = 0
clinical_glm$classification[which(clinical_glm$sample == 'P-0002130-T02-IM3')] = 0
clinical_glm$classification[which(clinical_glm$sample == 'P-0013100-T01-IM5')] = 0
clinical_glm$classification[which(clinical_glm$sample == 'P-0028409-T01-IM6')] = 1
clinical_glm$classification[which(clinical_glm$sample == 'P-0064933-T02-IM7')] = 0
clinical_glm = merge(clinical_glm, clinical[,c('SAMPLE_ID', 'CANCER_TYPE', 'MSI_TYPE', 
                                               'MSI_SCORE', 'AGE_AT_WHICH_SEQUENCING_WAS_REPORTED_(YEARS)', 
                                               'IMPACT_TMB_SCORE', 'MUTATION_COUNT')],
                     by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
clinical_glm = merge(clinical_glm, cohort$IMPACT_cohort[,c('SAMPLE_ID', 'SAMPLE_TYPE')],
                     by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
clinical_glm = merge(clinical_glm, arm_level[,c('id', 'genome_doubled', 'fraction_cna')],
                     by.x = 'sample', by.y = 'id', all.x = T)
clinical_glm$classification = as.integer(as.character(clinical_glm$classification))
colnames(clinical_glm) = c('sample', 'ploidy', 'Y_call', 'purity', 'CANCER_TYPE', 'MSI_TYPE', 
                           'MSI_SCORE', 'Age', 'TMB', 'Mut_Count', 'SAMPLE_TYPE', 'WGD', 'FGA')

clinical_glm$Mut_Count = NULL
clinical_glm$TMB = NULL
# clinical_glm$MSI_TYPE = NULL
clinical_glm$MSI_SCORE = NULL
clinical_glm$WGD = NULL

##----------------+
## Merging
##----------------+
log_input_df = clinical_glm %>%
  left_join(mutation_matrix) %>%
  mutate(Age = as.numeric(Age))

# Change gene factors to numeric
log_input_df[,colnames(mutation_matrix)[1:ncol(mutation_matrix)-1]] = log_input_df[, colnames(mutation_matrix)[1:ncol(mutation_matrix)-1]] %>%
  mutate_if(is.factor, as.character) %>%
  mutate_if(is.character, as.integer)

# convert TP53 status as a fixed variable
log_input_df$TP53_Status = ifelse(log_input_df$TP53 == 1, 1, 0)
log_input_df$TP53 = NULL

##----------------+
## Adding CNA data
##----------------+
gam = load('~/Desktop/1.0.genomic_data_all.Rdata')
gam_cna = as.data.frame(data_gam_onc$cna)
gam_cna = spread(gam_cna, Var2, Freq)
colnames(gam_cna)[1] = 'sample'
row.names(gam_cna) = gam_cna$sample
gam_cna$sample = NULL
gam_cna[is.na(gam_cna)] = 0
gam_cna$sample = row.names(gam_cna)
row.names(gam_cna) = NULL


log_input_df = log_input_df %>% 
  left_join(gam_cna)


##----------------+
## run gene-wise glm
##----------------+
# Filter for gene then for cancer type
# Run logistic regression if:
## 1. gene has at least 1 mutant and at least 1 wild type
## 2. gene has alteration frequency of >= 0.03
logistic_reg_fun = function(data_frame, gene, cancer_type){
  try({
    data_frame <- data_frame[is.na(data_frame[, gene]) == F &
                               data_frame[,"CANCER_TYPE"] == cancer_type,]
    if (length(which(data_frame[,gene] == 1)) == 0 | 
        length(which(data_frame[,gene] == 0)) == 0){
      log_results_df <- data.frame(variable = "Not Tested",
                                   gene = gene, 
                                   cancer_type = cancer_type, 
                                   comments = "No Mutations in this gene") 
    } else if (length(which(data_frame[,gene] == 1))/nrow(data_frame) < 0.03) {
      log_results_df <- data.frame(variable = "Not Tested",
                                   gene = gene, 
                                   cancer_type = cancer_type, 
                                   comments = "Mutation frequency <3%") 
    } else { 
      
      formula <- as.formula(paste0(gene, "~ Y_call + SAMPLE_TYPE + FGA + purity + TP53_Status + MSI_TYPE"))
      log_results <- glm(formula, data = data_frame, family = binomial)
      log_results_df <- as.data.frame(summary(log_results)$coefficients)
      log_results_df$variable <- row.names(log_results_df)
      log_results_df$gene <- gene
      log_results_df$cancer_type <- cancer_type
      colnames(log_results_df)[1:4] <- c("estimate", "std_err", "z_value", "p_value")
      conf_df <- as.data.frame(confint.default(log_results))
      conf_df$variable <- row.names(conf_df)
      log_results_df <- left_join(log_results_df, conf_df, by = "variable")
      log_results_df <- log_results_df[c(5,7,6,1:4,8,9)]
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
colnames(log_input_df) = gsub("-", "_", colnames(log_input_df))
gene_list = gsub("-", "_", colnames(mutation_matrix)[1:ncol(mutation_matrix)-1])
gene_list = GOI
gene_list = gene_list[!gene_list %in% c('TP53', "HLA-A_Deletion", "HLA-B_Deletion", "NKX2-1_Amplification", "NKX3-1_Deletion")]
top20 = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/04_Loss/IMPACT_Y_loss_incidences_top20.txt', sep = '\t')
cancer_type_list = unique(as.character(top20$CancerType))


# Expand table for gene list and cancer list (get every combo)
gene_cancer_df = expand.grid(gene_list, cancer_type_list)
gene_cancer_df =  gene_cancer_df %>%
  mutate_if(is.factor, as.character)


# Run logistic regression
log_results_df = mapply(logistic_reg_fun,
                         gene = gene_cancer_df$Var1,
                         cancer_type = gene_cancer_df$Var2,
                         MoreArgs = list(data_frame = log_input_df),
                         SIMPLIFY = FALSE)
for(i in 1:length(log_results_df)){
  if(class(log_results_df[[i]]) != 'data.frame'){
    print(i)
  } else next
}

log_results_df = log_results_df[-c(2156, 4484)]
log_results_df = do.call("rbind.fill", log_results_df)
log_results_df$p_adj = p.adjust(log_results_df$p_value, method = "fdr")

write.table(x = log_results_df, file = 'Data/05_Association/gene_level/gene_level_full_out.txt', sep = '\t')





##----------------+
## Fishers exact test
##----------------+
top20 = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/04_Loss/IMPACT_Y_loss_incidences_top20.txt', sep = '\t')
oncoKB = read.csv('Data/signedOut/data_mutations_extended.oncokb.txt', sep = '\t')
cohort = readRDS('Data/signedOut/Cohort_07132022.rds')

cancer_type_list = unique(as.character(top20$CancerType))
cohort_samples = cohort$IMPACT_Y_classification_final$sample
onco_cohort = oncoKB[which(oncoKB$Tumor_Sample_Barcode %in% cohort_samples & 
                             oncoKB$Mutation_Status == 'SOMATIC'), ]

loy = merge(cohort$IMPACT_Y_classification_final, 
            cohort$IMPACT_clinicalAnnotation[,c('SAMPLE_ID', 'CANCER_TYPE')], 
            by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)

loy$Y_call[which(loy$classification %in% c('loss', 'relative_loss'))] = 'Y_chrom_loss'
loy$Y_call[which(loy$classification %in% c('gain', 'gain_loss', 'wt'))] = 'intact_Y_chrom'
loy$classification = NULL


##----------------+
## which genes to keep
##----------------+
muts_keep = as.data.frame(table(onco_cohort$Hugo_Symbol) / length(unique(onco_cohort$Tumor_Sample_Barcode)))
muts_keep = muts_keep[which(muts_keep$Freq > 0.03), ]
muts_keep = rbind(muts_keep, data.frame(Var1 = 'VHL', Freq = 0.02))
GOI = unique(as.character(muts_keep$Var1))
GOI = c(GOI, 'EIF1AX', 'KDM5C', 'PHF6', 'ZRSR2')

##----------------+
## populate alteration
## matrix;
##----------------+
PanCancer = data.frame()
for(i in unique(GOI)){
  loss_n = 5011
  wt_n = 9311
  gene = i
  cancer = 'PanCancer'
  noMuts = onco_cohort[which(onco_cohort$Hugo_Symbol != i), ]
  Muts = onco_cohort[which(onco_cohort$Hugo_Symbol == i), ]
  
  #' LOY with at least 1 mut
  LOY_oneMut = length(intersect(loy$sample[which(loy$Y_call == 'Y_chrom_loss')],
                                Muts$Tumor_Sample_Barcode))
  #' LOY with no Mutation
  LOY_wt = loss_n - LOY_oneMut
  
  #' WT tumor with at least 1 mutation
  WT_oneMut = length(intersect(loy$sample[which(loy$Y_call == 'intact_Y_chrom')],
                               Muts$Tumor_Sample_Barcode))
  #' WT tumor with no Mutation
  WT_noMut = wt_n - WT_oneMut
  
  out = data.frame(cohort = cancer,
                   gene = gene,
                   LOY_1Mut = LOY_oneMut,
                   LOY_wt = LOY_wt,
                   WT_oneMut = WT_oneMut,
                   WT_wt = WT_noMut)
  PanCancer = rbind(PanCancer, out)
}

write.table(PanCancer, file = 'Data/05_Association/gene_level/PanCancerFisher.txt', sep = '\t')


##----------------+
## individual cancer types
## and genes;
##----------------+
loy_main = loy[which(loy$CANCER_TYPE %in% cancer_type_list), ]

all_cancer_gene = data.frame()
for(i in unique(loy_main$CANCER_TYPE)){
  cancer = i
  loss_n = table(loy_main$Y_call[which(loy_main$CANCER_TYPE == i)])[['Y_chrom_loss']]
  wt_n = sum(table(loy_main$Y_call[which(loy_main$CANCER_TYPE == i)])) - loss_n
  
  for(j in unique(GOI)){
    gene = j
    Muts = onco_cohort[which(onco_cohort$Hugo_Symbol == j), ]
    
    #' LOY with at least 1 mut
    LOY_oneMut = length(intersect(loy_main$sample[which(loy_main$Y_call == 'Y_chrom_loss' & loy_main$CANCER_TYPE == i)],
                                  Muts$Tumor_Sample_Barcode))
    #' LOY with no Mutation
    LOY_wt = loss_n - LOY_oneMut
    
    #' WT tumor with at least 1 mutation
    WT_oneMut = length(intersect(loy_main$sample[which(loy_main$Y_call == 'intact_Y_chrom' & loy_main$CANCER_TYPE == i)],
                                 Muts$Tumor_Sample_Barcode))
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
}

write.table(all_cancer_gene, file = 'Data/05_Association/gene_level/CancerType_Gene_Fisher.txt', sep = '\t')


##----------------+
## Fisher Test on gene and 
## Cance type level
##----------------+
all_Fisher = rbind(PanCancer, all_cancer_gene)
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
all_Fisher$FDR_pass = ifelse(all_Fisher$FDR <0.1, 'plot', 'not')



##----------------+
## Plot the observations
##----------------+
cancerTypes = ggplot(all_Fisher[all_Fisher$cohort != 'PanCancer', ], aes(x = cohort, y = log10p, color = FDR_pass)) +
  geom_jitter(width = 0.2) +
  scale_color_manual(values = c('not' = 'grey55',
                                'plot' = 'red'),
                     name = '') +
  geom_text_repel(aes(label = ifelse(FDR_pass == 'plot', gene, '')), color = 'black') +
  coord_flip() +
  geom_vline(xintercept = seq(1.5, 20, 1), linetype = 'dashed', color = 'grey35', size = 0.4) +
  scale_y_continuous(position = 'right', expand = c(0.01, 0)) +
  theme_std(base_size = 14) +
  theme(panel.border = element_rect(fill = NA, size = 2, color = 'black'),
        legend.position = 'none') +
  labs(x = '', y = '-log10(p-value)')

#' PanCancer
PanCancer_plot = ggplot(all_Fisher[all_Fisher$cohort == 'PanCancer', ], 
                        aes(x = cohort, y = log10p, color = FDR_pass)) +
  geom_jitter(width = 0.2) +
  scale_color_manual(values = c('not' = 'grey55',
                                'plot' = 'red'),
                     name = '') +
  geom_text_repel(aes(label = ifelse(FDR_pass == 'plot', gene, '')), 
                  max.overlaps = 300, position = position_dodge(width = 0.2)) +
  coord_flip() +
  scale_y_continuous(position = 'right', expand = c(0.01, 0)) +
  theme_std(base_size = 14) +
  theme(panel.border = element_rect(fill = NA, size = 2, color = 'black'),
        legend.position = 'none',
        axis.text.y = element_text(angle = 90, size = 16, hjust = 0.5),
        axis.ticks.y = element_blank()) +
  labs(x = '', y = '-log10(p-value)')


geneLevel_all = PanCancer_plot/cancerTypes + plot_layout(heights = c(0.15, 1))


#' out