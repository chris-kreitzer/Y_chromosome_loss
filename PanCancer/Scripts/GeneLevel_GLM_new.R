##----------------+
## Fisher's exact test
## for mutations in genes
## associated with LOY
##----------------+

clean()
gc()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')

cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
cohort_samples = cohort$IMPACT_Y_classification_final$sample
oncoKB = read.csv('Data/signedOut/data_mutations_extended.oncokb.txt', sep = '\t')
onco_cohort = oncoKB[which(oncoKB$Tumor_Sample_Barcode %in% cohort_samples), ]
onco_cohort = onco_cohort[which(onco_cohort$Mutation_Status == 'SOMATIC'), ]

##----------------+
## which genes to keep
##----------------+
muts_keep = as.data.frame(table(onco_cohort$Hugo_Symbol) / length(unique(onco_cohort$Tumor_Sample_Barcode)))
muts_keep = muts_keep[which(muts_keep$Freq > 0.03), ]
muts_keep = rbind(muts_keep, data.frame(Var1 = 'VHL', Freq = 0.02))
GOI = unique(as.character(muts_keep$Var1))
GOI = c(GOI, 'EIF1AX', 'KDM5C', 'PHF6', 'ZRSR2')

##----------------+
## set-up mutation matrix
##----------------+
mutation_matrix = setNames(data.frame(matrix(ncol = length(unique(GOI)), nrow = 0)), unique(GOI))
for(id in 1:length(cohort_samples)){
  print(id)
  data.mut.sub = onco_cohort[which(onco_cohort$Tumor_Sample_Barcode == cohort_samples[id]), ]
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
mutation_matrix$Sample.ID = NULL
mutation_matrix = mutation_matrix[, which(colnames(mutation_matrix) %in% GOI)]
mutation_matrix[is.na(mutation_matrix)] = 0
mutation_matrix$sample = row.names(mutation_matrix)
row.names(mutation_matrix) = NULL
mutation_matrix = as.data.frame.matrix(mutation_matrix)


##----------------+
## Merge Clinical and gene
## data information
##----------------+
cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
clinical = cohort$IMPACT_clinicalAnnotation
loy = cohort$IMPACT_Y_classification_final
loy_clinical = cohort$IMPACT_binaryY_call
arm_level = cohort$IMPACT_ARM_level_changes
clinical_glm = merge(loy[,c('sample', 'ploidy', 'classification')], loy_clinical[,c('sample.id', 'purity')],
                     by.x = 'sample', by.y = 'sample.id', all.x = T)

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
clinical_glm$MSI_TYPE = NULL
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

##----------------+
## run gene-wise glm
##----------------+
# Filter for gene then for cancer type
# Run logistic regression if:
## 1. gene has at least 1 mutant and at least 1 wild type
## 2. gene has alteration frequency of >= 0.03
logistic_reg_fun = function(data_frame, gene, cancer_type){
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
    
    formula <- as.formula(paste0(gene, "~ Y_call + SAMPLE_TYPE + FGA + purity"))
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
}

##----------------+
## Set gene list and 
## cancer type list
##----------------+
colnames(log_input_df) = gsub("-", "_", colnames(log_input_df))
gene_list = gsub("-", "_", colnames(mutation_matrix)[1:ncol(mutation_matrix)-1])

top20 = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/04_Loss/IMPACT_Y_loss_incidences_top20.txt', sep = '\t')
cancer_type_list = unique(as.character(top20$CancerType))

# Expand table for gene list and cancer list (get every combo)
gene_cancer_df = expand.grid(gene_list, cancer_type_list)
gene_cancer_df <-  gene_cancer_df %>%
  mutate_if(is.factor, as.character)


# Run logistic regression
log_results_df = mapply(logistic_reg_fun,
                         gene = gene_cancer_df$Var1,
                         cancer_type = gene_cancer_df$Var2,
                         MoreArgs = list(data_frame = log_input_df),
                         SIMPLIFY = FALSE)
log_results_df <- do.call("rbind.fill", log_results_df)
log_results_df$p_adj <- p.adjust(log_results_df$p_value, method = "fdr")





View(log_results_df)


gam = load('~/Desktop/1.0.genomic_data_all.Rdata')


gam
head(data_MAF$mut)
str(data_gam_onc$gene)
aa = dimnames(data_gam_onc$gene)[1]
length(unique(aa))
aa = as.character(unlist(aa))


length(intersect(cohort_samples, aa))

gam


head(data_samples)

TeaTasting <-
  matrix(c(3, 1, 1, 3),
         nrow = 2,
         dimnames = list(Guess = c("Milk", "Tea"),
                         Truth = c("Milk", "Tea")))


fisher.test(TeaTasting, alternative = "greater")
fisher.test(matrix(c(100, 147, 1290, 2841), ncol = 2), alternative = 'greater')

