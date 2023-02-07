##----------------+
## Paralogue chrX - chrY
## interactions; 
## biallelic VS monoallelic 
## inactivations;
## enrichment for biallelic loss might be 
## an indication of TSG of the Y chromosome
## 
## Tendency for two-hit event (indicative of TSG)
## male VS female
##----------------+
##
## start: 01/24/2023
## revision: 01/25/2023
## revision: 01/26/2023
## revision: 01/27/2023
## revision: 02/06/2023
## 
## chris-kreitzer


clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/MSKCC/Scripts/Allelic_Status_by_gene.R')

cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
allFemaleSamples_df = read.csv('Data/06_Sex_Disparity/AllFemale.tsv', sep = '\t')
allFemaleSamples = unique(allFemaleSamples_df$Sample.ID)
downsample = union(unique(cohort$SAMPLE_ID), unique(allFemaleSamples))

dmp_facets_gene = data.table::fread('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/msk_impact_facets_annotated.gene_level.txt.gz')
dmp_facets_gene$sample = substr(x = dmp_facets_gene$sample, start = 1, stop = 17)
dmp_facets_gene = dmp_facets_gene[which(dmp_facets_gene$sample %in% cohort$SAMPLE_ID), ]
dmp_muts = data.table::fread('Data/signedOut/data_mutations_extended.oncokb.txt.gz', sep = '\t') 
dmp_muts = dmp_muts[which(dmp_muts$Tumor_Sample_Barcode %in% cohort$SAMPLE_ID), ]
dmp_cna = data.table::fread('Data/signedOut/data_CNA.oncokb.txt.gz', sep = '\t')
data_gam = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = data_gam$GAM_Analysis
gc()

femaleReno = read.csv('Data/06_Sex_Disparity/FemaleRenoSamples.tsv', sep = '\t')
femaleReno = unique(femaleReno$Sample.ID)


##----------------+
## EXITS genes:
## and non EXITS genes
## Loss definitions
##----------------+
EXIT = c('ATRX', 'EIF1AX', 'KDM6A', 'KDM5C', 'ZRSR2')
nonExit = c('AMER1', 'AR', 'ARAF', 'BCOR', 'BTK', 'CRLF2', 'GATA1', 'MED12', 'RBM10', 'SH2D1A', 'STAG2', 'XIAP', 'PHF6')
loss = c('complete_loss', 'gain_loss', 'relative_loss', 'partial_loss')
wt = c('wt', 'gain', 'partial_gain')



##----------------+
## Focus on KDM5C;
## - focus on Renal-Cell Carcinoma
## - focus an female first
## - decipher whether we can disentangle 
## - mono-, biallelic and wt cases
## 
## - then move on to male samples and see 
## - whether the same pattern exits for
## - individual genes first; then all EXITs genes
##----------------+
##
## Hypothesis:
## - female with EXIT gene mutations are more likely to have lost
## - the remaining X-chromosome, suggesting enrichment for biallelic loss
##
## - for Kidney Cancer there is no CNA event recorded 
## - for FEMALES on the X-chromosome
##----------------+

##----------------+
## Automatization 
##----------------+
clean()
gc()
.rs.restartR()
setup('~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/MSKCC/Scripts/Allelic_Status_by_gene.R')

allFemaleSamples_df = read.csv('Data/06_Sex_Disparity/AllFemale.tsv', sep = '\t')
dmp_facets_gene = data.table::fread('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/msk_impact_facets_annotated.gene_level.txt.gz')
dmp_facets_gene$sample = substr(x = dmp_facets_gene$sample, start = 1, stop = 17)
dmp_muts = data.table::fread('Data/signedOut/data_mutations_extended.oncokb.txt.gz', sep = '\t') 

EXIT = c('ATRX', 'EIF1AX', 'KDM6A', 'KDM5C', 'ZRSR2')
nonExit = c('AMER1', 'AR', 'ARAF', 'BCOR', 'BTK', 'CRLF2', 'GATA1', 'MED12', 'RBM10', 'SH2D1A', 'STAG2', 'XIAP', 'PHF6')
chrX_genes = union(EXIT, nonExit)

female_out = data.frame()
for(i in unique(ctypes_keep)){
  try({
    cancer = i
    print(i)
    female_samples = unique(allFemaleSamples_df[which(allFemaleSamples_df$Cancer.Type == i), 'Sample.ID'])
    female_cna = as.data.frame(dmp_facets_gene[which(dmp_facets_gene$sample %in% female_samples), ])
    female_muts = as.data.frame(dmp_muts[which(dmp_muts$Tumor_Sample_Barcode %in% female_samples), ])
    
    for(j in unique(chrX_genes)){
      gene = j
      print(j)
      allelic_status_f = allelic_status(samples = female_samples,
                                        gene = j,
                                        copy_number_data = female_cna,
                                        mutation_data = female_muts)
      allelic_status_f = allelic_status_f[!duplicated(allelic_status_f), ]
      allelic_status_f = allelic_status_f[!allelic_status_f$allelic_call %in% c('ambiguous:FacetsFilter', 'check Facets fit'), ]
      
      nowt = length(allelic_status_f$id[which(allelic_status_f$mutation == 'none' & 
                                                allelic_status_f$cna_AI_n == 0)])
      
      nomut = length(allelic_status_f$id[which(allelic_status_f$mutation == 'none' & 
                                                 allelic_status_f$cna_AI_n != 0)])
      
      yeswt = length(allelic_status_f$id[which(allelic_status_f$mutation != 'none' & 
                                                 allelic_status_f$cna_AI_n == 0)])
      
      yesmut = length(allelic_status_f$id[which(allelic_status_f$mutation != 'none' & 
                                                  allelic_status_f$cna_AI_n != 0)])
      
      ##' assess significance
      ODDS = fisher.test(matrix(c(nowt, nomut, yeswt, yesmut), ncol = 2))$estimate[[1]]
      p.value = fisher.test(matrix(c(nowt, nomut, yeswt, yesmut), ncol = 2))$p.value
      
      out = data.frame(Cancer = cancer,
                       cohort = 'female',
                       gene = gene,
                       biallelic_n = yesmut,
                       mono_mutation = yeswt,
                       mono_cna = nomut,
                       total = length(unique(allelic_status_f$id)),
                       ODDS = ODDS,
                       p.value = p.value)
      
      female_out = rbind(female_out, out)
    }
    rm(female_samples, female_cna, female_muts, allelic_status_f)
  })
}


##-------
## adjust p value and
## save the output
female_out$p_adj = p.adjust(p = female_out$p.value, method = 'fdr')
write.table(female_out, file = 'Data/06_Sex_Disparity/Female_biallelicEnrichment.txt', sep = '\t', quote = F, row.names = F)




##----------------+
## same approach for male
## tumor samples;
##----------------+
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
ctypes_keep
male_out = data.frame()
for(i in unique(ctypes_keep)){
  try({
    print(i)
    cancer = i
    male_samples = unique(cohort[which(cohort$CANCER_TYPE == i), 'SAMPLE_ID'])
    male_muts = dmp_muts[which(dmp_muts$Tumor_Sample_Barcode %in% male_samples), ]
    male_cna = dmp_cna[which(dmp_cna$SAMPLE_ID %in% male_samples), ]
    male_dmp_cna = dmp_facets_gene[which(dmp_facets_gene$sample %in% male_samples), ]
    
    
    
    for(j in unique(male_samples)){
      for(k in unique(chrX_genes)){
        print(k)
        if(j %in% male_muts$Tumor_Sample_Barcode &
           length(male_muts$Tumor_Sample_Barcode[which(male_muts$Tumor_Sample_Barcode == j & male_muts$Hugo_Symbol == k)]) != 0){
          sample = j
          chrX_mut = male_muts$ONCOGENIC[which(male_muts$Tumor_Sample_Barcode == j & male_muts$Hugo_Symbol == k)]
          gene = k
        } else {
          sample = j
          chrX_mut = 'wt'
          gene = k
        }
        
        if(j %in% male_cna$SAMPLE_ID &
           length(male_cna$SAMPLE_ID[which(male_cna$SAMPLE_ID == j & male_cna$HUGO_SYMBOL == k)]) != 0){
          sample = j
          chrX_cna = male_cna$ALTERATION[which(male_cna$SAMPLE_ID == j & male_cna$HUGO_SYMBOL == k)]
          gene = k
        } else {
          sample = j
          chrX_cna = 'wt'
          gene = k
        }
        
        if(j %in% male_dmp_cna$sample &
           length(male_dmp_cna$sample[which(male_dmp_cna$sample == j & male_dmp_cna$gene == k)]) != 0){
          sample = j
          chrX_gene = male_dmp_cna$cn_state[which(male_dmp_cna$sample == j & male_dmp_cna$gene == k)]
          chrX_gene_filter = male_dmp_cna$filter[which(male_dmp_cna$sample == j & male_dmp_cna$gene == k)]
          gene = k
        } else {
          sample = j
          chrX_gene = 'wt'
          chrX_gene_filter = 'NA'
          gene = k
        }
        
        out = data.frame(cancer = i,
                         sample = sample,
                         chrX_mut = chrX_mut,
                         chrX_cna = chrX_cna,
                         chrX_gene = chrX_gene,
                         chrX_gene_filter = chrX_gene_filter,
                         gene = gene)
        
        male_out = rbind(male_out, out)
      }
    }
    rm(male_dmp_cna, male_cna, male_muts, male_samples)
  })
}

# write.table(x = male_out, file = 'Data/05_Mutation/gene_level/Male_chrX_mutations.txt', sep = '\t', row.names = F, quote = F)

male_out$chrX_gene = ifelse(male_out$chrX_gene_filter %in% c('suppress_segment_too_large', 'suppress_large_homdel', 'suppress_likely_unfocal_large_gain'), 'wt', male_out$chrX_gene)
male_out = male_out[!duplicated(male_out), ]
male_out$chrX_gene = ifelse(male_out$chrX_gene == 'HOMDEL', 'HOMDEL', 'wt')
male_out$chrX = ifelse(male_out$chrX_mut == male_out$chrX_cna & male_out$chrX_mut == male_out$chrX_gene, 'wt', 'mono')
male_out = male_out[!duplicated(male_out), ]
male_out = male_out[,c(1,2,3,4,5,7,8)]
male_out = male_out[!duplicated(male_out), ]
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]

male_out = merge(male_out, cohort[,c('SAMPLE_ID', 'classification')], by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
colnames(male_out)[8] = 'chrY'
male_out$chrY = ifelse(male_out$chrY %in% c('complete_loss', 'relative_loss', 'partial_loss', 'gain_loss'), 'loss', 'wt')


male_summary = data.frame()
for(i in unique(male_out$cancer)){
  data_sub = male_out[which(male_out$cancer == i), ]
  for(j in unique(data_sub$gene)){
    Cancer = i
    print(Cancer)
    gene = j
    
    nowt = length(data_sub$sample[which(data_sub$gene == j & data_sub$chrX == 'wt' &
                                          data_sub$chrY == 'wt')])
    
    nomut = length(data_sub$sample[which(data_sub$gene == j & data_sub$chrX == 'wt' &
                                           data_sub$chrY != 'wt')])
    
    yeswt = length(data_sub$sample[which(data_sub$gene == j & data_sub$chrX == 'mono' &
                                           data_sub$chrY == 'wt')])
    
    yesmut = length(data_sub$sample[which(data_sub$gene == j & data_sub$chrX == 'mono' &
                                            data_sub$chrY == 'loss')])
    
    ##' assess significance
    ODDS = fisher.test(matrix(c(nowt, nomut, yeswt, yesmut), ncol = 2))$estimate[[1]]
    p.value = fisher.test(matrix(c(nowt, nomut, yeswt, yesmut), ncol = 2))$p.value
    
    out = data.frame(Cancer = Cancer,
                     cohort = 'male',
                     gene = gene,
                     biallelic_n = yesmut,
                     mono_mutation = yeswt,
                     mono_cna = nomut,
                     total = length(unique(data_sub$sample)),
                     ODDS = ODDS,
                     p.value = p.value)
    
    male_summary = rbind(male_summary, out)
  }
  rm(data_sub, nowt, nomut, yeswt, yesmut, ODDS, p.value)
}


##-------
## adjust p-value and save output
##-------
male_summary$p_adj = p.adjust(p = male_summary$p.value, method = 'fdr')
write.table(male_summary, file = 'Data/06_Sex_Disparity/Male_biallelicEnrichment.txt', sep = '\t', row.names = F, quote = F)



##----------------+
## Visualization and post
## EDA; especially for CRLF2
##----------------+
female_biallelic = read.csv('Data/06_Sex_Disparity/Female_biallelicEnrichment.txt', sep = '\t')
female_biallelic = female_biallelic[!is.infinite(female_biallelic$ODDS), ]
male_biallelic = read.csv('Data/06_Sex_Disparity/Male_biallelicEnrichment.txt', sep = '\t')
male_biallelic = male_biallelic[!is.infinite(male_biallelic$ODDS), ]

fema = rbind(female_biallelic, male_biallelic)

fema_long = data.frame()
for(i in unique(fema$Cancer)){
  try({
    dataset = fema[which(fema$Cancer == i), ]
    dataset = dataset[!dataset$gene %in% 'PHF6', ]
    data_long = spread(dataset[,c('Cancer', 'cohort', 'gene', 'ODDS')], key = cohort, value = ODDS)
    data_long = do.call(data.frame,lapply(data_long, function(x) replace(x, is.infinite(x), 0)))
    print(head(data_long))
    fema_long = rbind(fema_long, data_long)
  })
}

fema_long = fema_long[complete.cases(fema_long),]
fema_long = merge(fema_long, fema[which(fema$cohort == 'male'), c('Cancer', 'gene', 'p_adj')],
                  by.x = c('Cancer', 'gene'), by.y = c('Cancer', 'gene'), all.x = T)
fema_long$plot = ifelse(fema_long$p_adj < 0.05, 'plot', 'noplot')
write.table(fema_long, file = 'Data/06_Sex_Disparity/Female_Male_summary.txt', sep = '\t', row.names = F, quote = F)


plot_fema_all = list()
for(i in unique(fema_long$Cancer)){
  max_m = ceiling(max(fema_long$male[which(fema_long$Cancer == i)]))
  
  plot_fema = ggplot(fema_long[which(fema_long$Cancer == i), ], 
                     aes(x = male, y = female, 
                         label = gene, color = plot)) +
    geom_point() +
    geom_text_repel(size = 2.5, max.overlaps = 50) +
    scale_color_manual(values = c('plot' = 'red',
                                  'noplot' = 'darkgrey')) +
    scale_y_continuous(limits = c(0, max_m)) +
    scale_x_continuous(limits = c(0, max_m)) +
    geom_abline(slope = 1, intercept = 0, 
                linetype = 'dashed', color = 'darkgrey', size = theme_get()$line$size) +
    theme_std(base_size = 14) +
    theme(aspect.ratio = 1,
          legend.position = 'none') +
    labs(title = i, x = 'male [ODDS]', y = 'female [ODDS]')
  plot_fema_all[[i]] = plot_fema
}

selected = plot_fema_all$`Bladder Cancer` + plot_fema_all$Glioma + plot_fema_all$`Pancreatic Cancer` + plot_fema_all$`Renal Cell Carcinoma` +
  plot_layout(ncol = 4)


ggsave_golden(filename = 'Figures_original/ODDS_fema_SoftSarcoma.pdf', plot = plot_fema_all$`Soft Tissue Sarcoma`, width = 6)




##----------------+
## Deep dive Bladder Cancer
##----------------+
female = read.csv('Data/06_Sex_Disparity/Female_biallelicEnrichment.txt', sep = '\t')
female = female[which(female$Cancer == 'Bladder Cancer'), ]

















##----------------+
## check CRLF2 in Renal Cell Carcinoma
##----------------+
## CRLF2 deletions are enriched in 
## - Soft Tissue Sarcoma
## - Bladder Cancer
## - Prostate Cancer
## 
## We also see strong enrichment of biallelic
## losses (enrichment) in males compared to females
## - maybe that's indicative of a TSG?
##
## Start with Renal Cell Carcinomas:
setwd('~/Documents/MSKCC/10_MasterThesis/')
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
male_Reno = unique(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE == 'Renal Cell Carcinoma')])
dmp_muts = data.table::fread('Data/signedOut/data_mutations_extended.oncokb.txt.gz', sep = '\t') 
EXIT = c('EIF1AX', 'KDM6A', 'KDM5C', 'ZRSR2')
nonExit = c('ATRX', 'AMER1', 'AR', 'ARAF', 'BCOR', 'BTK', 'CRLF2', 'GATA1', 'MED12', 'RBM10', 'SH2D1A', 'STAG2', 'XIAP', 'PHF6')
Y_homologue = c('EIF1AX', 'KDM5C', 'KDM6A')
chrX_genes = union(EXIT, nonExit)



##' manual check of chrX alterations (n=18 genes)
#' mono mutation on any X
reno_X_muts = dmp_muts[which(dmp_muts$Tumor_Sample_Barcode %in% male_Reno), ]
reno_X_muts = reno_X_muts[which(reno_X_muts$Hugo_Symbol %in% chrX_genes), ]
reno_X_muts_samples = unique(reno_X_muts$Tumor_Sample_Barcode)


#' mono cna on any X
remaining_males = male_Reno[!male_Reno %in% reno_X_muts_samples]

dmp_facets_gene = data.table::fread('~/Documents/MSKCC/10_MasterThesis/Data/signedOut/msk_impact_facets_annotated.gene_level.txt.gz')
dmp_facets_gene$sample = substr(x = dmp_facets_gene$sample, start = 1, stop = 17)

reno_cna = dmp_facets_gene[which(dmp_facets_gene$sample %in% remaining_males), ]
reno_cna = reno_cna[which(reno_cna$gene %in% chrX_genes), ]
reno_cna = reno_cna[!reno_cna$filter %in% c('suppress_segment_too_large', 'suppress_large_homdel'), ]

reno_x_cna_samples = unique(reno_cna$sample[which(reno_cna$cn_state == 'HOMDEL')])


#' wt on any X
reno_X_wt = male_Reno[!male_Reno %in% union(reno_X_muts_samples, reno_x_cna_samples)]


#' compile
reno_X_muts_samples
chrX_genes
Y_homologue

muts_out = data.frame()
for(i in unique(reno_X_muts_samples)){
  mu = dmp_muts[which(dmp_muts$Tumor_Sample_Barcode == i), ]
  for(j in unique(chrX_genes)){
    if(j %in% mu$Hugo_Symbol){
      gene = j
      oncogenic = paste(mu$ONCOGENIC[which(mu$Hugo_Symbol == j)], collapse = '')
      n = length(mu$Hugo_Symbol[which(mu$Hugo_Symbol == j)])
    } else {
      gene = j
      oncogenic = NA
      n = NA
    }
    out = data.frame(sample = i,
                     gene = gene,
                     alteration = oncogenic,
                     n = n)
    muts_out = rbind(muts_out, out)
  }
}

all_out = data.frame()
for(i in unique(muts_out$sample)){
  exit = muts_out[which(muts_out$sample == i & muts_out$gene %in% EXIT), ]
  if(any(!is.na(exit$alteration))){
    EXIT_mut = 1
  } else {
    EXIT_mut = 0
  }
  homo = muts_out[which(muts_out$sample == i & muts_out$gene %in% Y_homologue), ]
  if(any(!is.na(homo$alteration))){
    homo_mut = 1
  } else {
    homo_mut = 0
  }
  out = rbind(muts_out[which(muts_out$sample == i), ],
              data.frame(sample = i,
                         gene = c('EXIT', 'Y_homologue'),
                         alteration = c(EXIT_mut, homo_mut),
                         n = c(EXIT_mut, homo_mut)))
  all_out = rbind(all_out, out)
              
}


##----------------+
## CNA; n = 13 samples
##----------------+
reno_cna_short = reno_cna[which(reno_cna$sample %in% reno_x_cna_samples), c('sample', 'gene', 'cn_state', 'filter')]
reno_cna_short$cn_state = ifelse(reno_cna_short$cn_state == 'HOMDEL', 'Oncogenic', NA)
colnames(reno_cna_short)[3] = 'alteration'
colnames(reno_cna_short)[4] = 'n'

all_cna = data.frame()
for(i in unique(reno_cna_short$sample)){
  exit = reno_cna_short[which(reno_cna_short$sample == i & reno_cna_short$gene %in% EXIT), ]
  if(any(!is.na(exit$alteration))){
    EXIT_mut = 1
  } else {
    EXIT_mut = 0
  }
  homo = reno_cna_short[which(reno_cna_short$sample == i & reno_cna_short$gene %in% Y_homologue), ]
  if(any(!is.na(homo$alteration))){
    homo_mut = 1
  } else {
    homo_mut = 0
  }
  out = rbind(reno_cna_short[which(reno_cna_short$sample == i), ],
              data.frame(sample = i,
                         gene = c('EXIT', 'Y_homologue'),
                         alteration = c(EXIT_mut, homo_mut),
                         n = c(EXIT_mut, homo_mut)))
  all_cna = rbind(all_cna, out)
}



##----------------+
## wt cases
##----------------+
wt_out = data.frame()
for(i in unique(reno_X_wt)){
  out = data.frame(sample = i,
                   gene = c(chrX_genes, 'EXIT', 'Y_homologue'),
                   alteration = NA,
                   n = NA)
  wt_out = rbind(wt_out, out)
}


##----------------+
## male reno cohort X
maleRenoX = rbind(all_out, all_cna, wt_out)
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')

maleRenoX = merge(maleRenoX, cohort[,c('SAMPLE_ID', 'classification')], by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
maleRenoX$classification = ifelse(maleRenoX$classification %in% c('complete_loss', 'gain_loss', 'partial_loss'), 'loss', 'wt')
colnames(maleRenoX)[5] = 'chrY'

reno_out = data.frame()
for(i in unique(maleRenoX$gene)){
  gene = i
  
  nowt = length(maleRenoX$sample[which(maleRenoX$gene == i & is.na(maleRenoX$alteration) &
                                         maleRenoX$chrY == 'wt')])
  
  nomut = length(maleRenoX$sample[which(maleRenoX$gene == i & is.na(maleRenoX$alteration) &
                                          maleRenoX$chrY == 'loss')])
  
  yeswt = length(maleRenoX$sample[which(maleRenoX$gene == i & !is.na(maleRenoX$alteration) &
                                          maleRenoX$chrY == 'wt')])
  
  yesmut = length(maleRenoX$sample[which(maleRenoX$gene == i & !is.na(maleRenoX$alteration) &
                                           maleRenoX$chrY == 'loss')])
  
  ##' assess significance
  ODDS = fisher.test(matrix(c(nowt, nomut, yeswt, yesmut), ncol = 2))$estimate[[1]]
  p.value = fisher.test(matrix(c(nowt, nomut, yeswt, yesmut), ncol = 2))$p.value
  
  out = data.frame(gene = gene,
                   nowt = nowt,
                   nomut = nomut,
                   yeswt = yeswt,
                   yesmut = yesmut,
                   ODDS = ODDS,
                   p.value = p.value)
  reno_out = rbind(reno_out, out)
}

write.table(reno_out, file = 'Data/06_Sex_Disparity/Male_RenalCellCarcinoma_Fisher.txt', sep = '\t', row.names = F, quote = F)












##----------------+
## Investigate gene-level GLM;
## specifically: CRLF2 and Prostate;
##----------------+
clean()
gc()
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
IGV = read.csv('Data/04_Loss/010523/IGV_out.seg', sep = '\t')
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

##-------
## Prostate cancer with LOY
##-------
Prostate_LOY = clinical_glm$sample[which(clinical_glm$CANCER_TYPE == 'Prostate Cancer' & clinical_glm$Y_call == 1)]
Prostate_nonLOY = clinical_glm$sample[which(clinical_glm$CANCER_TYPE == 'Prostate Cancer' & clinical_glm$Y_call == 0)]
Prostate_LOY_IGV = IGV[which(IGV$ID %in% Prostate_LOY), ]
Prostate_LOY_IGV$group = 'LOY'
Prostate_nonLOY_IGV = IGV[which(IGV$ID %in% Prostate_nonLOY), ]
Prostate_nonLOY_IGV$group = 'nonLOY'

prostate = rbind(Prostate_LOY_IGV, Prostate_nonLOY_IGV)

write.table(Prostate_LOY_IGV, file = '~/Desktop/Prostate_IGV.seg', sep = '\t', quote = F, row.names = F)
write.table(prostate[,c(1:6)], file = '~/Desktop/ProstateA_IGV.seg', sep = '\t', quote = F, row.names = F)

write.table(x = prostate[,c('ID', 'group')], file = '~/Desktop/atributes.txt', sep = '\t', row.names = F, quote = F)


dim(Prostate_LOY_IGV)






#' out