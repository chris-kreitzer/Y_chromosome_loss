##----------------+
## Paralogue chrX - chrY
## interactions; 
## biallelic VS monoallelic 
## inactivations;
## enrichment for biallelic loss might be 
## an indication of TSG of the Y chromosome
##----------------+
##
## start: 01/24/2023
## revision: 01/25/2023
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
dmp_facets_gene = dmp_facets_gene[which(dmp_facets_gene$sample %in% downsample), ]
dmp_muts = data.table::fread('Data/signedOut/data_mutations_extended.oncokb.txt.gz', sep = '\t') 
dmp_muts = dmp_muts[which(dmp_muts$Tumor_Sample_Barcode %in% downsample), ]
dmp_cna = data.table::fread('Data/signedOut/data_CNA.oncokb.txt.gz', sep = '\t')
data_gam = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = data_gam$GAM_Analysis


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

xx = male_out
xx$chrX_gene = ifelse(xx$chrX_gene_filter %in% c('suppress_segment_too_large', 'suppress_large_homdel', 'suppress_likely_unfocal_large_gain'),  'wt', xx$chrX_gene)
xx = xx[!duplicated(xx), ]
xx$chrX_gene = ifelse(xx$chrX_gene == 'HOMDEL', 'HOMDEL', 'wt')
xx$chrX = ifelse(xx$chrX_mut == xx$chrX_cna & xx$chrX_mut == xx$chrX_gene, 'wt', 'mono')

yy = merge(xx, cohort[,c('SAMPLE_ID', 'classification')], by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
yy = yy[!duplicated(yy), ]
colnames(yy)[9] = 'chrY'
yy$chrY = ifelse(yy$chrY %in% c('complete_loss', 'relative_loss', 'partial_loss', 'gain_loss'), 'loss', 'wt')


male_summary = data.frame()
for(i in unique(yy$cancer)){
  data_sub = yy[which(yy$cancer == i), ]
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


