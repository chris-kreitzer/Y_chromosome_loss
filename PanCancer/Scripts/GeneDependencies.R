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
## Automatisation 
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
ctypes_keep





female_cna = as.data.frame(dmp_facets_gene[which(dmp_facets_gene$sample %in% femaleReno), ])
female_muts = as.data.frame(dmp_muts[which(dmp_muts$Tumor_Sample_Barcode %in% femaleReno), ])

dmp_facets_gene = as.data.frame(dmp_facets_gene)
dmp_muts = as.data.frame(dmp_muts)

xx = allelic_status(samples = femaleReno, 
                    gene = 'EIF1AX', 
                    copy_number_data = female_cna,
                    mutation_data = female_muts,
                    print_progress = NULL)

AI_female_Reno = xx

AI_female_Reno = AI_female_Reno[!duplicated(AI_female_Reno), ]
AI_female_Reno = AI_female_Reno[!AI_female_Reno$allelic_call %in% c('ambiguous:FacetsFilter', 'check Facets fit'), ]

nowt = length(allelic_status_f$id[which(allelic_status_f$mutation == 'none' & 
                                        allelic_status_f$cna_AI_n == 0)])

nmut = length(allelic_status_f$id[which(allelic_status_f$mutation == 'none' & 
                                        allelic_status_f$cna_AI_n != 0)])

ywt = length(allelic_status_f$id[which(allelic_status_f$mutation != 'none' & 
                                       allelic_status_f$cna_AI_n == 0)])

ymut = length(allelic_status_f$id[which(allelic_status_f$mutation != 'none' & 
                                        allelic_status_f$cna_AI_n != 0)])
a = fisher.test(matrix(c(nowt,nmut,ywt,ymut), ncol = 2))

a$p.value
a$estimate[[1]]

##----------------+
## male Renal Cell Carcinoma
##----------------+
maleReno = unique(cohort$SAMPLE_ID[which(cohort$CANCER_TYPE == 'Renal Cell Carcinoma')])
maleReno_cna = dmp_cna[which(dmp_cna$SAMPLE_ID %in% maleReno), ]


##' male mutations
maleReno_muts = dmp_muts[which(dmp_muts$Tumor_Sample_Barcode %in% maleReno & dmp_muts$Hugo_Symbol == 'KDM5C'), c('Hugo_Symbol', 'ONCOGENIC', 'Tumor_Sample_Barcode')]
maleReno_cna = dmp_cna[which(dmp_cna$SAMPLE_ID %in% maleReno & dmp_cna$HUGO_SYMBOL == 'KDM5C'), c('SAMPLE_ID', 'HUGO_SYMBOL', 'ALTERATION')]

maleReno_gene = dmp_facets_gene[which(dmp_facets_gene$sample %in% maleReno & dmp_facets_gene$gene == 'KDM5C'), c('sample', 'gene', 'cn_state', 'filter')]
maleReno_gene = maleReno_gene[which(maleReno_gene$filter %in% c('PASS', 'RESCUE')), ]
head(maleReno_gene)

alll = data.frame()
for(i in unique(maleReno)){
  print(i)
  if(i %in% dmp_muts$Tumor_Sample_Barcode & 
     length(dmp_muts$Tumor_Sample_Barcode[which(dmp_muts$Tumor_Sample_Barcode == i & dmp_muts$Hugo_Symbol == 'KDM5C')]) != 0){
    sample = i
    chrX_mut = dmp_muts$ONCOGENIC[which(dmp_muts$Tumor_Sample_Barcode == i & dmp_muts$Hugo_Symbol == 'KDM5C')]
    gene = 'KDM5C'
  } else {
    sample = i
    chrX_mut = 'wt'
    gene = 'KDM5C'
  }
  
  if(i %in% dmp_cna$SAMPLE_ID &
     length(dmp_cna$SAMPLE_ID[which(dmp_cna$SAMPLE_ID == i & dmp_cna$HUGO_SYMBOL == 'KDM5C')]) != 0){
    sample = i
    chrX_cna = dmp_cna$ALTERATION[which(dmp_cna$SAMPLE_ID == i & dmp_cna$HUGO_SYMBOL == 'KDM5C')]
    gene = 'KDM5C'
  } else {
    sample = i
    chrX_cna = 'wt'
    gene = 'KDM5C'
  }
  
  if(i %in% dmp_facets_gene$sample &
     length(dmp_facets_gene$sample[which(dmp_facets_gene$sample == i & dmp_facets_gene$gene == 'KDM5C')]) != 0){
    sample = i
    chrX_gene = dmp_facets_gene$cn_state[which(dmp_facets_gene$sample == i & dmp_facets_gene$gene == 'KDM5C')]
    chrX_gene_filter = dmp_facets_gene$filter[which(dmp_facets_gene$sample == i & dmp_facets_gene$gene == 'KDM5C')]
    gene = 'KDM5C'
  } else {
    sample = i
    chrX_gene = 'wt'
    chrX_gene_filter = 'NA'
    gene = 'KDM5C'
  }
  out = data.frame(sample = sample,
                   chrX_mut = chrX_mut,
                   chrX_cna = chrX_cna,
                   chrX_gene = chrX_gene,
                   chrX_gene_filter = chrX_gene_filter,
                   gene = gene)
  alll = rbind(alll, out)
}

alll$chrX_gene = ifelse(alll$chrX_gene_filter == 'suppress_segment_too_large', 'wt', alll$chrX_gene)
alll$chrX_gene = ifelse(alll$chrX_gene != 'HOMDEL', 'wt', alll$chrX_gene)
alll$chrX = ifelse(alll$chrX_mut == alll$chrX_cna & alll$chrX_mut == alll$chrX_gene, 'wt', 'mono')

yy = merge(alll, cohort[,c('SAMPLE_ID', 'classification')],
           by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
                   
View(yy)

xx = female_cna[which(female_cna$sample %in% EXIT_female_reno), ]



n_EXIT_female_reno = length(unique(female_muts$Tumor_Sample_Barcode[which(female_muts$Hugo_Symbol %in% EXIT)]))
EXIT_female_reno = unique(female_muts$Tumor_Sample_Barcode[which(female_muts$Hugo_Symbol %in% EXIT)])
n_nonEXIT_female_reno = length(unique(female_muts$Tumor_Sample_Barcode[which(female_muts$Hugo_Symbol %in% nonExit)]))
nonEXIT_female_reno = unique(female_muts$Tumor_Sample_Barcode[which(female_muts$Hugo_Symbol %in% nonExit)])

length(intersect(femaleReno, unique(female_cna$sample)))



fcna = dmp_cna[which(dmp_cna$SAMPLE_ID %in% femaleReno), ]
View(fcna)


female_CNA = dmp_facets_gene[which(dmp_facets_gene$sample %in% femaleReno), ]
female_CNA = female_CNA[which(female_CNA$gene == 'KDM5C'), ]
female_CNA = female_CNA[,c('sample', 'gene', 'seg_start', 'seg_end', 'tcn.em', 'lcn.em', 'mcn', 'cn_state', 'filter')]
female_muts = female_muts[which(female_muts$Hugo_Symbol == 'KDM5C'), ]
rm(dmp_facets_gene, dmp_muts)
gc()


##-------
## sample intersect:
##-------
female_samples = union(unique(female_CNA$sample), unique(female_muts$Tumor_Sample_Barcode))
female_allelic_status_KDM5C = allelic_status(samples = female_samples, 
                                             gene = 'KDM5C',
                                             copy_number_data = female_CNA, 
                                             mutation_data = female_muts)

female_allelic_status_KDM5C = female_allelic_status_KDM5C[which(female_allelic_status_KDM5C$allelic_call %in% c('wildtype', 'monoallelic', 'biallelic')), ]

ftable = as.data.frame(table(female_allelic_status_KDM5C$allelic_call))

ggplot(ftable, aes(x = 1, y = Freq, fill = Var1)) +
  geom_bar(position = 'fill', stat = 'identity')

  