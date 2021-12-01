## Lineage specific aspects; 
## Is there any specific SCNA associated with Y-chromosome loss
## 
## 11/27/21
## chris-kreitzer


setup(working.path = '~/Documents/GitHub/Y_chromosome_loss/PanCancer/')
clean()


## Data
CNA_OncoKB = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/data_CNA.oncokb.txt', sep = '\t')
Facets_Annotated = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
Gene_Facets = vroom::vroom('~/Documents/MSKCC/05_IMPACT40K/Data/msk_impact_facets_annotated.gene_level.txt')
ArmLevel_Facets = vroom::vroom('~/Documents/MSKCC/05_IMPACT40K/Data/msk_impact_facets_annotated.arm_level.txt.gz')
IMPACT468 = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/impact468_gene_panel.txt', sep = '\t', skip = 3)
Mutations = vroom::vroom('~/Documents/MSKCC/05_IMPACT40K/Data/data_mutations_extended.somatic.oncokb.vep.maf')
cohort = readRDS('Data_out/cohort_data.rds')


#' gene panel
IMPACT468 = as.character(colnames(IMPACT468))
IMPACT468 = as.character(IMPACT468[-1])
IMPACT468 = sub(pattern = '\\.', replacement = '-', x = IMPACT468)
IMPACT468 = sub(pattern = 'MLL4', replacement = 'KMT2D', x = IMPACT468)
IMPACT468 = sub(pattern = 'MLL3', replacement = 'KMT2C', x = IMPACT468)
IMPACT468 = sub(pattern = 'MLL2', replacement = 'KMT2B', x = IMPACT468)
IMPACT468 = sub(pattern = 'MLL', replacement = 'KMT2A', x = IMPACT468)


#' Firstly, we test whether certain cancer types tend to be dominated by either
#' copy number aberration or mutations (analogues to the publication added to the google drive document)
ECC = cohort$IMPACT.cohort$SAMPLE_ID[which(cohort$IMPACT.cohort$CANCER_TYPE == 'Esophagogastric Cancer')]

mutation_assignment = function(x){
  sample = x
  print(sample)
  cnas_oncokb = nrow(CNA_OncoKB[grepl(pattern = sample, x = CNA_OncoKB$SAMPLE_ID), ])
  cnas_all = Gene_Facets[grepl(pattern = sample, x = Gene_Facets$sample), ]
  cnas_all = nrow(cnas_all[!is.na(cnas_all$cf.em & cnas_all$lcn.em), ])
  mutations = nrow(Mutations[grepl(pattern = sample, x = Mutations$Tumor_Sample_Barcode), ])
  out = data.frame(sample = sample,
                   mutations = mutations,
                   cnas_oncokb = cnas_oncokb,
                   cnas_all = cnas_all)
  out
}

#' run lapply on all ECC samples (n=541)
ECC_full = lapply(unique(ECC), function(x) mutation_assignment(x))
ECC_full_merged = data.table::rbindlist(ECC_full)
write.table(ECC_full_merged, '../Data/Mut_CNA_tmp.txt', sep = '\t', row.names = F)

#' Visualization of Esophagastric cancer;
m = ggplot(ECC_full_merged, aes(x = cnas_oncokb, y = mutations)) +
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE
  ) + scale_fill_viridis_c() +
  geom_jitter(size = 0.2) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 20)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 20))
  


## check for the overall selection for Y-chromosome loss
Impact_cohort = unique(cohort$IMPACT.cohort$SAMPLE_ID)

cohort_out = data.frame()
for(i in unique(Impact_cohort)){
  da = Facets_Annotated[grep(pattern = i, x = Facets_Annotated$sample_id), c('sample_id', 'wgd', 'fga')]
  out = data.frame(sample = i,
                   sample_long = da$sample_id,
                   wgd = da$wgd,
                   fga = da$fga)
  cohort_out = rbind(cohort_out, out)
}


#' merge the cancer types:
Cancers = merge(cohort$IMPACT.cohort, cohort_out, by.x = 'SAMPLE_ID', by.y = 'sample', all.x = T)

allCancers = data.frame()
for(i in unique(Cancers$CANCER_TYPE)){
  try({
    print(i)
    da = Cancers[which(Cancers$CANCER_TYPE == i), ]
    mu = mean(da$fga)
    sd = sd(da$fga)
    n = nrow(da)
    prop = table(da$Y_call)[[2]] / sum(table(da$Y_call))
    out = data.frame(CancerType = i,
                     meanFGA = mu,
                     sd = sd,
                     n = n,
                     propLoss = prop)
    allCancers = rbind(allCancers, out)
  })
  
}

#' only work with Cancer where its individual contribution is bigger than 1%
allCancers = allCancers[which(allCancers$n > sum(allCancers$n) * 0.01), ]

#' Visualization 
ggplot(allCancers, aes(x = propLoss, y = meanFGA, label = CancerType)) +
  geom_jitter() +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  geom_abline(intercept = 0, slope = 1) +
  #geom_pointrange(aes(ymin = meanFGA - sd, ymax = meanFGA + sd), size = 0.2) +
  geom_text(size = 2, vjust = -1, hjust = 0) +
  theme(aspect.ratio = 1)


#' specifically look into certain cancer types
# Germ Cell
# Small Cell Lung Cancer
# Renal Cell Carcinoma
# Glioma
# Prostate Cancer
# Gastrointestinal STromal Tumor
# Mesothelioma
cancer.specific = c('Germ Cell Tumor', '^Small Cell Lung Cancer', 'Renal Cell Carcinoma',
                    'Glioma', 'Prostate Cancer', 'Gastrointestinal Stromal Tumor', 'Mesothelioma')

CancerSpecific = data.frame()

for(i in unique(cancer.specific)){
  print(i)
  try({
    da = Cancers[grep(pattern = i, Cancers$CANCER_TYPE), ]
    mu.prim = mean(da$fga[which(da$SAMPLE_TYPE == 'Primary')])
    sd.prim = sd(da$fga[which(da$SAMPLE_TYPE == 'Primary')])
    n.prim = nrow(da[which(da$SAMPLE_TYPE == 'Primary'), ])
    prop.prim = table(da$Y_call[which(da$SAMPLE_TYPE == 'Primary')])[[2]] / sum(table(da$Y_call[which(da$SAMPLE_TYPE == 'Primary')]))
    
    #' look into Metastasis
    mu.met = mean(da$fga[which(da$SAMPLE_TYPE == 'Metastasis')])
    n.met = nrow(da[which(da$SAMPLE_TYPE == 'Metastasis'), ])
    sd.met = sd(da$fga[which(da$SAMPLE_TYPE == 'Metastasis')])
    prop.met = table(da$Y_call[which(da$SAMPLE_TYPE == 'Metastasis')])[[2]] / sum(table(da$Y_call[which(da$SAMPLE_TYPE == 'Metastasis')]))
    
    out = data.frame(CancerType = rep(i, 2),
                     site = c('Primary', 'Metastasis'),
                     fga = c(mu.prim, mu.met),
                     sd = c(sd.prim, sd.met),
                     n = c(n.prim, n.met),
                     prop = c(prop.prim, prop.met))
    CancerSpecific = rbind(CancerSpecific, out)
  })
  
}

CancerSpecific

#' Visualization 
ggplot(CancerSpecific, aes(x = prop, y = fga, label = CancerType)) +
  geom_point(aes(color = CancerType, shape = site)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  geom_abline(intercept = 0, slope = 1) +
  #geom_pointrange(aes(ymin = meanFGA - sd, ymax = meanFGA + sd), size = 0.2) +
  geom_text(size = 2, vjust = -1, hjust = 0) +
  theme(aspect.ratio = 1)
  



plot =  ggplot(data, aes(x=ntrunc, y=beta_best, group=INDEX, colour=INDEX)) +
  geom_point(aes(shape=detectable), na.rm=TRUE, position="dodge") +
  geom_errorbar(aes(x=ntrunc, ymax=beta_high, ymin=beta_low), na.rm=TRUE, position="dodge")










