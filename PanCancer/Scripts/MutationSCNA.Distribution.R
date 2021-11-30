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


#' Visualization of Esophagastric cancer;
ggplot(ECC_full_merged, aes(x = cnas_oncokb, y = mutations)) +
  stat_density_2d( 
                  geom = "raster", contour = T) +
  scale_fill_distiller(palette = 4, direction = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
  theme(
    legend.position='none'
  )
  

stat_density_2d(geom = "polygon", aes(alpha = ..level..))



ex1 = CNA_OncoKB[grepl(pattern = 'P-0042825-T01-IM6.*', x = CNA_OncoKB$SAMPLE_ID), ]
ex2 = Facets_Annotated[grepl(pattern = 'P-0042825-T01-IM6.*', x = Facets_Annotated$sample_id), ]
ex3 = ArmLevel_Facets[grepl(pattern = 'P-0042825-T01-IM6.*', x = ArmLevel_Facets$sample), ]
ex4 = Gene_Facets[grepl(pattern = 'P-0042825-T01-IM6.*', x = Gene_Facets$sample), ]
ex5 = Mutations[grepl(pattern = 'P-0042825-T01-IM6.*', x = Mutations$Tumor_Sample_Barcode), ]
View(ex5)


head(IMPACT468)
head(ex4)

setdiff(IMPACT468, ex4$gene)
setdiff(ex4$gene, IMPACT468)
View(IMPACT468)
