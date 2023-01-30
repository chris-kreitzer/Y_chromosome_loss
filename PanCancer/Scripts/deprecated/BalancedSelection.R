##----------------+
## Balanced selection
## for LOY; consequence
## of CIN
##----------------+

library(tidyr)
library(dplyr)
clean()
gc()
gam = load('~/Desktop/1.0.genomic_data_all.Rdata')
cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
CancerTypes = c('Skin Cancer, Non-Melanoma', 'Glioma', 'Soft Tissue Sarcoma',
                'Prostate Cancer', 'Bone Cancer', 'Bladder Cancer')

loy_balanced = cohort$IMPACT_clinicalAnnotation[which(cohort$IMPACT_clinicalAnnotation$SAMPLE_ID %in% cohort$IMPACT_Y_classification_final$sample), ]
loy_balanced = loy_balanced[which(loy_balanced$CANCER_TYPE %in% CancerTypes), c('SAMPLE_ID', 'CANCER_TYPE')]

GenomeScores = read.csv('Data/05_Association/gene_level/GenomeScores.txt', sep = '\t')
Pathway_alterations = as.data.frame(data_gam_onc$OncoPath)
Pathway_alterations$SAMPLE_ID = row.names(Pathway_alterations)
row.names(Pathway_alterations) = NULL
GenomeScores_Selected = GenomeScores[which(GenomeScores$CancerType %in% CancerTypes), c('CancerType', 'fraction_LOY')]
balanced_pathway = left_join(loy_balanced, Pathway_alterations)


##----------------+
## populate the matrix
##----------------+
Pathway_out = data.frame()
for(i in unique(balanced_pathway$CANCER_TYPE)){
  cancer = i
  for(j in c('RTK_RAS', 'TP53', 'Epigenetic')){
    pathway = j
    n_altered = table(balanced_pathway[which(balanced_pathway$CANCER_TYPE == i), j])[[2]]
    n_all = sum(table(balanced_pathway[which(balanced_pathway$CANCER_TYPE == i), j]))
    frac = n_altered / n_all
    out = data.frame(cohort = cancer,
                     pathway = pathway,
                     fraction = frac)
    Pathway_out = rbind(Pathway_out, out)
  }
}

Pathway_out = merge(Pathway_out, GenomeScores_Selected, by.x = 'cohort', by.y = 'CancerType', all.x = T)


ggplot(Pathway_out[!Pathway_out$cohort %in% c('Bladder Cancer', 'Glioma'), ], aes(x = fraction, y = fraction_LOY, color = pathway, shape = cohort)) +
  geom_point() +
  geom_line(group = 2)



dim(Pathway_out)
TP53_out = merge(TP53_out, GenomeScores_Selected, by.x = 'cohort', by.y = 'CancerType', all.x = T)
TP53_out
cor.test(TP53_out$frac_TP53, TP53_out$fraction_LOY, method = 'spearman')