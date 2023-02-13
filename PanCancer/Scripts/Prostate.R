##----------------+
## Prostate; LOY investigations
##----------------+


clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')


cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
prostate = cohort[which(cohort$CANCER_TYPE == "Prostate Cancer"), ]
prostate = prostate[!is.na(prostate$classification), ]
data_gam = readRDS('Data/05_Mutation/data_gam.rds')
data_gam = data_gam$GAM_Analysis
prad_muts = data_gam[which(data_gam$sample %in% prostate$SAMPLE_ID), ]

full_prostate = merge(prostate, prad_muts, by.x = 'SAMPLE_ID', by.y = 'sample', all.x = T)

full_prostate$TP53_Deletion

table(prostate$classification)
sum(table(prostate$classification))
dim(prostate)
which(is.na(prostate$classification))
View(prostate)

cbio(prostate$SAMPLE_ID[which(prostate$classification == 'complete_loss')])

pro_out = data.frame()
for(i in c('complete_loss', 'wt')){
  for(j in 42:ncol(full_prostate)){
    print(j)
    alteration = colnames(full_prostate)[j]
    freq = sum(full_prostate[which(full_prostate$classification == i), j], na.rm = T)
    rel_freq = freq / length(full_prostate$SAMPLE_ID[which(full_prostate$classification == i)])
    category = i
    out = data.frame(category = category,
                     alteration = alteration,
                     frequency = rel_freq)
    pro_out = rbind(pro_out, out)
  }
}

pro_summary = data.frame()
for(i in unique(pro_out$alteration)){
  complete_loss = pro_out$frequency[which(pro_out$alteration == i & pro_out$category == 'complete_loss')]
  wt = pro_out$frequency[which(pro_out$alteration == i & pro_out$category == 'wt')]
  alteration = i
  out = data.frame(alteration = alteration,
                   wt = wt,
                   complete_loss = complete_loss)
  pro_summary = rbind(pro_summary, out)
}
pro_summary$display = ifelse(abs(pro_summary$wt - pro_summary$complete_loss) > 0.05, 'plot', 'notplot')

alteration_frequency = ggplot(pro_summary, aes(x = wt, y = complete_loss, color = display, label = alteration)) +
  geom_text_repel(aes(label = ifelse(pro_summary$display == 'plot', pro_summary$alteration, ''))) +
  geom_jitter() +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.45)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 0.45)) +
  scale_color_manual(values = c('plot' = 'darkred',
                                'notplot' = 'grey65',
                                name = '')) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = 'dashed', color = 'grey35') +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1,
        legend.position = 'none') +
  labs(x = 'wild type', y = 'complete loss')

ggsave_golden(filename = 'Figures_original/Prostate_AlterationFreq.pdf', plot = alteration_frequency, width = 6)





