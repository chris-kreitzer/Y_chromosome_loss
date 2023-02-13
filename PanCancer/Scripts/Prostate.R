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
copynumber = read.csv('Data/04_Loss/010523/CopyNumberStates.txt', sep = '\t')
copynumber = copynumber[which(copynumber$id %in% prostate$SAMPLE_ID), ]
IGV = read.csv('Data/04_Loss/010523/IGV_out.seg', sep = '\t')

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


##----------------+
## Univariate and multivariate Cox model
##----------------+
gain_loss = full_prostate$SAMPLE_ID[which(full_prostate$classification == 'gain_loss')]
gain_loss_prad = copynumber[which(copynumber$id %in% gain_loss & copynumber$chrom == 24), ]
gain_loss_prad = gain_loss_prad[which(gain_loss_prad$tcn.em == 0), ]


pGL = c('P-0003040-T01-IM5',
'P-0013100-T01-IM5',
'P-0016288-T01-IM6',
'P-0025571-T01-IM6',
'P-0026808-T01-IM6',
'P-0029183-T01-IM6',
'P-0037069-T01-IM6',
'P-0042930-T01-IM6',
'P-0053240-T01-IM6',
'P-0060161-T01-IM7',
'P-0067186-T01-IM7',
'P-0075758-T01-IM7')

partial_loss = full_prostate$SAMPLE_ID[which(full_prostate$classification == 'partial_loss')]
partial_loss = copynumber[which(copynumber$id %in% partial_loss & copynumber$chrom == 24 & copynumber$tcn.em == 0), ]
partial_loss[,c('id', "start", 'end')]


pLoss = c('P-0001698-T01-IM3',
'P-0016450-T01-IM6',
'P-0017964-T01-IM6',
'P-0018217-T01-IM6',
'P-0018798-T01-IM6',
'P-0021792-T01-IM6',
'P-0028784-T01-IM6',
'P-0029339-T01-IM6',
'P-0033767-T01-IM6',
'P-0053666-T01-IM6',
'P-0058347-T01-IM6',
'P-0063708-T01-IM7',
'P-0069465-T01-IM7',
'P-0074727-T01-IM7',
'P-0077871-T01-IM7')

pLoss_density = partial_loss[which(partial_loss$id %in% pLoss), ]
pLoss_density = pLoss_density[-1,]
pLoss_vector = data.frame()
for(i in 1:nrow(pLoss_density)){
  print(i)
  sample = pLoss_density$id[i]
  value = seq(from = pLoss_density$start[i], to = pLoss_density$end[i], by = 1000)
  out = data.frame(sample = sample,
                   vec = value)
  pLoss_vector = rbind(pLoss_vector, out)
}

dim(pLoss_vector)
head(pLoss_vector)
head(pLoss_density)
plot(density(pLoss_vector$vec))



IGV_pLoss = IGV[which(IGV$ID %in% c(pGL) & IGV$chrom == 24), ]
breakpoints = c()
for(i in unique(IGV_pLoss$ID)){
  n = length(IGV_pLoss$seg.mean[which(IGV_pLoss$ID == i)])
  m = n-1
  if(m == 0) next
  loc.end = IGV_pLoss$loc.end[which(IGV_pLoss$ID == i)][seq(1,m,1)]
  breakpoints = c(breakpoints, loc.end)
}

plot(density(breakpoints), 
     ylab = '',
     xlab = '',
     xaxt = 'n',
     yaxt = 'n',
     main = '', lwd = 2)
axis(side = 1, at = c(1e+07, 1.5e+07, 2e+07), labels = c('10Mb', '15Mb', '20Mb'))
mtext(text = 'Genomic Position', side = 1, line = 1.8)
mtext(text = 'Density', side = 2, line = 0.8)
mtext(text = 'Gain/Loss; Chromosome Y', side = 3, line = 0.8, adj = 0, cex = 1.5)
text(x = 1.98e+07, y = 0.2e-06, label = 'KDM5D')
box(lwd = 2)



##----------------+
## univariate model
##----------------+
LOY = c(full_prostate$SAMPLE_ID[which(full_prostate$classification == 'complete_loss')],
        pLoss, pGL)

full_prostate$Y_call = ifelse(full_prostate$SAMPLE_ID %in% LOY, 'loss', 'intact')
os_prostate = full_prostate[,c(colnames(full_prostate)[1:32], 'TP53_mut', 'TP53_Deletion')]
osp = os_prostate
tmperg = read.csv('Data/05_Mutation/TMPRSS2_ETS.tsv', sep = '\t')
tmperg = as.character(tmperg$Sample.ID)
osp$SAMPLE_TYPE[osp$SAMPLE_TYPE == 'Local Recurrence'] = 'Primary'
osp$TMPERG = NA
for(i in 1:nrow(osp)){
  if(osp$SAMPLE_ID[i] %in% tmperg){
    osp$TMPERG[i] = 1
  } else {
    osp$TMPERG[i] = 0
  }
}
osp$TMPERG = as.integer(as.numeric(osp$TMPERG))
osp$TP53_Deletion = as.numeric(as.integer(osp$TP53_Deletion))

osp$TP53_Status = ifelse(osp$TP53_Deletion == 1, 1,
                               ifelse(osp$TP53_mut == 1, 1, 0))
osp$TP53_Deletion = NULL
osp$TP53_mut = NULL
osp$TP53_Status = as.integer(as.numeric(osp$TP53_Status))

summary(coxph(Surv(bcr.time, bcr.bin) ~ ZNRF3, data = znrf3.mvcox.table))













