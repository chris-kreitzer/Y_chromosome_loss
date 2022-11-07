data = readRDS('Data/signedOut/Cohort_07132022.rds')

loy = data$IMPACT_Y_classification_final
clinics = data$IMPACT_clinicalAnnotation

eso = clinics[which(clinics$CANCER_TYPE == 'Esophagogastric Cancer'), 'SAMPLE_ID']
pan = clinics[which(clinics$CANCER_TYPE == 'Pancreatic Cancer'), 'SAMPLE_ID']
pro = clinics[which(clinics$CANCER_TYPE == 'Prostate Cancer'), 'SAMPLE_ID']
gli = clinics[which(clinics$CANCER_TYPE == 'Glioma'), 'SAMPLE_ID']
table(loy$classification[which(loy$sample %in% gli)])
wgd = read.csv('Data/04_Loss/QC_metrics.txt', sep = '\t')
wgd = wgd[which(wgd$ID %in% loy$sample), ]
head(wgd)


wgd_loy = merge(loy, wgd, by.x = 'sample', by.y = 'ID', all.x = T)


head(wgd_loy)
sum(table(wgd_loy$classification[which(wgd_loy$wgd == T)]))

wgd_plot = data.frame(category = rep(c('gain', 'gain_loss', 'loss', 'relative_loss', 'wt'), 2),
                      wgd = c(rep('NONE', 5), rep(1, 5)),
                      value = c(1693, 43, 2848, 17, 5055, 1281, 75, 1967, 178, 1161))
wgd_plot$wgd = factor(wgd_plot$wgd, levels = c('NONE', '1'))
wgd_plot$new_value[1:5] = (wgd_plot$value[1:5] / 9656) * 100
wgd_plot$new_value[6:10] = (wgd_plot$value[6:10] / 4662) * 100

source('Scripts/plot_theme.R')
WGD_LOY = ggplot(wgd_plot, aes(x = wgd, y = new_value, fill = category)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'loss' = '#0E3F7C',
                               'relative_loss' = '#00AEC8',
                               'gain' = '#D53833',
                               'gain_loss' = '#E3CC98'), name = '') +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_x_discrete(labels = c('none', '1x')) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(aspect.ratio = 2) +
  labs(x = '# genome doublings', y = 'Fraction of tumors')
  
ggsave(filename = 'Figures_original/Fraction_Y_loss_WGD.pdf', plot = WGD_LOY, device = 'pdf', width = 4)
  

##----------------+
## chromosome arms lost
##----------------+
arms = read.csv('Data/04_Loss/IMPACT_arm_change_out.txt', sep = '\t')
arms = arms[which(arms$ID %in% loy$sample), c('fcna', 'ID')]
arms = arms[!duplicated(t(apply(arms, 1, sort))), ]


arms_loss = arms[which(arms$ID %in% loy$sample[which(loy$classification == 'loss')]), ]
arms_relative_loss = arms[which(arms$ID %in% loy$sample[which(loy$classification == 'relative_loss')]), ]
arms_gain = arms[which(arms$ID %in% loy$sample[which(loy$classification == 'gain')]), ]
arms_wt = arms[which(arms$ID %in% loy$sample[which(loy$classification == 'wt')]), ]
arms_gain_loss = arms[which(arms$ID %in% loy$sample[which(loy$classification == 'gain_loss')]), ]

rm(FGA)
FGA = data.frame(category = c(rep('loss', length(arms_loss$fcna)),
                              rep('relative_loss', length(arms_relative_loss$fcna)),
                              rep('gain', length(arms_gain$fcna)),
                              rep('wt', length(arms_wt$fcna)),
                              rep('gain_loss', length(arms_gain_loss$fcna))),
                 value = c(arms_loss$fcna, arms_relative_loss$fcna, arms_gain$fcna, arms_wt$fcna, arms_gain_loss$fcna))

FGA$category = factor(FGA$category, levels = c('wt', 'loss', 'relative_loss', 'gain_loss', 'gain'))


FGA_Y_loss = ggplot(FGA, aes(x = category, y = value, fill = category)) +
  geom_boxplot() +
  scale_fill_manual(values = c('wt' = '#D7D8DA',
                               'loss' = '#0E3F7C',
                               'relative_loss' = '#00AEC8',
                               'gain' = '#D53833',
                               'gain_loss' = '#E3CC98'), name = '') +
  scale_y_continuous(expand = c(0.01,0)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(aspect.ratio = 0.5) +
  labs(x = '', y = 'FGA')

ggsave(filename = 'Figures_original/FGA_Y_loss.pdf', plot = FGA_Y_loss, device = 'pdf', width = 7)

t.test(FGA$value[which(FGA$category == 'loss')], 
            FGA$value[which(FGA$category == 'gain')])


##----------------+
## Focality of loss
##----------------+
copynumber = read.csv('Data/04_Loss/IMPACT_copynumber_out.txt', sep = '\t')
copynumber = copynumber[which(copynumber$ID %in% loy$sample), ]

copy24 = copynumber[which(copynumber$chrom == 24), ]

focality = data.frame()
for(i in unique(copy24$ID)){
  n = length(copy24$chrom[which(copy24$ID == i)])
  out = data.frame(id = i,
                   n = n)
  focality = rbind(focality, out)
}

barplot(height = c(14029, 269, 24), width = c(0.2, 0.2, 0.2), 
        space = c(0.1, 0.1, 0.1), names.arg = c('1', '2', '3'), col = 'white', 
        xlab = '# segments called', yaxt = 'n')
axis(side = 2, at = seq(0, 14200, 2000), las = 2, cex = 1.2, line = -0.5, lwd = 1.5)


##----------------+
## Breakpoints
##----------------+
foc2 = focality$id[which(focality$n > 1)]
copymore = copynumber[which(copynumber$ID %in% foc2), ]
copymore = copymore[which(copymore$chrom == 24), ]
copymore = copymore[,c('start', 'end', 'ID', 'ploidy')]

all_breakpoints = data.frame()
for(i in unique(copymore$ID)){
  if(length(copymore$start[which(copymore$ID == i)]) < 3){
    breakpoint = copymore$end[which(copymore$ID == i)][1]
  } else {
    breakpoint = c(copymore$end[which(copymore$ID == i)][1],
                   copymore$end[which(copymore$ID == i)][2])
  }
  out = data.frame(id = i, 
                   breakpoint = breakpoint)
  all_breakpoints = rbind(all_breakpoints, out)
}

breakpoint_plot = ggplot(all_breakpoints) +
  geom_density(aes(x = breakpoint)) +
  scale_x_continuous(limits = c(1400000, 25000000),
                     expand = c(0,0),
                     breaks = c(1460000, 7000000, 13000000, 19000000, 25000000),
                     labels = c(1460000, 7000000, 13000000, 19000000, 25000000) / 1e6) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_std(base_size = 14, base_line_size = 0.001) +
  theme(panel.border = element_rect(fill = NA, size = 1.5),
        axis.ticks = element_line(size = 0.45)) +
  labs(x = 'breakpoints [genomic coordinates in Mb]',
       y = 'density')
ggsave_golden(filename = 'Figures_original/breakpoints_Y.pdf', plot = breakpoint_plot, width = 7.5)  

##----------------+
## multiple segments
##----------------+
barplot(height = c(129, 118, 39, 3), width = c(0.2, 0.2, 0.2, 0.2), 
        space = c(0.1, 0.1, 0.1, 0.1), names.arg = c('loss', 'gain_loss', 'gain', 'relative_loss'), 
        col = 'white', 
        xlab = '', yaxt = 'n')
axis(side = 2, at = seq(0, 130, 20), las = 2, cex = 1.2, line = -0.5, lwd = 1.5)

## focal_gain_loss
focal_gain_loss = copy24[which(copy24$ID %in% loy$sample[which(loy$sample %in% foc2 & loy$classification == 'gain_loss')]), ]

gain_loss_all = data.frame()
for(i in unique(focal_gain_loss$ID)){
  if(focal_gain_loss$tcn.em[which(focal_gain_loss$ID == i)][1] > 
     focal_gain_loss$tcn.em[which(focal_gain_loss$ID == i)][2]){
    p_arm = 'yes'
  } else {
    p_arm = 'no'
  }
  out = data.frame(id = i,
                   p_arm = p_arm)
  gain_loss_all = rbind(gain_loss_all, out)
}

table(gain_loss_all$p_arm)
barplot(height = c(97, 21), width = c(0.2, 0.2), 
        space = c(0.1, 0.1), names.arg = c('Yp', 'Yq'), 
        col = 'white', 
        xlab = '', yaxt = 'n',
        ylim = c(0, 100))
axis(side = 2, at = seq(0, 120, 20), las = 2, cex = 1.2, line = 0, lwd = 1.5)


##----------------+
## Incidences across 
## tumor types
##----------------+
head(loy)
clinical = data$IMPACT_cohort
head(clinical)
head(loy)
prim_mets = merge(loy, clinical[,c('SAMPLE_ID', 'SAMPLE_TYPE', 'CANCER_TYPE', 'CANCER_TYPE_DETAILED')],
                  by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)
for(i in unique(prim_mets$CANCER_TYPE)){
  print(i)
  print(xtabs(~prim_mets$Y_call[which(prim_mets$CANCER_TYPE == i)] + prim_mets$SAMPLE_TYPE[which(prim_mets$CANCER_TYPE == i)]))
}

