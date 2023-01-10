##----------------+
## B1) There is no unique segment/probe/locus (marker) specifically 
## enriched in the Y chromosome loss cases in solid tumor tissue.
##----------------+

## start: 11/18/21
## extension: 11/22/21
## revision: 08/26/2022
## 
## chris-kreitzer

clean()
gc()
.rs.restartR()
setup(working.path =  '~/Documents/MSKCC/10_MasterThesis/')
source('Scripts/UtilityFunctions.R')



## Data input:
data_cnlr = vroom::vroom('Data_out/IMPACT/Cnlr_out.txt')
cohort = readRDS('Data_out/cohort_data.rds')
CN_IMPACT = read.csv('Data_out/IMPACT/IMPACT_copynumber_out.txt', sep = '\t')



#' determine which positions are covered most among all the MSK-IMPACT samples
FreqTable = as.data.frame(table(data_cnlr$maploc))
FreqTable$Var1 = as.numeric(as.character(as.factor(FreqTable$Var1)))

#' keep only postions which are shared at least among 10% of samples
FreqTable = FreqTable[which(FreqTable$Freq > 0.1 * length(unique(data_cnlr$ID))), ]
FreqTable$indx = seq(1, nrow(FreqTable), 1)
FreqTable$rel = FreqTable$Freq / length(unique(data_cnlr$ID))
FreqTable$label = ifelse(FreqTable$rel >= 0.5, 'label', 'not')

#' Visualization
postions_holistic_cohort = ggplot(FreqTable, aes(x = indx, y = Freq, size = rel)) +
  geom_point(color = "red", fill = alpha("red", 0.3), alpha = 0.5, shape = 21, stroke = 1) +
  geom_segment(aes(x = indx, xend = indx, y = 0, yend = Freq), size = 0.3) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     breaks = break.x, 
                     labels = labels.x) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, length(unique(data_cnlr$ID))),
                     breaks = c(0.25*length(unique(data_cnlr$ID)), 0.5*length(unique(data_cnlr$ID)), 0.75*length(unique(data_cnlr$ID))),
                     labels = c('25%', '50%', '75%')) +
  scale_size("Markers shared\namong samples [%]", 
             limits = c(0.1, 1), 
             labels = c('10%', '25%', "50%", '75%', '100%'),
             range = c(0.1,10)) +
  geom_hline(yintercept = c(0.25*length(unique(data_cnlr$ID)), 0.5*length(unique(data_cnlr$ID)), 0.75*length(unique(data_cnlr$ID))),
             linetype = 'dashed', size = 0.2, color = 'grey55') +
  theme_Y(base_size = 14, panel_border = T, base_line_size = 1) +
  geom_text(aes(label = ifelse(label == 'label', Var1, '')),
            hjust = 0, nudge_x = 5, size = 4) +
  labs(x = '', y = '')
                  
ggsave_golden(filename = 'Figures/Markers_allSamples.pdf', plot = postions_holistic_cohort, width = 12)

break.x = c(1, 50, 100, 150, 200, 250, 300, 350)
labels.x = c()
for(i in break.x){
  var = FreqTable$Var1[which(FreqTable$indx == i)]
  labels.x = c(labels.x, var)
}

labels.x = paste0(round(labels.x / 1000000, 1), ' Mb')



#' look only into those samples which have lost the Y-chromosome
IMPACT_cohort = cohort$IMPACT.cohort
Impact_Y_loss = IMPACT_cohort$SAMPLE_ID[which(IMPACT_cohort$Y_call == 'Y_chrom_loss')]
data_cnlr_loss = data_cnlr[which(data_cnlr$ID %in% Impact_Y_loss), ]

FreqTableLoss = as.data.frame(table(data_cnlr_loss$maploc))
FreqTableLoss$Var1 = as.numeric(as.character(as.factor(FreqTableLoss$Var1)))

FreqTableLoss = FreqTableLoss[which(FreqTableLoss$Freq > 0.1 * length(unique(data_cnlr_loss$ID))), ]
FreqTableLoss$indx = seq(1, nrow(FreqTableLoss), 1)
FreqTableLoss$rel = FreqTableLoss$Freq / length(unique(data_cnlr_loss$ID))
FreqTableLoss$label = ifelse(FreqTableLoss$rel >= 0.5, 'label', 'not')

#' Visualization:
postions_loss_cohort = ggplot(FreqTableLoss, aes(x = indx, y = Freq, size = rel)) +
  geom_point(color = "red", fill = alpha("red", 0.3), alpha = 0.5, shape = 21, stroke = 1) +
  geom_segment(aes(x = indx, xend = indx, y = 0, yend = Freq), size = 0.3) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     breaks = break.x, 
                     labels = labels.x) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, length(unique(data_cnlr_loss$ID))),
                     breaks = c(0.25*length(unique(data_cnlr_loss$ID)), 0.5*length(unique(data_cnlr_loss$ID)), 0.75*length(unique(data_cnlr_loss$ID))),
                     labels = c('25%', '50%', '75%')) +
  scale_size("Markers shared\namong samples [%]", 
             limits = c(0.1, 1), 
             labels = c('10%', '25%', "50%", '75%', '100%'),
             range = c(0.1,10)) +
  geom_hline(yintercept = c(0.25*length(unique(data_cnlr_loss$ID)), 0.5*length(unique(data_cnlr_loss$ID)), 0.75*length(unique(data_cnlr_loss$ID))),
             linetype = 'dashed', size = 0.2, color = 'grey55') +
  theme_Y(base_size = 14, panel_border = T, base_line_size = 1) +
  geom_text(aes(label = ifelse(label == 'label', Var1, '')),
            hjust = 0, nudge_x = 5, size = 4) +
  labs(x = '', y = '')


break.x = c(1, 50, 100, 150, 200, 250, 300, 350)
labels.x = c()
for(i in break.x){
  var = FreqTableLoss$Var1[which(FreqTableLoss$indx == i)]
  labels.x = c(labels.x, var)
}

labels.x = paste0(round(labels.x / 1000000, 1), ' Mb')


ggsave_golden(filename = 'Figures/Markers_LossSamples.pdf', plot = postions_loss_cohort, width = 12)

