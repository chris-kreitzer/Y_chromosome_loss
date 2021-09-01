## Look into Y-chromosome incidences among different cancer types
## Make an overview graph; as Bastien initially did (back in 2019)
## 
## 09/01/2021
## chris-kreitzer
## 


set.seed(112)
rm(list = ls())
.rs.restartR()


## Libraries and Input
cohortData = readRDS('Data_out/cohort_data.rds')
IMPACT_data = cohortData$IMPACT.cohort
WES_data = cohortData$WES.cohort
TCGA_data = cohortData$male_TCGA_cohort

#' make overview of available data resources (IMPACT, WES and TCGA)
#' IMPACT
top_20_cancer_type = table(IMPACT_data$CANCER_TYPE)
top_20_cancer_type = top_20_cancer_type[order(top_20_cancer_type, decreasing = T)]
top_20_cancer_type = top_20_cancer_type[1:20]

pdf(paste0('Figures/IMPACT_Top20_cancer_type.pdf'), width = 9, height = 6)
par(oma = c(0,10,0,0))
barplot(top_20_cancer_type[20:1], 
        horiz = T, 
        border = NA, 
        axes = F, 
        las = 2, 
        xlim = c(0, 1800), 
        cex.axis = 0.75,
        main = paste0('MSK-IMPACT (n = ', dim(IMPACT_data)[1], ')'))
axis(side = 1, at = c(0,500,1000,1500,2000), labels = c(0,500,1000,1500,2000))
dev.off()

#' WES
top_20_cancer_type = table(WES_data$Parental_Tumor_Type)
top_20_cancer_type = top_20_cancer_type[order(top_20_cancer_type, decreasing = T)]
top_20_cancer_type = top_20_cancer_type[1:20]

pdf(paste0('Figures/WES_Top20_cancer_type.pdf'), width = 9, height = 6)
par(oma = c(0,10,0,0))
barplot(top_20_cancer_type[20:1], 
        horiz = T, 
        border = NA, 
        axes = F, 
        las = 2, 
        xlim = c(0, 150), 
        cex.axis = 0.75,
        main = paste0('MSK-WES (n = ', dim(WES_data)[1], ')'))
axis(side = 1, at = c(0, 50, 100, 150), labels = c(0, 50, 100, 150))
dev.off()

#' TCGA
top_20_cancer_type = table(TCGA_data$type)
top_20_cancer_type = top_20_cancer_type[order(top_20_cancer_type, decreasing = T)]
top_20_cancer_type = top_20_cancer_type[1:20]

pdf(paste0('Figures/TCGA_Top20_cancer_type.pdf'), width = 9, height = 6)
par(oma = c(0,10,0,0))
barplot(top_20_cancer_type[20:1], 
        horiz = T, 
        border = NA, 
        axes = F, 
        las = 2, 
        xlim = c(0, 500), 
        cex.axis = 0.75,
        main = paste0('TCGA-WES (n = ', dim(TCGA_data)[1], ')'))
axis(side = 1, at = seq(0, 500, 100), labels = seq(0, 500, 100))
dev.off()



#' Y-chromosome loss among different cancer types; stratified according to primary and metastatic samples;
#' start with IMPACT
IMPACT_data = cohortData$IMPACT.cohort
IMPACT_top_20 = table(IMPACT_data$CANCER_TYPE)
IMPACT_top_20 = IMPACT_top_20[order(IMPACT_top_20, decreasing = T)]
IMPACT_top_20 = IMPACT_top_20[1:20]

IMPACT_incidence = data.frame()
for(i in unique(names(IMPACT_top_20))){
  n.all = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE %in% c('Primary', 'Metastasis') & IMPACT_data$CANCER_TYPE == i)])
  n.primary = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Primary' & IMPACT_data$CANCER_TYPE == i)])
  n.metastasis = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$SAMPLE_TYPE == 'Metastasis' & IMPACT_data$CANCER_TYPE == i)])
  
  primary.frequency = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$SAMPLE_TYPE == 'Primary' & IMPACT_data$CANCER_TYPE == i)]) / n.primary
  metastatic.frequency = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$SAMPLE_TYPE == 'Metastasis' & IMPACT_data$CANCER_TYPE == i)]) / n.metastasis
  cancer_median = length(IMPACT_data$SAMPLE_ID[which(IMPACT_data$Y_call == 'Y_chrom_loss' & IMPACT_data$CANCER_TYPE == i)]) / n.all
  cancer_summary = data.frame(CancerType = i,
                              Type = c('Primary', 'Metastasis'),
                              value = c(primary.frequency, metastatic.frequency),
                              median_cancerType = cancer_median)
  
  IMPACT_incidence = rbind(IMPACT_incidence, cancer_summary)
  rm(n.all)
  rm(primary.frequency)
  rm(metastatic.frequency)
  rm(cancer_median)
}

IMPACT_incidence

IMPACT_incidence_plot = ggplot(IMPACT_incidence, 
                               aes(x = reorder(CancerType, median_cancerType), 
                                   y = value, 
                                   fill = Type)) + 
                                 
  geom_bar(stat = 'identity', position = position_dodge(0.82), width = 0.8) +
  coord_flip()

IMPACT_incidence_plot
