## Investigate IMPACT panel data; binary loss of Y-chromosome
## Data was fetched from Juno; downstream analysis will take place here
## I have already discussed the parameters we used;
## 
## Start: 08/31/2021
## chris-kreitzer
## 


set.seed(112)
rm(list = ls())
.rs.restartR()


## Libraries and Input
library(data.table)
IMPACT_copyNumber = read.csv('Data_out/IMPACT/IMPACT_copynumber_out.txt', sep = '\t')
WES_copyNumber = read.csv('Data_out/WES/WES_copyNumber_out.txt', sep = '\t')
ID_matches = read.csv('Data_in/WES_cBio_ID.match.txt', sep = '\t')
cohortData = readRDS(file = 'Data_out/cohort_data.rds')


## Processing
Y_loss_call = function(data, sample.id){
  data.sub = data[which(data$ID == sample.id), ]
  data.sub = data.sub[which(data.sub$chrom == 24), ]
  if(nrow(data.sub) == 0){
    return(data.frame(cbind(Y_call = NA, 
                            sample.id,
                            length = NA,
                            ploidy = NA,
                            purity = NA,
                            min.segment = NA,
                            max.segment = NA)))
  } else {
    data.sub$length = data.sub$end - data.sub$start
    Y_ratio = sum(data.sub$length[which(data.sub$tcn.em == 0)], na.rm = T) / sum(data.sub$length, na.rm = T)
    
    # return values:
    purity = unique(data.sub$purity)
    length = unique(sum(data.sub$length, na.rm = T))
    ploidy = unique(data.sub$ploidy)
    min.segment = unique(min(data.sub$start))
    max.segment = unique(max(data.sub$end))
    
    Y_call = ifelse(Y_ratio >= 0.5, 'Y_chrom_loss', 'intact_Y_chrom')
    Y_call = unique(Y_call)
    
    sample.id = unique(sample.id)
    
    return(data.frame(cbind(Y_call, 
                            sample.id,
                            length,
                            ploidy,
                            purity,
                            min.segment = min.segment,
                            max.segment = max.segment)))
  }
}


#' WES binary calls out
WES_out = lapply(unique(WES_copyNumber$ID), function(x) Y_loss_call(data = WES_copyNumber, sample.id = x))
WES_out = rbindlist(WES_out, fill = T)

#' IMPACT binary calls out
IMPACT_out = lapply(unique(IMPACT_copyNumber$ID), function(x) Y_loss_call(data = IMPACT_copyNumber, sample.id = x))
IMPACT_out = rbindlist(IMPACT_out, fill = T)


#' update the cohort data;
#' IMPACT cohort
IMPACT_cohort = cohortData$IMPACT.cohort
IMPACT_cohort = merge(IMPACT_cohort, IMPACT_out[, c('Y_call', 'sample.id', 'ploidy', 'purity')],
                      by.x = 'SAMPLE_ID', by.y = 'sample.id', all.x = T)

IMPACT_cohort$TUMOR_PURITY = as.numeric(as.character(IMPACT_cohort$TUMOR_PURITY))
IMPACT_cohort$purity = as.numeric(as.character(IMPACT_cohort$purity))
IMPACT_cohort$ploidy = as.numeric(as.character(IMPACT_cohort$ploidy))

#' WES data
WES_cohort = cohortData$WES.cohort
WES_out$sample.id = substr(WES_out$sample.id, start = 3, stop = 17)
WES_cohort = merge(WES_cohort, WES_out[, c('Y_call', 'sample.id', 'ploidy', 'purity')],
                   by.x = 'CMO_Sample_ID', by.y = 'sample.id', all.x = T)

WES_cohort$purity = as.numeric(as.character(WES_cohort$purity))
WES_cohort$ploidy = as.numeric(as.character(WES_cohort$ploidy))

#' ID matches
ID_matches = ID_matches

cohortData = list(IMPACT.cohort = IMPACT_cohort,
                  WES.cohort = WES_cohort,
                  ID_matches = ID_matches)


saveRDS(cohortData, file = 'Data_out/cohort_data.rds')

#' out


