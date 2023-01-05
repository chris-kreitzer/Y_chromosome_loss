##-----------------
## Tumor classification: 
## LOY:
## - complete
## - partial
## - relative
## - gain
## - gain-loss
##-----------------
##
## start: 09/21/2022
## revision: 01/05/2023
## chris-kreitzer


clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/')

#' data
CNA = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/04_Loss/010523/CopyNumberStates.txt', sep = '\t')


##-----------------
## FUNCTION:
##-----------------
Y_loss_call = function(data, sample.id){
  try({
    print(sample.id)
    data.sub = data[which(data$id == sample.id), ]
    data.sub = data.sub[which(data.sub$chrom == 24), ]
    if(nrow(data.sub) == 0) next
    else {
      #' expected copy number Y
      exp_Y = ifelse(unique(round(data.sub$ploidy)) %in% seq(0, 20, 2), 
                     unique(round(data.sub$ploidy / 2)), 
                     unique(floor(data.sub$ploidy / 2)))
      
      #' observed copy number Y
      if(nrow(data.sub) == 1){
        data.sub$classification = ifelse(data.sub$tcn.em > exp_Y, 'gain',
                                         ifelse(data.sub$tcn.em == exp_Y, 'wt', 'loss'))
      } else {
        for(i in 1:nrow(data.sub)){
          data.sub$classification[i] = ifelse(data.sub$tcn.em[i] > exp_Y, 'gain',
                                              ifelse(data.sub$tcn.em[i] == exp_Y, 'wt', 'loss'))
        }
      }
      
      data.sub$length = data.sub$end - data.sub$start
      Y_ratio = sum(data.sub$length[which(data.sub$tcn.em == 0)], na.rm = T) / sum(data.sub$length, na.rm = T)
      Y_call = ifelse(Y_ratio >= 0.5, 'Y_chrom_loss', 'intact_Y_chrom')
      Y_call = unique(Y_call)
      
      #' return values
      data.sub$Y_call = Y_call
      data.sub$ID = sample.id
      data.sub$expected = exp_Y
      data.sub = as.data.frame(data.sub, row.names = NULL)
      return(data.sub)
    }
  })
}

Y_evaluation = lapply(unique(CNA$id),
                      function(x) Y_loss_call(data = CNA, sample.id = x))




#' exclude samples with no call
Y_evaluation = Y_evaluation[!grepl(pattern = 'Error*', Y_evaluation$chrom), ]
Y_evaluation = Filter(function(x) length(x) > 1, Y_evaluation)
Y_evaluation = data.table::rbindlist(Y_evaluation)



##-----------------
## sample classification:
##-----------------
n_segments = table(Y_evaluation$id)
mono_segments = names(n_segments[which(n_segments == 1)])
double_segments = names(n_segments[which(n_segments == 2)])
tripe_segments = names(n_segments[which(n_segments == 3)])
quatro_segments = names(n_segments[which(n_segments == 4)])

##----------------+
## classification of 
## mono segments
##----------------+
mono_Y = Y_evaluation[which(Y_evaluation$id %in% mono_segments), ]
Y_out_mono = data.frame()
for(i in unique(mono_Y$id)){
  print(i)
  data.sub = mono_Y[which(mono_Y$id == i), c('cnlr.median', 'tcn.em', 
                                   'lcn.em', 'ID', 'purity', 
                                   'ploidy', 'classification',
                                   'Y_call', 'expected')]
  classi = data.sub$classification
  Y_expected = data.sub$expected
  
  out = ifelse(data.sub$tcn.em == 0, 'complete_loss',
               ifelse(data.sub$tcn.em < Y_expected & data.sub$tcn.em != 0, 'relative_loss',
                      ifelse(data.sub$tcn.em == Y_expected, 'wt',
                             ifelse(data.sub$tcn.em > Y_expected, 'gain', NA))))
    Y_out = data.frame(sample = i,
                       ploidy = data.sub$ploidy,
                       Y_expected = Y_expected,
                       Y_call = data.sub$Y_call,
                       classification = out)
    Y_out_mono = rbind(Y_out_mono, Y_out)
}


##----------------+
## double segments
##----------------+
double_Y = Y_evaluation[which(Y_evaluation$id %in% double_segments), ]
Y_out_double = data.frame()
for(i in unique(double_Y$id)){
  print(i)
  data.sub = mono_Y[which(mono_Y$id == i), c('cnlr.median', 'tcn.em', 
                                             'lcn.em', 'ID', 'purity', 
                                             'ploidy', 'classification',
                                             'Y_call', 'expected')]
  for(j in 1:nrow(data.sub)){
    Y_expected = unique(data.sub$expected)
    dominant = 
    
    out = ifelse(data.sub$tcn.em == 0, 'complete_loss',
                 ifelse(data.sub$tcn.em < Y_expected & data.sub$tcn.em != 0, 'relative_loss',
                        ifelse(data.sub$tcn.em == Y_expected, 'wt',
                               ifelse(data.sub$tcn.em > Y_expected, 'gain', NA))))
    Y_out = data.frame(sample = i,
                       ploidy = data.sub$ploidy,
                       Y_expected = Y_expected,
                       Y_call = data.sub$Y_call,
                       classification = out)
    Y_out_mono = rbind(Y_out_mono, Y_out)
  }
  }
  classi = data.sub$classification
  Y_expected = data.sub$expected
  
  out = ifelse(data.sub$tcn.em == 0, 'complete_loss',
               ifelse(data.sub$tcn.em < Y_expected & data.sub$tcn.em != 0, 'relative_loss',
                      ifelse(data.sub$tcn.em == Y_expected, 'wt',
                             ifelse(data.sub$tcn.em > Y_expected, 'gain', NA))))
  Y_out = data.frame(sample = i,
                     ploidy = data.sub$ploidy,
                     Y_expected = Y_expected,
                     Y_call = data.sub$Y_call,
                     classification = out)
  Y_out_mono = rbind(Y_out_mono, Y_out)
}





write.table(x = Y_out_all, file = '~/Documents/MSKCC/10_MasterThesis/Data/04_Loss/Categorial_Classification_Y.txt', sep = '\t', row.names = F)



#' out