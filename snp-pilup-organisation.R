clean()
.rs.restartR()
library(stringr)
library(dplyr)

DMP_bam = read.csv('~/Documents/MSKCC/dmp-2021/Genomics/key.txt', sep = '\t', header = F)
sample_pairing = read.csv('~/Documents/MSKCC/Subhi/CSF/Data/FINAL_samples/sample_pairing.txt', sep = '\t', header = F)
research = read.csv('~/Documents/MSKCC/Subhi/CSF/Data/files.txt', sep = '\t', header = F)
research$V1 = NULL

ids = read.csv('~/Documents/MSKCC/Subhi/CSF/Data/FINAL_samples/data_clinical_sample.txt', sep = '\t', skip = 3)
colnames(ids) = ids[1, ]
ids = ids[-1, c('SAMPLE_ID', 'PATIENT_ID')]
ids$PATIENT_ID = ifelse(grepl('C-', ids$PATIENT_ID), ids$PATIENT_ID, paste0('C-', ids$PATIENT_ID))


dmp_path = '/juno/res/dmpcollab/dmpshare/share/irb12_245/'
research_path = '/juno/res/ci/share/schultz/'

all_out = data.frame()
for(i in unique(ids$PATIENT_ID)){
  sample_subset = ids[which(ids$PATIENT_ID == i), ]
  for(j in 1:nrow(sample_subset)){
    id = sample_subset$SAMPLE_ID[j]
    if(grepl('_N90', id)) next
    else {
      print(id)
      if(any(grepl(id, DMP_bam$V1))){
        normal = grep(substr(id, start = 1, stop = 9), DMP_bam$V1, value = T)
        normal = grep('-N0', normal, value = T)
        normal = ifelse(length(normal) == 2, normal[1], normal)
        normal_id = str_split_fixed(normal, pattern = ',', 2)[,1]
        normal_path = str_split_fixed(normal, pattern = ',', 3)[,2]
        normal_path = paste0(dmp_path, substr(x = normal_path, start = 0, stop = 1), '/', 
                             substr(x = normal_path, start = 2, stop = 2), '/', normal_path, '.bam')
        
        
        path = grep(id, DMP_bam$V1, value = T)
        path = str_split_fixed(path[1], pattern = ',', 3)[,2]
        path = paste0(dmp_path, substr(x = path, start = 0, stop = 1), '/', 
                      substr(x = path, start = 2, stop = 2), '/', path, '.bam')
      } else if (any(grepl(pattern = id, x = research$V2))){
        path = grep(pattern = id, research$V2, value = T)
        path = grep(pattern = '.bam$', path, value = T)
        path = paste0(research_path, substring(path, 3))
        path = ifelse(length(path) == 2, path[1], path)
        normal_id = NA
        normal_path = NA
        
      } else {
        path = NA
        normal_id = NA
        normal_path = NA
      }
      
      out = data.frame(PATIENT_ID = i,
                       sample = c(id, normal_id),
                       path = c(path, normal_path))
      all_out = rbind(all_out, out)
    }
  }
}

all_out = all_out[!with(all_out, is.na(sample) & is.na(path)), ]
CSF_out = all_out %>% distinct(sample, path, .keep_all = TRUE)


##-----------------
## SAMPLE pairing if 
## no DMP is available
##-----------------
ID_match = data.frame()
for(i in unique(CSF_out$PATIENT_ID)){
  if(any(grepl(pattern = 'P-0', CSF_out$sample[which(CSF_out$PATIENT_ID == i)]))) {
    out = CSF_out[which(CSF_out$PATIENT_ID == i), ]
  }
  else {
    CSF_sub = CSF_out[which(CSF_out$PATIENT_ID == i), ]
    CSF_sample = as.character(CSF_sub[1, 'sample'])
    normal_csf = sample_pairing$V1[grepl(pattern = CSF_sample, sample_pairing$V2)]
    normal_csf = ifelse(length(normal_csf) > 1, normal_csf[1], normal_csf)
    normal_csf_path = grep(pattern = normal_csf, x = research$V2, value = T)
    normal_csf_path = grep(pattern = '.bam$', x = normal_csf_path, value = T)
    normal_csf_path = paste0(research_path, substring(normal_csf_path, first = 3))
    
    out = rbind(CSF_sub, data.frame(PATIENT_ID = i,
                                    sample = normal_csf,
                                    path = normal_csf_path))
  }
  ID_match = rbind(ID_match, out)
  rm(out)
}

ID_match = ID_match %>% distinct(sample, .keep_all = TRUE)

write.table(x = ID_match, file = '~/Documents/MSKCC/Subhi/CSF/Data/FINAL_samples/sample_match.txt', sep = '\t', row.names = F, quote = F)

#' out

