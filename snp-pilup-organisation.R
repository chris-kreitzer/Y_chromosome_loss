DMP_bam = read.csv('~/Documents/MSKCC/dmp-2021/Genomics/key.txt', sep = '\t', header = F)
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
        path = grep(id, DMP_bam$V1, value = T)
        path = str_split_fixed(path[1], pattern = ',', 3)[,2]
        path = paste0(dmp_path, substr(x = path, start = 0, stop = 1), '/', 
                      substr(x = path, start = 2, stop = 2), '/', path, '.bam')
      } else if (any(grepl(pattern = id, x = research$V2))){
        path = grep(pattern = id, research$V2, value = T)
        path = grep(pattern = '.bam$', path, value = T)
        path = paste0(research_path, substring(path, 3))
        path = ifelse(length(path) == 2, path[1], path)
        normal = NA
        
        
      } else {
        path = NA
        normal = NA
      }
      out = data.frame(PATIENT_ID = i,
                       sample = id,
                       path = path,
                       normal = normal)
      all_out = rbind(all_out, out)
    }
  }
}


