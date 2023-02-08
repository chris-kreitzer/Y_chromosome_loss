library(data.table)
low_purity_samples = read.csv('~/Master/Low_purity[0:0.29]_samples.txt', sep = '\t')
low_purity_samples = low_purity_samples[1:nrow(low_purity_samples), ]

depth_Y = function(samples){
  try({
    sample_out = data.frame()
    data.in = data.table::fread(input = samples, sep = ',')
    data.in = data.in[which(data.in$Chromosome == 'Y'), ]
    data.in$NOR.DP = data.in$File1A + data.in$File1R
    data.in$TUM.DP = data.in$File2A + data.in$File2R
    data.in = data.in[which(data.in$NOR.DP >= 20), ]
    sample = basename(samples)
    print(sample)
    panel = ifelse(grepl(pattern = 'IM3', sample), 'IM3',
                   ifelse(grepl(pattern = 'IM5', sample), 'IM5',
                          ifelse(grepl(pattern = 'IM6', sample), 'IM6',
                                 ifelse(grepl(pattern = 'IM7', sample), 'IM7', 'NA'))))
    
    out = data.frame(Sample = sample,
                     panel = panel,
                     average_depth_TUM = mean(data.in$TUM.DP, na.rm = T),
                     average_depth_NOR = mean(data.in$NOR.DP, na.rm = T))
    out  
  })
}

xx = lapply(unique(low_purity_samples), function(x) depth_Y(samples = x))
xx = Filter(function(x) length(x) > 1, xx)
xx = data.table::rbindlist(xx)
write.table(xx, file = '~/Master/Average_Depth_lowPurity.txt', sep = '\t', row.names = F, quote = F)


##-- using MADSEQ output in normals
cohort = readRDS('Data/00_CohortData/Cohort_071322.rds')
cohort = cohort[which(cohort$Study_include == 'yes'), ]
lowp = cohort$SAMPLE_ID[which(cohort$purity <= 0.3)]
highp = cohort$SAMPLE_ID[which(cohort$purity >= 0.8)]

files = list.files('Data/03_Mosaicism/NormalizedDepth/', full.names = T)

pp = function(samples){
  try({
    print(samples)
    if(any(grepl(pattern = samples, x = files))){
      data.in = data.table::fread(input = files[grep(pattern = samples, x = files)], sep = '\t')
      data.in = data.in[which(data.in$seqnames == 'chrY')]
      data.in$sample = samples
    } else next
  })
  data.in
}


yy = lapply(lowp, function(x) pp(samples = x))
yy = Filter(function(x) length(x) > 1, yy)
yy = data.table::rbindlist(yy)


xx = lapply(highp, function(x) pp(samples = x))
xx = Filter(function(x) length(x) > 1, xx)
xx = data.table::rbindlist(xx)

summary(xx$normed_depth)
summary(yy$normed_depth)



##-- low purity
count = read.csv('Data/01_Coverage_Depth/Average_Depth_lowPurity.txt', sep = '\t')
count$Sample = substr(x = count$Sample, start = 17, stop = 33)
lowpf = data.frame()
for(i in unique(yy$sample)){
  print(i)
  sample = i
  mean_p = mean(yy$normed_depth[which(yy$sample == i)], na.rm = T)
  out = data.frame(sample = sample,
                   mean_depth_f = mean_p)
  lowpf = rbind(lowpf, out)
}

low = merge(count, lowpf, by.x = 'Sample', by.y = 'sample', all.x = T)
low = low[!is.na(low$mean_depth_f), ]

cor.test(low$average_depth_NOR, low$mean_depth_f)
ggplot(low, aes(x = average_depth_NOR, y = mean_depth_f)) +
  geom_jitter() +
  scale_y_continuous(limits = c(0,600))




