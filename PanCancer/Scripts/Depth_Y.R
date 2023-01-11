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
