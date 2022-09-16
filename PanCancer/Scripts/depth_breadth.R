##-----------------
## Chromosome Y coverage
## starting from IMPACT pileups
## which genes (features) are covered
## assess the depth and breadth
##-----------------

## start: 09/14/2022
## chris-kreitzer

## juno version 
## n=22,320 IMPACT cases


library(data.table)
library(dplyr)

files = read.csv('/home/kreitzec/Master/Data/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')
files = as.character(files$counts_file)
annotations = read.csv('/home/kreitzec/Master/Data/Y_chromosome_annotation_hg19.txt', sep = '\t')

depth_breadth = function(x, annotation){
  try({
    sample_out = data.frame()
    data.in = data.table::fread(input = x, sep = ',')
    data.in = data.in[which(data.in$Chromosome == 'Y'), ]
    data.in$NOR.DP = data.in$File1A + data.in$File1R
    data.in = data.in[which(data.in$NOR.DP >= 20), ]
    sample = basename(x)
    print(sample)
    panel = ifelse(grepl(pattern = 'IM3', sample), 'IM3',
                   ifelse(grepl(pattern = 'IM5', sample), 'IM5',
                          ifelse(grepl(pattern = 'IM6', sample), 'IM6',
                                 ifelse(grepl(pattern = 'IM7', sample), 'IM7', 'NA'))))
 
    for(i in 1:nrow(annotation)){
      if(any(between(x = data.in$Position, left = annotation$Start_hg19[i], right = annotation$End_hg19[i]))){
        gene = annotation$Gene_hgnc_symbol[i]
      } else  next
      out = data.frame(Sample = sample,
                       panel = panel,
                       gene = gene,
                       ave_depth = mean(data.in$NOR.DP[which(between(x = data.in$Position, 
                                                                     left = annotation$Start_hg19[i], 
                                                                     right = annotation$End_hg19[i]))]),
                       breadth = (max(data.in$Position[which(between(x = data.in$Position, 
                                                                     left = annotation$Start_hg19[i], 
                                                                     right = annotation$End_hg19[i]))]) -
                                    min(data.in$Position[which(between(x = data.in$Position, 
                                                                       left = annotation$Start_hg19[i], 
                                                                       right = annotation$End_hg19[i]))])) /
                         (annotation$End_hg19[i] - annotation$Start_hg19[i]))
      sample_out = rbind(sample_out, out)
    }
    return(sample_out)
    rm(sample, data.in, panel)
  })
}

breath_out = lapply(files, function(x) depth_breadth(x = x, annotation = annotations))
breath_out = data.table::rbindlist(breath_out)

write.table(x = breath_out, file = '/home/kreitzec/Master/Data/breath_out.txt', sep = '\t', row.names = F)

#' out
 
