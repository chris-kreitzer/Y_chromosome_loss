# Investigate the overall coverage of the Y-chromosome from from IMPACT-Seq data

library(data.table)

bin.number = 5900 
start.chromosome = 1
end.chromosome = 59000000
Y.chromo.MasterFile_df = data.frame(bin = NA,
                                    loc = seq(start.chromosome, end.chromosome, 1))                        
bins.out = lapply(seq(1, bin.number, 1), 
                  function(x) rep(paste0('bin:', x), 
                                  nrow(Y.chromo.MasterFile_df) / bin.number))

bins.out = unlist(bins.out)
Y.chromo.MasterFile_df$bin = bins.out

Panel.Prostate.samples = read.csv('/juno/home/kreitzec/WES_Prostate/Panel.prostate.samplepath.txt', header = F)
Panel.Prostate.samples = as.character(Panel.Prostate.samples$V1)

min.function = function(data, bin){min(data$loc[which(data$bin == bin)])}
max.function = function(data, bin){max(data$loc[which(data$bin == bin)])}

#' count how many positions (bases) or covered > threshold in given bin
coverage.count.function = function(data, bin){
  n = length(data$Tumor_Coverage[which(data$Tumor_Coverage != 0 & data$bin == bin)])
}

#' select probes which overcome default Facets parameter; see github/facets
Facets.preprocess = function(data){
  # default input parameters
  err.thresh = Inf 
  del.thresh = Inf
  ndepth = 35
  ndepthmax = 1000
  
  data$Normal_Coverage = data$File1A + data$File1R
  data$Tumor_Coverage = data$File2A + data$File2R
  
  # discard non-fintite values
  ii = which(data$File1E <= err.thresh & 
               data$File1D <= del.thresh & 
               data$File2E <= err.thresh & 
               data$File2D <= del.thresh)
  
  data = data[ii, ]
  
  # discard low coverage positions
  depthN.keep = data$Normal_Coverage >= ndepth & 
    data$Normal_Coverage < ndepthmax
  
  data = data[depthN.keep, ]
}  


Y_coverage = data.frame()
for(i in sample(Panel.Prostate.samples, 10, replace = F)){
  print(i)
  input = data.table::fread(i)
  input = input[which(input$Chromosome == 'Y'), ]
  input.modi = Facets.preprocess(data = input)
  data.merge = merge(Y.chromo.MasterFile_df, 
                            input.modi[, c('Position', 'Normal_Coverage', 'Tumor_Coverage')], 
                            by.x = 'loc', 
                            by.y = 'Position',
                            all.x = T)
  
  data.merge[is.na(data.merge)] = 0
  message('Analyse bins now')
  data.merge.out = data.frame(bin = unique(data.merge$bin),
                              start = unlist(lapply(unique(data.merge$bin), function(x) min.function(data = data.merge, bin = x))),
                              end = unlist(lapply(unique(data.merge$bin), function(x) max.function(data = data.merge, bin = x))),
                              SNPs.covered = unlist(lapply(unique(data.merge$bin), function(x) coverage.count.function(data = data.merge, bin = x))))
  
  Y_coverage = rbind(Y_coverage, data.merge.out)
}
  
write.table(Y_coverage, file = '/juno/home/kreitzec/WES_Prostate/IMPACT-Y.coverage.txt', sep = '\t', row.names = F)

