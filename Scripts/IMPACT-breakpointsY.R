## investigate the average breakpoints from IMPACT samples on the Y-chromosome

library(data.table)

## Input
IMPACT.prostate.samples = read.csv('/juno/home/kreitzec/WES_Prostate/WES.prostate.samplepath.txt', header = F)
IMPACT.prostate.samples = as.character(IMPACT.prostate.samples$V1)

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


breakpoints.panel.out_df = data.frame()
for(i in unique(IMPACT.prostate.samples)){
  input = data.table::fread(i)
  input = input[which(input$Chromosome == 'Y'), ]
  input.raw = nrow(input)
  input.modi = nrow(Facets.preprocess(input))
  ave.breaks = data.frame(id = i,
                          seq.accessed = input.raw,
                          breakpoints.used_facets = input.modi)
  
  breakpoints.panel.out_df = rbind(breakpoints.panel.out_df, ave.breaks)
  
}

write.table(breakpoints.panel.out_df, file = '/juno/home/kreitzec/WES_Prostate/IMPACT.breakpoints_out.txt', 
            row.names = F, sep = '\t')

