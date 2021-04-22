## Create IGV file for every FacetsY run, both on WES and IMPACT samples:
## 
## 09/04/2021:

# install local FacetsY and pctGCdata

require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('dplyr', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('data.table', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')

# load FACETS countsfiles for Prostate WES samples
WES.prostate.samples = read.csv('/juno/home/kreitzec/WES_Prostate/WES.prostate.samplepath.txt', header = F)
WES.prostate.samples = as.character(WES.prostate.samples$V1)

# function to convert FacetsY output into IGV
format_igv_seg = function(facets_output, sample_id, normalize = TRUE) {
  
  if (!all(c('snps', 'segs', 'dipLogR') %in% names(facets_output))) {
    stop(paste0('Input is missing segs, snps or dipLogR ojbect.'), call. = FALSE)
  }
  
  seg = group_by(facets_output$snps, chrom, seg) %>% 
    summarize(loc.start = min(maploc),
              loc.end = max(maploc)) %>% 
    ungroup() %>% 
    left_join(., select(facets_output$segs, chrom, seg, num.mark, seg.mean = cnlr.median),
              by = c('chrom', 'seg')) %>% 
    mutate(ID = sample_id) %>% 
    select(ID, chrom, loc.start, loc.end, num.mark, seg.mean)
  
  if (normalize) { seg = mutate(seg, seg.mean = seg.mean - facets_output$dipLogR) }
  data.frame(seg)
}


# loop through every file and run FacetsY

# Facets parameters:
cval.preprocess = 25
cval.postprocess = 150

WES.prostate.IGV_out = data.frame()
for(i in unique(WES.prostate.samples)){
  try({
    print(i)
    data.in = facetsY::readSnpMatrix(i)
    data.pre = facetsY::preProcSample(data.in, cval = cval.preprocess, gbuild = 'hg19', snp.nbhd = 250)
    data.process = facetsY::procSample(data.pre, cval = cval.postprocess)
    data.out = facetsY::emcncf(data.process)
    
    ID = substr(i, start = 72, stop = 107)
    
    #' create list
    FacetsY_output = list(snps = data.process$jointseg,
                          segs = data.out$cncf,
                          dipLogR = data.out$dipLogR)
    
    IGV_out = format_igv_seg(facets_output = FacetsY_output,
                             sample_id = ID)
    
    WES.prostate.IGV_out = rbind(WES.prostate.IGV_out, IGV_out)
    
  })
}

write.table(x = WES.prostate.IGV_out, file = '/juno/home/kreitzec/WES_Prostate/WES.prostate.IGV.out.seg', row.names = F, sep = '\t', quote = F)

## out