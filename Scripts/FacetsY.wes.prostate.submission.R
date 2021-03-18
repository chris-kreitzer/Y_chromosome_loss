# install local FacetsY and pctGCdata

require('pctGCdata', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('facetsY', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')

# load FACETS countsfiles for Prostate WES samples
WES.prostate.samples = read.csv('/juno/home/kreitzec/WES_Prostate/WES.prostate.samplepath.txt', header = F)
WES.prostate.samples = as.character(WES.prostate.samples$V1)

# loop through every file and run FacetsY
# rbind output into comprehensive cncf file (downstream analysis on local PC)

# Facets parameters:
cval.preprocess = 25
cval.postprocess = 150

Prostate.out_df = data.frame()
flags.out = c()
Cnlr.prostate.out_df = data.frame()

for(i in unique(WES.prostate.samples)){
  try({

	data.in = facetsY::readSnpMatrix(i)
  	data.pre = facetsY::preProcSample(data.in, cval = cval.preprocess, gbuild = 'hg19')
  	data.process = facetsY::procSample(data.pre, cval = cval.postprocess)
  	data.out = facetsY::emcncf(data.process)
        
        ID = substr(i, start = 72, stop = 107)
	purity = data.out$purity
	purity = ifelse(is.na(purity), 0, purity)

  	data.return = data.out$cncf
  	data.return$ID = ID
  	data.return$purity = purity
  	data.return$ploidy = data.out$ploidy
  	
	# catch flags, if there are some
	facets.flags = data.out$emflags
	if(!is.null(facets.flags)){
		flags.out = c(flags.out, paste0(facets.flags, ';', ID))
	}

  	Prostate.out_df = rbind(Prostate.out_df, data.return)

       # fetch cnlr values for chrY; orthogonal validation
       data.cnlr = data.process$jointseg
       data.cnlr = data.cnlr[which(data.cnlr$chrom == 24), ]
       data.cnlr$ID = ID

       Cnlr.prostate.out_df = rbind(Cnlr.prostate.out_df, data.cnlr)
     })
}

write.table(Prostate.out_df, file = '/juno/home/kreitzec/WES_Prostate/WES.Prostate.processed.out.txt', row.names = F, sep = '\t')
write.table(flags.out, file = '/juno/home/kreitzec/WES_Prostate/WES_errorflags.txt', row.names = F)
write.table(Cnlr.prostate.out_df, file = '/juno/home/kreitzec/WES_Prostate/Cnlr.WES.prostate.out.txt', row.names = F, sep = '\t')
