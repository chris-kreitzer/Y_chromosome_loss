library(readr)
library(data.table)


readsnpmatrix = function(path, err.thresh = 10, del.thresh = 10){
  pileup = readr::read_csv(file = path, col_types = rep(c("c", "n","c", "n"), c(1,1,2,8)))
  # remove chr if present in Chrom
  if (grepl("chr", pileup$Chromosome[1])) {
    pileup$Chromosome = gsub("chr", "", pileup$Chromosome)
  }
  # remove loci where errors and deletions exceeded thresholds
  ii = which(pileup$File1E <= err.thresh & pileup$File1D <= del.thresh & 
               pileup$File2E <= err.thresh & pileup$File2D <= del.thresh & 
               !pileup$Chromosome %in% c('MT', 'chrM', 'Y', 'chrY'))
  
  rcmat = pileup[ii, 1:2]
  rcmat$NOR.DP = pileup$File1R[ii] + pileup$File1A[ii]
  rcmat$NOR.RD = pileup$File1R[ii]
  rcmat$TUM.DP = pileup$File2R[ii] + pileup$File2A[ii]
  rcmat$TUM.RD = pileup$File2R[ii]
  
  
  return(rcmat)
}
