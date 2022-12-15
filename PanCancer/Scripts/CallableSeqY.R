##----------------+
## Callable sequence; 
## DNA elements from the Y-chromosome
## Read coverage, and average
##----------------+
#' 
#' start: 06/29/2022
#' revision: 06/30/2022
#' revision: 07/06/2022
#' revision: 07/13/2022
#' revision: 08/17/2022
#' revision: 12/15/2022
#'
#' chris-kreitzer


clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/Data/')
library(scales)
library(ggrepel)


##-----------------
## Seq. coverage across Y-chromosome; n=50
files = list.files(path = '~/Desktop/mnt/ATMcountdata/', full.names = T)

master = c()
for(i in unique(files)){
  print(i)
  data.in = facetsY::readSnpMatrix(filename = i)
  data.in = data.in[which(data.in$Chromosome == 'Y' & data.in$NOR.DP > 5), ]
  master = c(master, unique(data.in$Position))
}

master = unique(master)

##-------
Master = data.frame(Chromosome = 'Y',
                    loci = master)
for(i in unique(files)){
  print(i)
  data.in = facetsY::readSnpMatrix(filename = i)
  data.in = data.in[which(data.in$Chromosome == 'Y' & data.in$NOR.DP > 5), ]
  Master = merge(Master, data.in[, c('Position', 'NOR.DP')], 
                 by.x = 'loci', 
                 by.y = 'Position', 
                 all.x = T)
  rm(data.in)
}

##-------
Master$mean = rowSums(Master[,3:ncol(Master)], na.rm = TRUE) / rowSums(!is.na(Master[ ,3:ncol(Master)]))

#' exclude PCDH11Y and centromeric region
Master = Master[which(Master$loci >= 2654550 & 
                        Master$loci <= 28000000), ]
Master = Master[with(Master, !((loci %in% 4922131:5612269))), ]
Master = Master[with(Master, !((loci %in% 10500000:14000000))), ]

##-------
Master = read.csv('Data/01_Coverage_Depth/Coverage_n50.txt', sep = '\t')
Master$seq = seq(1, nrow(Master), 1)

## Visualization
plot(Master$seq, Master$mean, 
     lwd = 1, 
     col = 'black', 
     type = 'S', 
     pch = 19,
     ylim = c(20, 600),
     yaxt = 'n',
     xaxt = 'n',
     xlab = '',
     ylab = '')
axis(side = 2, at = seq(0, 600, 100), lwd = 2, las = 2, line = 0.3)
axis(side = 1, at = seq(0, 2500, 200), las = 2)
for(i in seq(3, ncol(Master)-2, 1)){
  lines(Master$seq, Master[, i], col = 'grey75')
}
lines(Master$seq, Master$mean, lwd = 1.5, col = 'black', type = 'S', pch = 19)
box(lwd = 2)
mtext(text = 'Coverage', side = 2, line = 3, cex = 1.2)
mtext(text = 'Y chromosome coordinates [2.7-27.8 Mb]', side = 1, line = 1, cex = 1.2)


##----------------+
## per base coverage; 
##----------------+ 
library(facets)
paths = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
cohort = read.csv('ProstateCohort_cbio.tsv', sep = '\t')
cohort = cohort[,c('Patient.ID', 'Sample.ID')]
cohort = cohort[sample(nrow(cohort), 1000, replace = F), ]
Regions = read.csv('Y_features_all_hg38.txt', sep = '\t')
Regions = as.data.frame(Regions)


all_out = data.frame()
for(patient in 1:nrow(cohort)){
  try({
    path_count = grep(pattern = cohort$Sample.ID[patient], x = Facets_paths$counts_file, value = T)
    file_in = facets::readSnpMatrix(filename = path_count)
    file_in = file_in[which(file_in$Chromosome == 'Y'), ]
    
    for(gene in 1:nrow(Regions)){
      gene_selected = file_in[which(file_in$Position >= Regions$Start_hg19[gene] & file_in$Position <= Regions$End_hg19[gene]), ]
      gene_name = Regions$Gene_name[gene]
      
      if(any(gene_selected$TUM.DP > 30)){
        n = length(gene_selected$TUM.DP[which(gene_selected$TUM.DP > 30)])
        n_sum = sum(gene_selected$TUM.DP[which(gene_selected$TUM.DP > 30)])
        out = data.frame(Gene = gene_name,
                         Patient = cohort$Sample.ID[patient],
                         mean_cov = n_sum / n,
                         delta_5 = Regions$Start_hg19[gene] - min(gene_selected$Position[which(gene_selected$TUM.DP > 30)]),
                         delta_3 = Regions$End_hg19[gene] - max(gene_selected$Position[which(gene_selected$TUM.DP > 30)]))
      } else {
        out = data.frame(Gene = gene_name,
                         Patient = cohort$Sample.ID[patient],
                         mean_cov = 'N/A',
                         delta_5 = 'N/A',
                         delta_3 = 'N/A')
      }
      all_out = rbind(all_out, out)
    }
  })
}
