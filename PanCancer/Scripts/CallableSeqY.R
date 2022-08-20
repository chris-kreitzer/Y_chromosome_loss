#' Callable sequence; DNA elements from the Y-chromosome
#' Read coverage (filtered and raw mapping) and per base coverage;
#' 
#' start: 06/29/2022
#' revision: 06/30/2022
#' revision: 07/06/2022
#' revision: 07/13/2022
#' revision: 08/17/2022
#' 
#' chris-kreitzer


clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/Data/')
library(scales)
library(ggrepel)


#' Import data
Alignments = read.csv(file = '01_Coverage_Depth/AlignmentStats.txt', sep = '\t')
Regions = read.csv(file = '01_Coverage_Depth/Y_features_annotated_hg19.txt', sep = '\t')
Regions = as.data.frame(Regions)

#' merge regions with gene name
Alignments = merge(Alignments, Regions[,c('hgnc_symbol', 'Start_hg19')], 
                              by.x = 'start',  by.y = 'Start_hg19', 
                   all.x = T)


Alignment_summary = data.frame()
for(i in unique(Alignments$hgnc_symbol)){
  gene_filtered = exp(mean(log(Alignments$records[which(Alignments$hgnc_symbol == i & Alignments$tag == 'filtered')])))
  gene_unfiltered = exp(mean(log(Alignments$records[which(Alignments$hgnc_symbol == i & Alignments$tag == 'unfiltered')])))
  out_f = data.frame(gene = i,
                     filtered_geo = gene_filtered,
                     filtered_median = median(Alignments$records[which(Alignments$hgnc_symbol == i & Alignments$tag == 'filtered')]),
                     unfiltered_geo = gene_unfiltered,
                     unfiltered_median = median(Alignments$records[which(Alignments$hgnc_symbol == i & Alignments$tag == 'unfiltered')]))
  
  Alignment_summary = rbind(Alignment_summary, out_f)
}


#' if more than 50% of the read drop; exclude gene from analysis
Alignment_summary$keep = ifelse(Alignment_summary$filtered_median > 20, 'yes', 'no')
Alignment_summary$keep[is.na(Alignment_summary$keep)] = 'no'
Alignment_summary = Alignment_summary[!is.na(Alignment_summary$gene), ]

#' save output:
Alignment_out = Alignment_summary[which(Alignment_summary$keep == 'yes'), ]
keep = as.character(unique(y$gene))
Alignment_out = Alignment_out[which(Alignment_out$gene %in% keep), ]
Alignment_out = merge(Alignment_out, Regions[,c('hgnc_symbol', 'gene_biotype')], by.x = 'gene', by.y = 'hgnc_symbol', all.x = T)
write.table(x = Alignment_out, file = '~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/DNA_elements_uniquely_covered.txt', sep = '\t', row.names = F, quote = F)
# write.table(Alignment_summary, file = 'signedOut/Alignment_summary.txt', sep = '\t', row.names = F, quote = F)



##-----------------
## Visualization
filter = Alignment_summary[,c('gene', 'filtered_median', 'keep')]
colnames(filter) = c('gene', 'value', 'keep')
filter$tag = 'filtered'
unfilter = Alignment_summary[,c('gene', 'unfiltered_median', 'keep')]
colnames(unfilter) = c('gene', 'value', 'keep')
unfilter$tag = 'unfiltered'

AlignmentsY = rbind(filter, unfilter)
AlignmentsY$subject = rep(seq(1, nrow(Alignment_summary), 1), 2)
AlignmentsY$tag = factor(AlignmentsY$tag, levels = c('unfiltered', 'filtered'))

ggplot(AlignmentsY, aes(x = tag, y = value, group = subject, color = keep, label = ifelse(keep == 'yes', gene, ""))) +
  geom_line(size = 0.65) +
  geom_text_repel(size = 3.5) +
  #geom_point(size = 0.3) +
  geom_hline(yintercept = c(1, 10, 100, 1000, 10000), color = 'grey30', size = 0.25, linetype = 'dashed') +
  scale_color_manual(values = c('no' = 'red', 'yes' = 'black'), name = '') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e0, 1e4)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        aspect.ratio = 1.8,
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black', size = 2),
        legend.position = 'none') +
  labs(y = '# sequence reads [median]', 
       x = paste0(length(unique(AlignmentsY$gene[which(AlignmentsY$keep == 'yes')])),
                  '/', length(unique(AlignmentsY$gene))),
       title = paste0('n = ', length(unique(Alignments$SampleID)), '; before and after filtering'))
                                                             



##-----------------
## Cross-Check

#' Pass


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
Master = read.csv('Coverage_n50.txt', sep = '\t')
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



###############################################################################
#' per base coverage of kept genes; from countfiles (snp-pileup)
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

##-----------------
## Aggregate length to be assessed with MSK-IMPACT sequencing
pileup = facetsY::readSnpMatrix(filename = '~/Desktop/mnt/ATMcountdata/countsMerged____P-0036764-T01-IM6_P-0036764-N01-IM6.dat.gz')
Y_chromosome = pileup[which(pileup$Chromosome == 'Y' & 
                               pileup$Position >= 2654550 & 
                               pileup$Position <= 28000000), ]

#' exclude PCDH11Y and centromeric region
Y_chromosome = Y_chromosome[with(Y_chromosome, !((Position %in% 4922131:5612269))), ]
Y_chromosome = Y_chromosome[with(Y_chromosome, !((Position %in% 10500000:14000000))), ]
Y_chromosome = Y_chromosome[which(Y_chromosome$NOR.DP > 20), ]
(max(Y_chromosome$Position) - min(Y_chromosome$Position)) / 59373566


##-----------------------------------------------------------------------------
## Investigate the output:
baseCoverage = read.csv('baseCoverage.txt', sep = '\t')
baseCoverageGene = baseCoverage[baseCoverage$mean_cov != 'N/A', ]
baseCoverageGene$mean_cov = as.numeric(as.character(baseCoverageGene$mean_cov))

hist(baseCoverageGene$mean_cov[which(baseCoverageGene$Gene == 'ZFY')], nclass = 50)
