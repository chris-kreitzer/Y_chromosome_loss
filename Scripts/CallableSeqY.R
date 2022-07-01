#' Callable sequence; genes from the Y-chromosome
#' Read coverage (unique sequencing) and per base coverage
#' 
#' start: 06/29/2022
#' revsion: 06/30/2022
#' chris-kreitzer


clean()
gc()
.rs.restartR()

setwd('~/Documents/MSKCC/10_MasterThesis/Data/')

#' discard double file; cases
Alignments = read.csv(file = 'AlignmentStatsY.txt', sep = '\t')
Regions = readxl::read_excel(path = 'Y_gene_table.xlsx')
Regions = as.data.frame(Regions)

# files = names(which(table(Alignments$file) > 84))
# nonredundancy = data.frame()
# for(i in unique(Alignments$file)){
#   if(nrow(Alignments[which(Alignments$file == i), ]) > 84){
#     data_sub = Alignments[which(Alignments$file == i), ]
#     data_sub = data_sub[1:84, ]
#   } else {
#     data_sub = Alignments[which(Alignments$file == i), ]
#   }
#   nonredundancy = rbind(nonredundancy, data_sub)
# }
# 
# Alignments = nonredundancy


#' merge regions with gene name
Alignments = merge(Alignments, Regions[,c('Gene_name', 'Start_hg19')], 
                              by.x = 'start',  by.y = 'Start_hg19', all.x = T)


Alignment_summary = data.frame()
for(i in unique(Alignments$Gene_name)){
  gene_filtered = exp(mean(log(Alignments$records[which(Alignments$Gene_name == i & Alignments$tag == 'tumor_filtered')])))
  gene_unfiltered = exp(mean(log(Alignments$records[which(Alignments$Gene_name == i & Alignments$tag == 'tumor_unfiltered')])))
  out_f = data.frame(gene = i,
                     filtered_geo = gene_filtered,
                     filtered_median = median(Alignments$records[which(Alignments$Gene_name == i & Alignments$tag == 'tumor_filtered')]),
                     unfiltered_geo = gene_unfiltered,
                     unfiltered_median = median(Alignments$records[which(Alignments$Gene_name == i & Alignments$tag == 'tumor_unfiltered')]))
  
  Alignment_summary = rbind(Alignment_summary, out_f)
}

#' if more than 60% of the read drop; exclude gene from analysis
Alignment_summary$keep = ifelse(Alignment_summary$filtered_median / Alignment_summary$unfiltered_median > 0.6, 'keep', 'discard')
Alignment_summary$keep[is.na(Alignment_summary$keep)] = 'discard'

#' Visualization:
filter = Alignment_summary[,c('gene', 'filtered_median', 'keep')]
colnames(filter) = c('gene', 'value', 'keep')
filter$tag = 'filtered'
unfilter = Alignment_summary[,c('gene', 'unfiltered_median', 'keep')]
colnames(unfilter) = c('gene', 'value', 'keep')
unfilter$tag = 'unfiltered'

AlignmentsY = rbind(filter, unfilter)
AlignmentsY$subject = rep(seq(1, 42, 1), 2)
AlignmentsY$tag = factor(AlignmentsY$tag, levels = c('unfiltered', 'filtered'))

ggplot(AlignmentsY, aes(x = tag, y = value, group = subject, color = keep)) +
  geom_line(size = 0.35) +
  geom_text(aes(label = gene), size = 3, vjust = 0.5, hjust = 0.5) +
  geom_point() +
  scale_color_manual(values = c('discard' = 'red', 'keep' = 'black')) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e0, 1e4)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        aspect.ratio = 1) +
  labs(y = '# sequence reads [median]', x = '11/42', title = 'n = 1,756; before after filtering')
  



###############################################################################
#' per base coverage of kept genes; from countfiles (snp-pileup)
library(facets)
paths = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
cohort = read.csv('ProstateCohort_cbio.tsv', sep = '\t')
cohort = cohort[,c('Patient.ID', 'Sample.ID')]
cohort = cohort[sample(nrow(cohort), 1000, replace = F), ]
Regions = readxl::read_excel(path = 'Y_gene_table.xlsx')
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



##-----------------------------------------------------------------------------
## Investigate the output:
baseCoverage = read.csv('../Data/baseCoverage.txt', sep = '\t')
baseCoverageGene = baseCoverage[baseCoverage$mean_cov != 'N/A', ]
baseCoverageGene$mean_cov = as.numeric(as.character(baseCoverageGene$mean_cov))

hist(baseCoverageGene$mean_cov[which(baseCoverageGene$Gene == 'NLGN4Y')])

Regions = readxl::read_excel(path = 'Y_gene_table.xlsx')
Regions = as.data.frame(Regions)


for(i in unique(baseCoverage$Gene)){
  
}


a = facets::readSnpMatrix(filename = '~/Desktop/mnt/ATMcountdata/countsMerged____P-0029087-T02-IM6_P-0029087-N01-IM6.dat.gz')
a = a[which(a$Chromosome == 'Y'), ]
a = a[which(a$Position >= 6778727 & a$Position <= 6970011), ]

a


library("seqinr")

y = read.fasta(file = '~/Desktop/TGIF2LY.fasta')
y = y[[1]]
x = read.fasta(file = '~/Desktop/TGIF2LX.fasta')
x = x[[1]]
dotPlot(y, x, wsize = 1, nmatch = T, wstep = 3)


Regions = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/proteinCodingY.txt', sep = '\t', skip = 1)
files = read.csv('~/Desktop/AlignmentStatsY.txt', sep = '\t')
out = merge(files, Regions[, c('Gene_name', 'Start')], 
            by.x = 'start', by.y = 'Start')



all_out = data.frame()
for(i in unique(out$Gene_name)){
  gene_filtered = exp(mean(log(out$records[which(out$Gene_name == i & out$tag == 'tumor_filtered')])))
  gene_unfiltered = exp(mean(log(out$records[which(out$Gene_name == i & out$tag == 'tumor_unfiltered')])))
  out_f = data.frame(gene = i,
                   filtered = gene_filtered,
                   unfiltered = gene_unfiltered)
  all_out = rbind(all_out, out_f)
}

