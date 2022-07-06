#' Callable sequence; genes and features from the Y-chromosome
#' Read coverage (filtered and raw mapping) and per base coverage;
#' 
#' start: 06/29/2022
#' revision: 06/30/2022
#' revision: 07/06/2022
#' chris-kreitzer


clean()
gc()
.rs.restartR()
setwd('~/Documents/MSKCC/10_MasterThesis/Data/')
library(scales)
library(ggrepel)


#' Import data
Alignments = read.csv(file = 'AlignmentStats.txt', sep = '\t')
Regions = read.csv(file = 'Y_features_all_hg38.txt', sep = '\t')
Regions = as.data.frame(Regions)



#' merge regions with gene name
Alignments = merge(Alignments, Regions[,c('hgnc_symbol', 'Start_hg19')], 
                              by.x = 'start',  by.y = 'Start_hg19', all.x = T)


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
Alignment_summary$keep = ifelse(Alignment_summary$filtered_median / Alignment_summary$unfiltered_median > 0.6 &
                                  Alignment_summary$filtered_median > 20, 'yes', 'no')
Alignment_summary$keep[is.na(Alignment_summary$keep)] = 'no'


##-----------------
## Visualization
filter = Alignment_summary[,c('gene', 'filtered_median', 'keep')]
colnames(filter) = c('gene', 'value', 'keep')
filter$tag = 'filtered'
unfilter = Alignment_summary[,c('gene', 'unfiltered_median', 'keep')]
colnames(unfilter) = c('gene', 'value', 'keep')
unfilter$tag = 'unfiltered'

AlignmentsY = rbind(filter, unfilter)
AlignmentsY$subject = rep(seq(1, 430, 1), 2)
AlignmentsY$tag = factor(AlignmentsY$tag, levels = c('unfiltered', 'filtered'))

ggplot(AlignmentsY, aes(x = tag, y = value, group = subject, color = keep, label = ifelse(value > 50, gene, ""))) +
  geom_line(size = 0.35) +
  geom_text_repel(size = 3.5) +
  geom_point(size = 0.3) +
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



##-----------------------------------------------------------------------------
## Investigate the output:
baseCoverage = read.csv('baseCoverage.txt', sep = '\t')
baseCoverageGene = baseCoverage[baseCoverage$mean_cov != 'N/A', ]
baseCoverageGene$mean_cov = as.numeric(as.character(baseCoverageGene$mean_cov))

hist(baseCoverageGene$mean_cov[which(baseCoverageGene$Gene == 'ZFY')], nclass = 50)

Regions = readxl::read_excel(path = 'Y_gene_table.xlsx')
Regions = as.data.frame(Regions)


for(i in unique(baseCoverage$Gene)){
  
}



