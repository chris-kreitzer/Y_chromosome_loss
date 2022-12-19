##----------------+
## Chromosome Y coverage
## starting from IMPACT pileups
## which genes (features) are covered
## assess the depth and breadth
##----------------+

## start: 09/14/2022
## revision: 12/19/2022
## chris-kreitzer

## juno version 
## n=22,320 IMPACT cases

clean()
gc()
setwd('~/Documents/MSKCC/10_MasterThesis/')
cohort = read.csv('Data/00_CohortData/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')

##----------------+
## Part1: 
## Feature coverage across 
## n=429 elements; 
## CountMatrices (snp-pileup) 
## average-depth; features to retain
##----------------+

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


##-----------------
## Part2:
## Downstream analysis;
##-----------------
## 09/16/2022
## 12/19/2022

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/10_MasterThesis/')
coverage = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/01_Coverage_Depth/breath_out.txt', sep = '\t')
coverage$Sample[which(unique(duplicated(coverage$Sample)))]

n = length(unique(coverage$Sample))
gene_coverage = data.frame()
for(i in unique(coverage$gene)){
  gene = i
  fraction = length(unique(coverage$Sample[which(coverage$gene == i)])) / n
  out = data.frame(gene = gene,
                   fraction = fraction * 100)
  gene_coverage = rbind(gene_coverage, out)
}

##-----------------
## only include genes which 
## are covered in at least 
## 75% of the cohort
##-----------------
gene_coverage_final = gene_coverage[which(gene_coverage$fraction >= 74), ]
gene_coverage_final = gene_coverage_final[order(gene_coverage_final$fraction, decreasing = T), ]
gene_coverage_final$gene = factor(gene_coverage_final$gene, levels = gene_coverage_final$gene)
source('~/Documents/MSKCC/10_MasterThesis/Scripts/plot_theme.R')
write.table(x = gene_coverage_final, file = 'Data/01_Coverage_Depth/Gene_Coverage_Final.txt', sep = '\t', row.names = F)

coverage_IMPACT = ggplot(gene_coverage_final, aes(x = gene, y = fraction)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = seq(25, 75, 25), color = 'white') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = 'Genomic features covered\n by MSK-IMPACT assay [%]')

ggsave(filename = '~/Documents/MSKCC/10_MasterThesis/Figures_original/Feature_Coverage_IMPACT.pdf', plot = coverage_IMPACT, width = 8)


##-----------------
## Stratify by gene panel
##-----------------
sample_ids = unique(coverage$Sample)
panel = data.frame(sample = sample_ids,
                   panel = NA)
panel$panel = ifelse(grepl('IM3', panel$sample), 'IM3',
                     ifelse(grepl('IM5', panel$sample), 'IM5',
                            ifelse(grepl('IM6', panel$sample), 'IM6',
                                   ifelse(grepl('IM7', panel$sample), 'IM7', NA))))
                     

IM3 = 1184
IM5 = 4453
IM6 = 13629
IM7 = 3046

gene_panel = data.frame()
for(i in unique(coverage$gene)){
  for(j in unique(coverage$panel)){
    gene = i
    if(j == 'IM3'){
      fraction = (length(coverage$Sample[which(coverage$gene == i & coverage$panel == j)]) / IM3) * 100
    } else if (j == 'IM5'){
      fraction = (length(coverage$Sample[which(coverage$gene == i & coverage$panel == j)]) / IM5) * 100
    } else if (j == 'IM6'){
      fraction = (length(coverage$Sample[which(coverage$gene == i & coverage$panel == j)]) / IM6) * 100
    } else {
      fraction = (length(coverage$Sample[which(coverage$gene == i & coverage$panel == j)]) / IM7) * 100
    }
    
    out = data.frame(gene = gene,
                     panel = j,
                     fraction = fraction)
    gene_panel = rbind(gene_panel, out)
  }
}

##-----------------
## consider 17 elements from before
##-----------------
gene_panel = gene_panel[which(gene_panel$gene %in% gene_coverage_final$gene), ]
gene_panel$gene = factor(gene_panel$gene, levels = gene_coverage_final$gene)
gene_panel$panel = factor(gene_panel$panel, levels = c('IM3', 'IM5', 'IM6', 'IM7'))
write.table(x = gene_panel, file = 'Data/01_Coverage_Depth/Feature_coverage_PANEL.txt', sep = '\t', row.names = F)


Feature_Panel = ggplot(gene_panel, aes(x = gene, y = fraction, fill = panel)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = c('IM3' = '#BABAB9',
                               'IM5' = '#898999',
                               'IM6' = '#5B5F75',
                               'IM7' = '#3B455F')) +
  geom_hline(yintercept = seq(25, 75, 25), color = 'white') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = 'Genomic features covered\n by MSK-IMPACT assay [%]')

ggsave(filename = '~/Documents/MSKCC/10_MasterThesis/Figures_original/Feature_Coverage_IMPACT_PANEL.pdf', plot = Feature_Panel, width = 8)



##-----------------
## Breadth of coverage
##-----------------
library(dplyr)
library(ggrepel)

coverage_final = coverage[which(coverage$gene %in% gene_coverage_final$gene), ]
depth_breadth_summary = coverage_final %>%
  group_by(gene) %>% 
  summarize(mean_depth = mean(ave_depth),
            mean_breadth = mean(breadth),
            sd_depth = sd(ave_depth),
            sd_breadth = sd(breadth))

write.table(depth_breadth_summary, file = 'Data/01_Coverage_Depth/Depth_Breadth_summary.txt', sep = '\t', row.names = F)

depth_plot = ggplot(depth_breadth_summary, aes(x = mean_breadth, y = mean_depth)) +
  geom_jitter(size = 2) +
  geom_pointrange(aes(ymin = mean_depth - sd_depth, ymax = mean_depth + sd_depth)) + 
  geom_pointrange(aes(xmin = mean_breadth - sd_breadth, xmax = mean_breadth + sd_breadth)) +
  geom_text_repel(aes(label = gene), hjust = 0, vjust = 1) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100)) +
  theme_std(base_size = 14, base_line_size = 1) +
  theme(aspect.ratio = 1) +
  labs(x = 'Breadth of feature coverage [%]',
       y = 'Depth of feature coverage')
ggsave(filename = 'Figures_original/Depth_Breadth.pdf', plot = depth_plot, device = 'pdf', width = 8)


