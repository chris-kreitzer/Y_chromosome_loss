##-- WorkFlow: SCNA alteration detection in CSF samples
## 03/11/2024

##-- Dependencies
library(patchwork)
library(dplyr)
library(data.table)
library(readr)
library(mclust)
library(ggpubr)
source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/CSF/Scripts/FacetsPlot.R')
source('~/Documents/GitHub/CSF/Scripts/cf_Plot.R')
source('~/Documents/GitHub/CSF/Scripts/ReadSnpMatrix.R')
source('~/Documents/GitHub/CSF/Scripts/gene_closeup.R')
source('~/Documents/GitHub/CSF/Scripts/GMM.R')
source('~/Documents/GitHub/CSF/Scripts/ExonicMapping_CDKN2A.R')
setwd('~/Documents/MSKCC/11_CSF/07_CNA_Refit/')

##-- Defining genes of interest
GOIs = c('CDKN2A','CDK4','PTEN','EGFR','PDGFRA',
         'KIT','KDR','MET', 'RB1', 'MDM2', 'CDK6')

##-- Loss and Gain threshold
loss_threshold = 0.3922108
gain_threshold =  0.3524098
seed = 100

facets_plot = function(fit){
  suppressPackageStartupMessages(library(patchwork))
  i = cnlr_plot(facets_data = fit)
  ii = valor_plot(facets_data = fit)
  iii = icn_plot(facets_data = fit)
  iv = cf_plot(facets_data = fit)
  i / ii / iii / iv
}


##-----------------
## Data and prerequisite
dirs_all = list.files(path = '.', pattern = 'compiled_SCNA_output.rds$', full.names = T, recursive = T)
dirs_any = list.files(path = '.', pattern = '.rds$', full.names = T, recursive = T)
countmatrices = list.files(path = '.', pattern = '.dat.gz$', full.names = T, recursive = T)
cBIO = read.csv('~/Documents/MSKCC/11_CSF/00_Data/CNA_Data_cBIO.txt', sep = '\t')
reserve = cBIO[which(cBIO$Hugo_Symbol %in% GOIs), c('Hugo_Symbol', 'P.0024073.T03.IM7')]
colnames(reserve)[1] = 'gene'
colnames(reserve)[2] = 'cBIO_call'
reserve$cBIO_call = 0
NicSocci = read.csv('../00_Data/Nic_Socci_binary.txt', sep = '\t')
NicSocci = NicSocci[which(row.names(NicSocci) %in% GOIs), ]
NicSocci$Hugo_Symbol = row.names(NicSocci)
reserve_n = reserve
colnames(reserve_n)[2] = 'Nic_call'
reserve_n$Nic_call = 0
reserve_f = reserve
colnames(reserve_f)[2] = 'Facets_call'
##-----------------


sample = 'P-0047557-T01-IM6'


##-----------------
## Facets
##-----------------
dirs_all[grep(pattern = sample, dirs_all)]
#dirs_any[grep(pattern = sample, dirs_any)]

fit = readRDS(dirs_all[grep(pattern = sample, dirs_all)][1])$Facets
if(is.null(fit)){
  gene_level = reserve_f
} else {
  colnames(fit$segs)[10] = ifelse(colnames(fit$segs)[10] != 'start', 'start', colnames(fit$segs)[10])
  colnames(fit$segs)[11] = ifelse(colnames(fit$segs)[11] != 'end', 'end', colnames(fit$segs)[11])
  diplogr = ifelse(facets_fit_qc(fit)$facets_qc, fit$dipLogR, mean(fit$snps$cnlr))
  gene_level = facetsSuite::gene_level_changes(facets_output = fit, genome = 'hg19', algorithm = 'em')
  gene_level = gene_level[which(gene_level$gene %in% GOIs), ]
  gene_level$wgd = facets_fit_qc(fit)$wgd
  gene_level$map = paste(gene_level$wgd, gene_level$tcn, gene_level$mcn, gene_level$lcn, sep = ':')
  gene_level = gene_level[,c('gene','tcn', 'cn_state', 'filter', 'map')]
  gene_level = merge(gene_level[,c('gene', 'tcn', 'cn_state', 'map')],
                     copy_number_states[,c('numeric_call', 'map_string')],
                     by.x = 'map', by.y = 'map_string', all.x = T)
  gene_level$numeric_call = ifelse(is.na(gene_level$numeric_call) & gene_level$tcn > 5, 2, gene_level$numeric_call)
  gene_level = gene_level[,c('gene', 'numeric_call')]
  colnames(gene_level)[2] = 'Facets_call'
  for(i in 1:nrow(gene_level)){
    gene_level$Facets_call[i] = ifelse(facets_fit_qc(fit)$facets_qc, gene_level$Facets_call[i], 0)
  }
}


##-----------------
## CnLR/GMM pipeline
##-----------------
countmatrices[grep(pattern = sample, countmatrices)]
datain = readsnpmatrix(path = countmatrices[grep(pattern = sample, countmatrices)][1])
datain = as.data.frame(datain)
dense_out = facetsSuite::run_facets(read_counts = datain,
                                    cval = 100, 
                                    snp_nbhd = 50, 
                                    ndepth = 25, 
                                    seed = seed)

diplogr = ifelse(is.null(fit), mean(dense_out$snps$cnlr),
                 ifelse(facets_fit_qc(fit)$facets_qc, fit$dipLogR, mean(dense_out$snps$cnlr)))

GMM_gene = data.frame()
for(i in GOIs){
  print(i)
  try({
    gg = gene_closeup(data = dense_out, gene = i, diplogr = diplogr)
    goi_snps = gg$snps
    goi_snps = goi_snps[!duplicated(goi_snps$cnlr), ]
    
    snps = goi_snps$cnlr
    classi = Mclust(data = snps, verbose = F)
    
    sample_gmm = data.frame(cnlr = classi$data,
                            cluster = classi$classification,
                            uncertainty = classi$uncertainty)
    sample_gmm = cbind(sample_gmm, classi$z)
    
    sample_gmm = merge(sample_gmm, goi_snps[,c('maploc', 'cnlr', 'gene')],
                       by.x = 'cnlr', by.y = 'cnlr', all.x = T)
    sample_gmm = sample_gmm[which(sample_gmm$uncertainty <= 0.2), ]
    
    
    components = ifelse(classi$G > 1, TRUE, FALSE)
    probes = ifelse(length(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) > 10, TRUE, FALSE)
    ttable = table(sample_gmm$cluster[which(sample_gmm$gene == 'color')])
    cluster_high = as.integer(names(which(ttable/sum(ttable) > 0.6)))
    c_sel = ifelse(probes == FALSE, 'no_distinct_cluster_assignment_possible',
                   ifelse(length(cluster_high) > 0, cluster_high,
                          ifelse(length(which(ttable/sum(ttable) > 0.4)) > 1 &
                                   classi$parameters$mean[names(classi$parameters$mean) == which.max(as.integer(names(which(ttable/sum(ttable) > 0.4))))][[1]] >= diplogr &
                                   classi$parameters$mean[names(classi$parameters$mean) == which.max(as.integer(names(which(ttable/sum(ttable) > 0.4)))) - 1][[1]] <= diplogr,
                                 'no_distinct_cluster_assignment_possible',
                                 which.max(as.integer(names(which(ttable/sum(ttable) > 0.4)))))))
    cluster_mean = ifelse(c_sel == 'no_distinct_cluster_assignment_possible', diplogr,
                          ifelse(is.null(names(classi$parameters$mean)), classi$parameters$mean,
                                 classi$parameters$mean[which(names(classi$parameters$mean) == c_sel)][[1]]))
    
    
    ##-- loss or gain?
    loss = ifelse(cluster_mean < diplogr &
                    mean(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) < diplogr &
                    abs(mean(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) - diplogr) >= loss_threshold &
                    t.test(sample_gmm$cnlr[which(sample_gmm$gene == 'color')], 
                           sample_gmm$cnlr[which(sample_gmm$gene == 'no_color')])$p.value <= 0.01, -2, as.integer(0))
    gain = ifelse(cluster_mean > diplogr &
                    mean(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) > diplogr &
                    abs(mean(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) - diplogr) >= gain_threshold &
                    t.test(sample_gmm$cnlr[which(sample_gmm$gene == 'color')], 
                           sample_gmm$cnlr[which(sample_gmm$gene == 'no_color')])$p.value <= 0.01, 2, as.integer(0))
    
    gene_call = ifelse(probes == FALSE, 'NP',
                       ifelse(c_sel == 'no_distinct_cluster_assignment_possible', 'NP',
                              ifelse(loss == 0, gain,
                                     ifelse(gain == 0, loss, 0))))
    
    out = data.frame(gene = i,
                     GMM_call = gene_call)
    GMM_gene = rbind(GMM_gene, out)
    
  })
  
  rm(components, ttable, cluster_high, cluster_mean, 
     gg, goi_snps, c_sel, classi, snps,
     out, gene_call, gain, loss, i, probes, sample_gmm)
  gc()
  
}

if(dim(GMM_gene)[1] != 11){
  ssdd = setdiff(GOIs, GMM_gene$gene)
  GMM_gene = rbind(GMM_gene,
                   data.frame(gene = ssdd,
                              GMM_call = 'NP'))
} else {
  GMM_gene = GMM_gene
}


##-----------------
## CnLR/GMM pipeline
##-----------------
CnLR = readRDS(dirs_all[grep(pattern = sample, dirs_all)][1])$CnLR
CnLR_gene = CnLR[,c('gene', 'call')]
colnames(CnLR_gene)[2] = 'CnLR_call'

out1 = merge(gene_level, GMM_gene, by = 'gene', all.x = T)
out2 = merge(out1, CnLR_gene, by = 'gene', all.x = T)

##-----------------
## cBIO; clinical annotations
##-----------------
sample_c = gsub(pattern = '-', replacement = '.\\', x = sample)

if(sample_c %in% colnames(cBIO)){
  cBIO_sel = cBIO[which(cBIO$Hugo_Symbol %in% GOIs), c('Hugo_Symbol', sample_c)]
  colnames(cBIO_sel)[1] = 'gene'
  colnames(cBIO_sel)[2] = 'cBIO_call'
} else {
  cBIO_sel = reserve
}

out3 = merge(out2, cBIO_sel, by = 'gene', all.x = T)


##-----------------
## Nic_Socci calls
##-----------------

if(sample_c %in% colnames(NicSocci)){
  Nic_Sel = NicSocci[,c('Hugo_Symbol', sample_c)]
  colnames(Nic_Sel)[2] = 'Nic_call'
  colnames(Nic_Sel)[1] = 'gene'
} else {
  Nic_Sel = reserve_n
}


out4 = merge(out3, Nic_Sel, by = 'gene', all.x = T)
write.table(x = out4, file = paste0(sample, '/', sample, '_numeric_assignment_raw.txt'), sep = '\t', row.names = F, quote = F)


##-----------------
## Final call
##-----------------
out = out4
out$i_pipe = ifelse(out$GMM_call == out$CnLR_call, out$GMM_call,
                    ifelse(out$GMM_call == 0 & out$CnLR_call != 0, out$CnLR_call,
                           ifelse(out$CnLR_call != out$GMM_call, out$GMM_call, 
                                  ifelse(out$GMM_call != out$CnLR_call, out$GMM_call, 0))))
out$i_pipe_nic = ifelse(out$Nic_call == 0, out$i_pipe, out$Nic_call)
out$call = ifelse(out$cBIO_call != 0, out$cBIO_call,
                  ifelse(out$Facets_call == out$i_pipe_nic, out$Facets_call,
                         ifelse(out$Facets_call == 0 & out$i_pipe_nic != 0, out$i_pipe_nic,
                                ifelse(out$Facets_call != 0 & out$i_pipe_nic == 0, out$Facets_call,
                                       ifelse(out$Facets_call < 0 & out$i_pipe_nic < 0, out$Facets_call,
                                              ifelse(out$Facets_call > 0 & out$i_pipe_nic > 0, out$Facets_call,
                                                     ifelse(out$Facets_call < 0 & out$i_pipe_nic > 0 |
                                                              out$Facets_call > 0 & out$i_pipe_nic < 0, 'review', '0')))))))

all_zeros = apply(out, 2, function(x) all(x == 0))
zero_columns = names(out)[all_zeros]
out = out[, !all_zeros]
out = out[,!names(out) %in% c('GMM_call', 'CnLR_call')]
out = out[,!names(out) %in% c('Nic_call', 'i_pipe')]



##-----------------
out$call[which(out$gene == goi)] = 0
out = out[,c('gene', 'call')]
colnames(out)[2] = sample
write.table(x = out, file = paste0(sample, '/', sample, '_numeric_assignment_final.txt'), sep = '\t', row.names = F, quote = F)
binary = merge(binary, out, by = 'gene', all.x = T)
#write.table(x = binary, file = '~/Desktop/BINARY_CSF_03222024.txt', sep = '\t', row.names = F, quote = F)
##-----------------

if(any(grepl(pattern = 'review', x = out))){
  stop()
} else {
  rm(cBIO_sel, CnLR, CnLR_gene, fit, gene_level, GMM_gene, Nic_Sel,
     out1, out2, out3, sample, sample_c, 
     datain, dense_out, out, out4, diplogr,
     all_zeros, zero_columns, goi)
  gc()
  dev.off()
}







##-----------------
##-- manual review
out4      # where does the signal come from
goi = 'PDGFRA'
gg = gene_closeup(dense_out, gene = goi, diplogr = diplogr)
gg$plot
facetsSuite::closeup_plot(facets_data = fit, highlight_gene = goi, genome = 'hg19')
goi_snps = gg$snps
goi_snps = goi_snps[!duplicated(goi_snps$cnlr), ]
snps = goi_snps$cnlr
classi = Mclust(data = snps, verbose = F)
sample_gmm = data.frame(cnlr = classi$data,
                        cluster = classi$classification,
                        uncertainty = classi$uncertainty)
sample_gmm = cbind(sample_gmm, classi$z)
sample_gmm = merge(sample_gmm, goi_snps[,c('maploc', 'cnlr', 'gene')],
                   by.x = 'cnlr', by.y = 'cnlr', all.x = T)
sample_gmm = sample_gmm[which(sample_gmm$uncertainty <= 0.2), ]
components = ifelse(classi$G > 1, TRUE, FALSE)
probes = ifelse(length(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) > 10, TRUE, FALSE)
ttable = table(sample_gmm$cluster[which(sample_gmm$gene == 'color')])
cluster_high = as.integer(names(which(ttable/sum(ttable) > 0.6)))
classi$parameters$mean
c_sel = ifelse(probes == FALSE, 'no_distinct_cluster_assignment_possible',
               ifelse(length(cluster_high) > 0, cluster_high,
                      ifelse(length(which(ttable/sum(ttable) > 0.4)) > 1 &
                               classi$parameters$mean[names(classi$parameters$mean) == which.max(as.integer(names(which(ttable/sum(ttable) > 0.4))))][[1]] >= diplogr &
                               classi$parameters$mean[names(classi$parameters$mean) == which.max(as.integer(names(which(ttable/sum(ttable) > 0.4)))) - 1][[1]] <= diplogr,
                             'no_distinct_cluster_assignment_possible',
                             which.max(as.integer(names(which(ttable/sum(ttable) > 0.4)))))))

cluster_mean = ifelse(c_sel == 'no_distinct_cluster_assignment_possible', diplogr,
                      ifelse(is.null(names(classi$parameters$mean)), classi$parameters$mean,
                             classi$parameters$mean[which(names(classi$parameters$mean) == c_sel)][[1]]))
summary(sample_gmm$cnlr[which(sample_gmm$gene == 'color')])
loss = ifelse(cluster_mean < diplogr &
                mean(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) < diplogr &
                abs(mean(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) - diplogr) >= loss_threshold &
                t.test(sample_gmm$cnlr[which(sample_gmm$gene == 'color')], 
                       sample_gmm$cnlr[which(sample_gmm$gene == 'no_color')])$p.value <= 0.01, -2, as.integer(0))
gain = ifelse(cluster_mean > diplogr &
                mean(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) > diplogr &
                abs(mean(sample_gmm$cnlr[which(sample_gmm$gene == 'color')]) - diplogr) >= gain_threshold &
                t.test(sample_gmm$cnlr[which(sample_gmm$gene == 'color')], 
                       sample_gmm$cnlr[which(sample_gmm$gene == 'no_color')])$p.value <= 0.01, 2, as.integer(0))

gene_call = ifelse(probes == FALSE, 'NP',
                   ifelse(c_sel == 'no_distinct_cluster_assignment_possible', 'NP',
                          ifelse(loss == 0, gain,
                                 ifelse(gain == 0, loss, 0))))

gene_call

##-----------------

rm(cBIO_sel, CnLR, CnLR_gene, fit, gene_level, GMM_gene, Nic_Sel,
   out1, out2, out3, sample_gmm, i, probes, sample, sample_c, 
   classi, datain, dense_out, gg, goi_snps, out, out4, diplogr, c_sel, cluster_high, cluster_mean, components,
   all_zeros, goi, snps, ttable, zero_columns)
gc()
dev.off()



a) visual: no
##--
b) >1 cluster: yes
c) >10 probes with probability of > 80% assigned to specific cluster: yes
d) >60% of probes belonging to one cluster (< OR > diplogr - gain, loss): no (waterfall, probably indifferent, diffuse or just diploid)
e) ABS difference between gene_mean and diplogr greater or less than thresholds** 
  f) AND significant (t.test; considering local variance with ALPHA < 0.01). 