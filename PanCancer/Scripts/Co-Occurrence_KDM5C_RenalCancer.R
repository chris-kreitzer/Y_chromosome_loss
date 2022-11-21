##----------------+
## gene-level GLM
## downstream analysis
## co-occurrence and mutual 
## exclusivity in Renal Cell
## Carcinoma
##----------------+

## start: 11/20/2022
## chris-kreitzer


library(ggrepel)
library(forcats)
library(tidyr)
library(plyr)
source('Scripts/plot_theme.R')
source('~/Documents/GitHub/MSKCC/Scripts/basicOncoprint.R')

gene_level_df = read.csv('Data/05_Association/gene_level/gene_level_full_out.txt', sep = '\t')
gene_level_df = gene_level_df[which(gene_level_df$variable == 'Y_call' & gene_level_df$p_adj <= 0.1), ]


ggplot(gene_level_df, aes(x = gene, y = estimate)) +
  geom_jitter(size = 1.5, position = position_dodge(width = 0.2)) +
  geom_pointrange(aes(ymin = estimate - std_err,
                      ymax = estimate + std_err),
                  size = 0.5) +
  geom_vline(xintercept = seq(1.5, 17, 1), linetype = 'dashed', color = 'grey35', size = 0.4) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey35', size = 0.4) +
  coord_flip() +
  geom_text_repel(aes(label = cancer_type)) +
  theme_std(base_size = 14) +
  theme(panel.border = element_rect(fill = NA, colour = 'black', size = 2)) +
  labs(x = '', y = 'log ODDS')
  



##----------------+
## Investigate KDM5C
## for Renal Cell Carcinoma
##----------------+
cohort = readRDS('Data/signedOut/Cohort_07132022.rds')
loy = cohort$IMPACT_Y_classification_final
clinical = cohort$IMPACT_clinicalAnnotation

loy = merge(loy, clinical[,c('SAMPLE_ID', 'CANCER_TYPE')],
            by.x = 'sample', by.y = 'SAMPLE_ID', all.x = T)

RenoCancer = loy[which(loy$CANCER_TYPE == 'Renal Cell Carcinoma'), ]
RenoCancer$Y_call[which(RenoCancer$classification %in% c('loss', 'relative_loss'))] = 'Y_chrom_loss'
RenoCancer$Y_call[which(RenoCancer$classification %in% c('gain', 'wt'))] = 'intact_Y_chrom'
RenoCancer$Y_expected = NULL
RenoCancer$classification = NULL
RenoLoss = unique(RenoCancer$sample[which(RenoCancer$Y_call == 'Y_chrom_loss')])

onco_cna = read.csv('Data/signedOut/data_CNA.oncokb.txt.gz', sep = '\t')
onco_cna_reno = onco_cna[which(onco_cna$SAMPLE_ID %in% RenoCancer$sample), ]

onco_mut = read.csv('Data/signedOut/data_mutations_extended.oncokb.txt', sep = '\t')
onco_mut_reno = onco_mut[which(onco_mut$Tumor_Sample_Barcode %in% RenoCancer$sample), ]

##----------------+
## Assign Status to Cohort
##----------------+
RenoCancer$TP53 = NA
RenoCancer$KDM5C = NA
RenoCancer$VHL = NA
RenoCancer$MTOR = NA

for(i in 1:nrow(RenoCancer)){
  for(j in c('TP53', 'KDM5C', 'VHL', 'MTOR')){
    if(j %in% onco_mut_reno$Hugo_Symbol[which(onco_mut_reno$Tumor_Sample_Barcode == RenoCancer$sample[i])]){
      RenoCancer[i, j] = 1
    } else {
      RenoCancer[i, j] = 0
    }
  }
}

#' Adding CNA info
RenoCancer[which(RenoCancer$sample == 'P-0059086-T01-IM7'), 'KDM5C'] = 1
RenoCancer[which(RenoCancer$sample == 'P-0019971-T01-IM6'), 'TP53'] = 1
RenoCancer[which(RenoCancer$sample == 'P-0048377-T01-IM6'), 'VHL'] = 1
RenoCancer[which(RenoCancer$sample == 'P-0012516-T01-IM5'), 'VHL'] = 1

#' convert data to long format
RenoCancer$CANCER_TYPE = NULL
RenoCancer$Y_call[which(RenoCancer$Y_call == 'Y_chrom_loss')] = 1
RenoCancer$Y_call[which(RenoCancer$Y_call == 'intact_Y_chrom')] = 0

#RenoCancer = gather(RenoCancer, gene, value, Y_call:MTOR, factor_key = T)
#RenoCancer = RenoCancer[which(RenoCancer$sample %in% RenoLoss), ]

## test
x = RenoCancer
x$ploidy = NULL
x = t(x)
colnames(x) = x[1, ]
x = x[-1, ]
x = as.data.frame.matrix(x, row.names = row.names(x))
y = sapply(x, as.numeric)
row.names(y) = row.names(x)

z = y[-4, ]


z = memoSort(M = z)


oncoPrint(z)
oncoPrint <- function(M, sort=TRUE) {
  if(sort) {
    alts <- memoSort(M);		
  } else {
    alts <- M;
  }
  
  ngenes <- nrow(alts);
  nsamples <- ncol(alts);
  coverage <- sum(rowSums(alts) > 0);
  
  ### OncoPrint
  numOfOncos <- ngenes*nsamples;
  oncoCords <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos );
  colnames(oncoCords) <- c("xleft", "ybottom", "xright", "ytop", "altered");
  
  xpadding <- .01;
  ypadding <- .01;
  cnt <- 1;
  for(i in 1:ngenes) {
    for(j in 1:nsamples) {
      xleft <- j-1 + xpadding;
      ybottom <- ((ngenes-i+1) -1) + ypadding;
      xright <- j - xpadding;
      ytop <- (ngenes-i+1) -ypadding;
      altered <- alts[i, j];
      
      oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
      cnt <- cnt+1;
    }
  }
  
  colors <- rep("white", cnt);
  colors[ which(oncoCords[, "altered"] == 1) ] <- "black";
  plot(c(0, nsamples), c(0, ngenes), type="n", main=sprintf("Gene set altered in %.2f%%: %d of %d cases", coverage/nsamples*100, coverage, nsamples), xlab="Samples", ylab="", yaxt="n");
  rect(oncoCords[, "xleft"], oncoCords[, "ybottom"],oncoCords[, "xright"], oncoCords[, "ytop"], col=colors, border="white");
  axis(2, at=(ngenes:1)-.5, labels=rownames(alts), las=2);
}

oncoPrint(y)

View(y)


#' make a plot - co-occurence; mutual exclusive
ggplot(RenoCancer, aes(x = sample, y = gene, fill = value)) +
  geom_tile() +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c('1' = 'red',
                                '0' = 'grey85'))
  
  

