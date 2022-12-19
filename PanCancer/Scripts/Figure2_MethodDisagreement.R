## Method validation; disagreement plot
## 08/21/2022


#' example:
#' P-0019114-T01-IM6

WES_cnLR = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/WES_cnLR_out.txt', sep = '\t')
IMPACT_cnLR = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/Cnlr_out.txt', sep = '\t')
WES.IMPACT.ids_df = read.csv('~/Documents/MSKCC/10_MasterThesis/Data/02_Method_Validation/WES_cBio_ID.match.txt', sep = '\t')


a = IMPACT_cnLR[which(IMPACT_cnLR$ID == 'P-0019114-T01-IM6'), ]
a = a[order(a$maploc, decreasing = F), ]
a$seq = seq(1, nrow(a), 1)
a$label = a$maploc / 1e6

dev.off()
par(mfrow = c(2,1),
    oma = c(1,0,0,0))
pdf(file = '~/Documents/MSKCC/10_MasterThesis/Figures_original/MethodValidation_DisagreementPlot.pdf',
    width = 10, height = 5, paper = 'a4r', onefile = T)
plot(a$seq, a$cnlr,
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlab = '',
     ylim = c(-5, 2))

axis(side = 1, 
     at = a$seq, 
     labels = round(a$maploc / 1e6, 3), 
     las = 2, 
     cex.lab = 0.2)

axis(side = 2, 
     at = seq(-5, 2, 1),
     las = 1)
box(lwd = 2)
mtext(text = 'CnLR', side = 2, line = 1.8)
mtext(text = 'MSK-IMPACT: P-0019114-T01-IM6', side = 3, line = 1, adj = 0)


## whole exome sequenced sample
b = WES_cnLR[which(WES_cnLR$sample_id == 'C_MA7WDN_M001_d'), ]
b = b[order(b$maploc, decreasing = F), ]
b$seq = seq(1, nrow(b), 1)

plot(b$seq, b$cnlr,
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlab = '',
     ylim = c(-5, 2))
axis(side = 1, 
     at = b$seq, 
     labels = round(b$maploc / 1e6, 3), 
     las = 2, 
     cex.lab = 0.2)

axis(side = 2, 
     at = seq(-5, 2, 1),
     las = 1)
box(lwd = 2)
mtext(text = 'CnLR', side = 2, line = 1.8)
mtext(text = 'MSK-WES: C_MA7WDN_M001_d', side = 3, line = 1, adj = 0)
dev.off()




setwd('~/Documents/MSKCC/10_MasterThesis/')
a = read.csv('Data/02_Method_Validation/mLRR_IMPACT_WES_n950.txt', sep = '\t')
dim(a)
cbio(a$DMP_Sample_ID)
