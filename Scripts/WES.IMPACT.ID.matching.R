## Match WES recapture sequencing ID's with normal P-XXX-XXX's listed in cbio
## data obtained through: https://delphi.mskcc.org/sample-tracker/
## Access to WES samples, via special permission (lab manager; Jene White)
## data is basically stored @ /juno/work/tempo/wes_repo/Results/v1.3.x/cohort_level/MSKWESRP/somatic
 
## data fetched: 24/11/2020
## script modified: 15/03/2021

rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/04_Y_chromo_loss/')


## Libraries
library(readxl)

## Input
Delphi.IDs = readxl::read_excel(path = '~/Documents/MSKCC/04_Y_chromo_loss/Delphi_IDs.xlsx')
Delphi.IDs$`CMO Sample ID` = gsub(pattern = '-', replacement = '_', x = Delphi.IDs$`CMO Sample ID`)
WES.Pathnames = read.csv('~/Documents/MSKCC/04_Y_chromo_loss/WES_sampleIDs.txt', sep = '\t', header = F)


# function to link WES and cBIO IDs
Delphi.IDs$Facet_Path = NA
for(i in 1:nrow(Delphi.IDs)){
  Delphi.IDs$Facet_Path[i] = grep(pattern = Delphi.IDs$`CMO Sample ID`[i], 
                                  x = WES.Pathnames$V1, value = T)
}

# data rearrangements:
# data handling on cluster (later)
Delphi.IDs = Delphi.IDs[,c(1, 2, 3, 5, 6, 7, 8, 9, 17)]
Delphi.IDs$Facet_Path = paste0(Delphi.IDs$Facet_Path, '/facets/', Delphi.IDs$Facet_Path,
                               '/', Delphi.IDs$Facet_Path, '.snp_pileup.gz')
Delphi.IDs$Facet_Absolute_Path = paste0('/juno/work/tempo/wes_repo/Results/v1.3.x/cohort_level/MSKWESRP/somatic/', 
                                        Delphi.IDs$Facet_Path)
colnames(Delphi.IDs)[10] = 'Facet_Countfile'
colnames(Delphi.IDs) = gsub(pattern = ' ', replacement = '_', x = colnames(Delphi.IDs))

write.table(Delphi.IDs, 
            file = '~/Documents/MSKCC/04_Y_chromo_loss/Data/WES_cBio_ID.match.txt', 
            sep = '\t', 
            quote = F, 
            row.names = F)


## Visualize the distribution of WES samples (as per 11/2020)
WES.IMPACT.ID.match_df = read.csv('~/Documents/MSKCC/04_Y_chromo_loss/Data/WES_cBio_ID.match.txt', sep = '\t')
Freq.WES = table(WES.IMPACT.ID.match_df$Parental_Tumor_Type)
Freq.WES = Freq.WES[order(Freq.WES, decreasing = T)]
Freq.WES = Freq.WES[1:15]


pdf(file = 'Figures/WES.cancertype.distribution.pdf', width = 9, height = 6)
par(oma = c(0,10,0,0))
barplot(Freq.WES[15:1], horiz = T, border = NA, axes = F, las = 2, xlim = c(0, 420))
axis(side = 1, at = c(0, 100, 200, 300, 400), labels = c(0, 100, 200, 300, 400))
dev.off()

## out
