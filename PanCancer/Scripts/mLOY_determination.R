##----------------+
## MADSEQ; downstream
## analysis; identification
## of mLOY
##----------------+

## start: 12/09/2022
## revision: 12/12/2022
## chris-kreitzer

clean()
gc()
.rs.restartR()

setwd('~/Documents/MSKCC/10_MasterThesis/')


##----------------+
## which samples where
## not considered for
## mLOY calculation
##----------------+
IMPACT = read.csv('Data/00_CohortData/IMPACT_dataFreeze_07.13.22.txt', sep = '\t')
IMPACT_mLOY = IMPACT[!is.na(IMPACT$Normal_bam), ]
mLOY_files = list.files('Data/03_Mosaicism/Normalized/', full.names = T)
mLOY_files_short = substr(x = basename(mLOY_files), start = 1, stop = 17)
length(intersect(mLOY_files_short, IMPACT_mLOY$SAMPLE_ID))
length(intersect(mLOY_files_short, IMPACT_mLOY$SAMPLE_ID)) == length(mLOY_files)

## TODO: number may be different
## because lymphoid and myeloid cancers
## were excluded; old masterfile used 
## for mLOY calculation



##----------------+
## read the data and 
## calculate read-depth 
## ratio
##----------------+
mLOY_files = list.files('Data/03_Mosaicism/Normalized/', full.names = T)
read_mLOY = function(data){
  print(data)
  name = substr(basename(data), start = 1, stop = 17)
  file = read.csv(data, sep = '\t')
  file$width = NULL
  file$strand = NULL
  file$id = name
  file
}

x = lapply(unique(mLOY_files), function(x) read_mLOY(data = x))
x = data.table::rbindlist(x)
xx = as.data.frame(x)

##----------------+
## Downstream anlysis
##----------------+
loy_ratio = function(id, data){
  name = id
  data = data[which(data$id == id), ]
  print(name)
  Y = median(data$normed_depth[which(data$seqnames == 'chrY')])
  A = median(data$normed_depth[!data$seqnames %in% c('chrX', 'chrY')])
  c22 = median(data$normed_depth[which(data$seqnames == 'chr22')])
  c22_ratio = c22 / A
  ratio = Y / A
  ratio_norm = Y / unique(data$ref_depth)
  out = data.frame(id = name,
                   ratio = ratio,
                   ratio_norm = ratio_norm,
                   c22 = c22_ratio)
  out
}

y = lapply(unique(xx$id), function(x) loy_ratio(id = x, data = xx))
y = data.table::rbindlist(y)


##----------------+
## expected ploidies:
## correction factor
##----------------+
par(mfrow = c(1,2))
plot(density(y$c22),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5)
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'DNA concentration [target chr./full genome]', side = 1, line = 2)
mtext(text = 'Chromosome 22', side = 3, line = 0.5, adj = 0)

plot(density(y$ratio),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5)
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'DNA concentration [target chr./full genome]', side = 1, line = 2)
mtext(text = 'Chromosome Y', side = 3, line = 0.5, adj = 0)
# save: (device size): 8.61 x 3.82

y$E22 = y$c22 * 2
density.chromo22 = density(y$E22, na.rm = T)
density.max22 = density.chromo22$x[which.max(density.chromo22$y)]
correction_factor22 = density.max22 - 2

y$EY = y$ratio * 2
density.chromoY = density(y$EY, na.rm = T)
density.maxY = density.chromoY$x[which.max(density.chromoY$y)]
correction_factorY = density.maxY - 1

y$O22 = abs(y$E22 - correction_factor22)
y$OY = abs(y$EY - correction_factorY)


## plot: observed vs corrected
par(mfrow = c(2,2))
plot(density(y$c22),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5)
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'DNA concentration [target chr./full genome]', side = 1, line = 2)
mtext(text = 'Chromosome 22: observed', side = 3, line = 0.5, adj = 0)

plot(density(y$ratio),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5)
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'DNA concentration [target chr./full genome]', side = 1, line = 2)
mtext(text = 'Chromosome Y: observed', side = 3, line = 0.5, adj = 0)

plot(density(y$O22),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5)
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'Ploidy', side = 1, line = 2)
mtext(text = 'Chromosome 22: corrected', side = 3, line = 0.5, adj = 0)

plot(density(y$OY),
     yaxt = 'n',
     ylab = '',
     xlab = '',
     main = '',
     lwd = 1.5)
box(lwd = 2)
mtext(text = 'Density', side = 2, line = 1)
mtext(text = 'Ploidy', side = 1, line = 2)
mtext(text = 'Chromosome Y: corrected', side = 3, line = 0.5, adj = 0)

# save: (device size): 8.61 x 5


##----------------+
## determine 99% CI cutoff
##----------------+
y$seq = seq(1, nrow(y), 1)

jitter = ggplot(y, aes(x = OY, y = seq)) +
  geom_jitter(shape = 17, size = 0.7) +
  geom_vline(xintercept = quantile(y$OY, probs = 0.01)[[1]], col = 'red', linewidth = 0.75) +
  geom_vline(xintercept = quantile(y$OY, probs = 0.99)[[1]], col = 'red', linewidth = 0.75) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10, colour = 'black'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(y = 'Individuals', x = 'Ploidy')

histo = ggplot(y, aes(x = OY, y = ..density..)) +
  geom_histogram(bins = 400, col = 'black', fill = 'black') +
  geom_vline(xintercept = quantile(y$OY, probs = 0.01)[[1]], col = 'red', linewidth = 0.75) +
  geom_vline(xintercept = quantile(y$OY, probs = 0.99)[[1]], col = 'red', linewidth = 0.75) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 1.5),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

plot_grid(histo, jitter, nrow = 2, rel_heights = c(1,5), align = 'hv')
# save: (custom device)


##----------------+
## cut-off threshold
##----------------+
lower_cutoff = quantile(y$OY, probs = 0.01)[[1]]
upper_cutoff = quantile(y$OY, probs = 0.99)[[1]]
