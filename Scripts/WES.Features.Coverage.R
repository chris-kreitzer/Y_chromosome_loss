## PanCancer Analysis on Y-chromosome loss:
## 
## Average number of Features covered on the Y-chromsome by WE-Sequencing
## directly covered by WES- and IMPACT-Seq.
## 
## Start: 14/04/2021


## Input
Features.Y_chromosome = read.csv('/juno/home/kreitzec/WES_Prostate/', sep = '\t')
Sample.path = read.csv('/juno/home/kreitzec/WES_Prostate/Panel.prostate.samplepath.txt', header = F)
Sample.path = as.character(Sample.path$V1)
print(Sample.path)

## Processing:
Protein.coding.genes = Features.Y_chromosome[which(Features.Y_chromosome$gene_biotype == 'protein_coding'), 
                                             c('hgnc_symbol', 'start_position', 'end_position', 'external_gene_name')]



## Analysis:
number.protein.coding.genes = nrow(Protein.coding.genes)

cval.preprocess = 25
cval.postprocess = 200
snp.nbhd = 250

Features.covered_out = data.frame()

for(i in unique(Sample.path)){
  try({
    print(i)
    ID = substr(i, start = 61, stop = 90)
    data.in = facetsY::readSnpMatrix(i)
    data.pre = facetsY::preProcSample(data.in, 
                                      cval = cval.preprocess, 
                                      gbuild = 'hg19', 
                                      snp.nbhd = snp.nbhd)
    data.process = facetsY::procSample(data.pre, 
                                       cval = cval.postprocess)
    
    data.analysis = data.process$jointseg
    data.analysis = data.analysis[which(data.analysis$chrom == 24), ]
    
    #' look into Feature coverage
    for(i in 1:nrow(Protein.coding.genes)){
      if(any(between(x = data.analysis$maploc, lower = Protein.coding.genes$start_position[i], 
                     upper = Protein.coding.genes$end_position[i]))){
        out = data.frame(Gene = Protein.coding.genes$external_gene_name[i],
                         Match = 1,
                         number = sum(between(x = data.analysis$maploc, 
                                              lower = Protein.coding.genes$start_position[i], 
                                              upper = Protein.coding.genes$end_position[i])),
                         ID = ID)
        
      } else next
      Features.covered_out = rbind(Features.covered_out, out)
    }
    
  })
  
}

write.table(Features.covered_out, file = '/juno/home/kreitzec/WES_Prostate/IMPACT.Features.covered.out.txt', row.names = F, sep = '\t')

