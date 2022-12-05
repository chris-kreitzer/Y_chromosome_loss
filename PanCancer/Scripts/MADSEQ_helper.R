##----------------+
## MADSEQ helper functions;
##----------------+

library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)
library(BSgenome)
library(diffloop)
library(BSgenome.Hsapiens.UCSC.hg19)


## given a GRange object and bam file, 
## calculate average coverage for each range
calculateSubCoverage = function(range, bam){
  ## read bam file from given ranges,
  ## filter out duplicated reads, secondary reads and unmapped reads
  ## exclude reads with mapQ==0
  param = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, 
                                        isSecondaryAlignment=FALSE, 
                                        isDuplicate=FALSE),
                       which=range,
                       mapqFilter=1)
  ## read alignment
  sub_alignment = GenomicAlignments::readGAlignments(bam, param = param)
  ## calculate coverage
  cov = GenomicAlignments::coverage(sub_alignment)
  cov = cov[range]
  ## return average coverage for each region
  round(mean(cov))
}


## given the path to targeted bed file and bam file, 
## create a GRanges object containing coverage for each targeted region
getCoverage = function(bam, target_bed, genome_assembly="hg19"){
  
  ## read in target bed table
  target_gr = rtracklayer::import(target_bed)
  if(nchar(seqlevels(target_gr)[1])>3){
    seqlevels(target_gr,pruning.mode="coarse")=c("chr1","chr2","chr3","chr4","chr5",
                                                 "chr6","chr7","chr8","chr9","chr10",
                                                 "chr11","chr12","chr13","chr14",
                                                 "chr15","chr16","chr17","chr18",
                                                 "chr19","chr20","chr21","chr22",
                                                 "chrX","chrY")
  }
  else{
    seqlevels(target_gr,pruning.mode="coarse")=c("1","2","3","4","5",
                                                 "6","7","8","9","10",
                                                 "11","12","13","14",
                                                 "15","16","17","18",
                                                 "19","20","21","22",
                                                 "X","Y")
  }
  target_gr = sort(target_gr)
  # remove regions overlapped with REs
  #target_gr = removeRE(target_gr,genome_assembly)
  target_gr = removeGap(target_gr,genome_assembly)
  target_gr = removeHLA(target_gr,genome_assembly)
  target_gr = removeAQP(target_gr,genome_assembly)
  nRegion = length(target_gr)
  cat(paste(nRegion, "non-repeats regions from", length(seqlevels(target_gr)),
            "chromosomes in the bed file.", sep=" "))
  ## use helper function to calculate average coverage for each region, 
  ## in order to handle large bam files, 
  ## process 1000 regions at a time to reduce memory usage
  message("calculating depth from BAM...")
  depth = NULL
  for (i in seq(1,nRegion,1000)){
    ## report progress
    if(i%%5000==1&i>1) cat(paste(i-1,"regions processed\n"))
    end = ifelse(i+999>nRegion,nRegion,i+999)
    sub_depth = calculateSubCoverage(target_gr[i:end],bam)
    depth = c(depth, sub_depth)
  }
  ## see if number of depth equals to number of regions
  if (length(depth) == nRegion) mcols(target_gr)$depth = depth
  else stop(paste("with", nRegion, "target regions, only", 
                  length(depth), "processed. Please check your input.",
                  sep=" "))
  target_gr
}


## helper function to prepare gc data
## given targets as a GRanges object, calculate gc content for each range
calculateGC = function(
    range,
    genome_assembly="hg19"){
  genome = BSgenome::getBS(genome_assembly)
  ## use alphabetFrequency function in biostring to calculate GC percent
  message("calculating GC content...")
  base_frequency = alphabetFrequency(BSgenomeViews(genome,range),
                                     as.prob = TRUE)[,c("C","G")]
  gc_content = apply(base_frequency,1,sum)
  gc_content
}

## helper function to correct coverage by GC content
## given a GRanges object, output corrected coverage
correctGCBias = function(
    object,
    plot=FALSE){
  ## correct coverage by GC content
  ## convert GRanges object to frame
  gc_depth = as.data.frame(object)
  name = names(gc_depth)
  
  ## check if data has been quantile normalized
  ## if data has been quantile normalized, following analysis is operated
  ## on quantiled depth
  quantiled = FALSE
  if (is.element("quantiled_depth",name)){
    ## exclude regions with 0 coverage
    quantiled = TRUE
    gc_depth = gc_depth[gc_depth$depth>0&gc_depth$quantiled_depth>0,]
    names(gc_depth) = sub("quantiled_depth","coverage",name)
  }
  else{
    gc_depth = gc_depth[gc_depth$depth>0,]
    names(gc_depth) = sub("depth","coverage",name)
  }
  
  ## round gc content to 0.001 increments
  gc_depth = data.frame(gc_depth,round_gc=round(gc_depth$GC,3))
  ## split data by GC content
  split_gc = split(gc_depth,gc_depth$round_gc)
  coverage_by_gc = sapply(split_gc,function(x)mean(x$coverage,na.rm=TRUE))
  gc_coverage = data.frame(round_gc = as.numeric(names(coverage_by_gc)),
                           mean_reads = coverage_by_gc)
  ## fit coverage and GC content by loess
  gc_coverage_fit = stats::loess(gc_coverage$mean_reads~gc_coverage$round_gc,
                                 span=0.5)
  ## the expected coverage is the mean of the raw coverage
  expected_coverage = mean(gc_depth[,"coverage"])
  
  ## plot GC vs. raw reads plot
  if(plot == TRUE){
    plot(x = gc_coverage$round_gc, 
         y = gc_coverage$mean_reads, 
         pch = 16, col = "blue", cex = 0.6, 
         ylim = c(0,1.5*quantile(gc_coverage$mean_reads,0.95,na.rm=TRUE)),
         xlab = "GC content", ylab = "raw reads", 
         main = "GC vs Coverage Before Norm",cex.main=0.8)
    graphics::lines(gc_coverage$round_gc, 
                    stats::predict(gc_coverage_fit, gc_coverage$round_gc), 
                    col = "red", lwd = 2)
    graphics::abline(h = expected_coverage, lwd = 2, col = "grey", lty = 3)
  }
  
  ## correct reads by loess fit
  normed_coverage = NULL
  for (i in 1:24){
    ## check if the coordinate is with "chr" or not
    if(nchar(as.character(seqnames(object)@values[1]))>3) {
      chr = paste("chr",i,sep="")
      if (i == 23) chr = "chrX"
      if (i == 24) chr = "chrY"
    }
    else{
      chr = i
      if (i == 23) chr = "X"
      if (i == 24) chr = "Y"
    }
    tmp_chr = gc_depth[gc_depth$seqnames == chr,]
    if(nrow(tmp_chr)==0) next
    chr_normed = NULL
    for (j in 1:nrow(tmp_chr)){
      tmp_coverage = tmp_chr[j,"coverage"]
      tmp_GC = tmp_chr[j,"GC"]
      # predicted read from the loess fit
      tmp_predicted = stats::predict(gc_coverage_fit, tmp_GC)
      # calculate the error biased from expected
      tmp_error = tmp_predicted - expected_coverage
      tmp_normed = tmp_coverage - tmp_error
      chr_normed = c(chr_normed, tmp_normed)
    }
    normed_coverage = c(normed_coverage,chr_normed)
  }
  gc_depth = cbind(gc_depth,normed_coverage = normed_coverage)
  gc_depth = gc_depth[!is.na(gc_depth$normed_coverage),]
  
  ## calculate and plot GC vs coverage after normalization
  split_gc_after = split(gc_depth,gc_depth$round_gc)
  coverage_by_gc_after = sapply(split_gc_after,
                                function(x)mean(x$normed_coverage,
                                                na.rm=TRUE))
  gc_coverage_after=data.frame(round_gc=
                                 as.numeric(names(coverage_by_gc_after)),
                               mean_reads=coverage_by_gc_after)
  gc_coverage_fit_after = stats::loess(gc_coverage_after$mean_reads
                                       ~gc_coverage_after$round_gc,span=0.5)
  
  ## plot GC vs coverage after normalization
  if (plot == TRUE){
    plot(x = gc_coverage_after$round_gc, 
         y = gc_coverage_after$mean_reads, 
         pch = 16, col = "blue", cex = 0.6, 
         ylim = c(0,1.5*quantile(gc_coverage_after$mean_reads,0.95,
                                 na.rm=TRUE)),
         xlab = "GC content", ylab = "normalized reads", 
         main = "GC vs Coverage After Norm",cex.main=0.8)
    graphics::lines(gc_coverage_after$round_gc, 
                    stats::predict(gc_coverage_fit_after, gc_coverage_after$round_gc),
                    col = "red", lwd = 2)
  }
  
  ## round normalized coverage to integer
  gc_depth$normed_coverage = round(gc_depth$normed_coverage)
  ## exclude regions with corrected coverage <0
  gc_depth = gc_depth[gc_depth$normed_coverage>0,]
  
  ## convert gc_depth into a GRanges object
  if (quantiled == TRUE){
    res = GRanges(seqnames = Rle(gc_depth$seqnames), 
                  ranges = IRanges(start=gc_depth$start,end=gc_depth$end),
                  strand = rep("*",nrow(gc_depth)),
                  depth = gc_depth$depth,
                  quantiled_depth = gc_depth$coverage,
                  GC = gc_depth$GC,
                  normed_depth = gc_depth$normed_coverage)
  }
  else{
    res = GRanges(seqnames = Rle(gc_depth$seqnames), 
                  ranges = IRanges(start=gc_depth$start,end=gc_depth$end),
                  strand = rep("*",nrow(gc_depth)),
                  depth = gc_depth$coverage,
                  GC = gc_depth$GC,
                  normed_depth = gc_depth$normed_coverage)
  }
  res = res[!is.na(mcols(res)$normed_depth)]
  res
}


## function to calculate mean coverage for each chromosome after normalization
## and could plot out the coverage before and after normalization
## input: a GRangesList object
calculateNormedCoverage = function(object, plot = FALSE){
  if(nchar(seqlevels(object)[1])>3){
    chr_name=c("chr1","chr2","chr3","chr4","chr5",
               "chr6","chr7","chr8","chr9","chr10",
               "chr11","chr12","chr13","chr14",
               "chr15","chr16","chr17","chr18",
               "chr19","chr20","chr21","chr22",
               "chrX","chrY")
  }
  else{
    chr_name = c("1","2","3","4","5","6","7","8","9","10","11","12",
                 "13","14","15","16","17","18","19","20","21","22","X","Y")
  }
  nSample = length(object)
  sample_name = names(object)
  split_object = sapply(object,function(x)split(x,seqnames(x)))
  
  ## calculate average coverage for each chromosome after normalization
  after_chr = NULL
  for (i in 1:nSample){
    sub_after_chr = sapply(split_object[[i]],
                           function(x)mean(mcols(x)$normed_depth))
    after_chr = rbind(after_chr,sub_after_chr[chr_name])
  }
  rownames(after_chr) = sample_name
  after_chr = replace(after_chr,is.nan(after_chr),0)
  
  ## if plot requested, then plot 
  if (plot == TRUE){
    graphics::par(mfrow=c(ifelse(nSample>1,3,2),1))
    ## calculate average coverage before normalization
    before_chr = NULL
    quantiled_chr = NULL
    for (i in 1:nSample){
      sub_before_chr = sapply(split_object[[i]],
                              function(x)mean(mcols(x)$depth))
      before_chr = rbind(before_chr,sub_before_chr[chr_name])
      if (nSample>1){
        sub_quantiled_chr = sapply(split_object[[i]],
                                   function(x)
                                     mean(mcols(x)$quantiled_depth))
        quantiled_chr=rbind(quantiled_chr,sub_quantiled_chr[chr_name])
      }
    }
    before_chr = replace(before_chr,is.nan(before_chr),0)
    quantiled_chr = replace(quantiled_chr,is.nan(quantiled_chr),0)
    ## plot
    cols = sample(grDevices::colors(),nSample,replace = TRUE)
    nChr = ncol(after_chr)
    ## 1. plot raw coverage
    plot(1:nChr,rep(1,nChr),type="n",
         ylim=c(0.5*min(before_chr,na.rm=TRUE),
                1.5*max(before_chr,na.rm=TRUE)),
         xlab="chromosome",ylab="average coverage",
         main = "raw data",xaxt="n")
    graphics::axis(1,at=seq(1,nChr),chr_name[1:nChr],las=2)
    for (i in 1:nSample){
      graphics::lines(1:nChr,before_chr[i,],type="b",pch=16,col=cols[i])
    }
    graphics::legend("topright",sample_name,pch=16,col=cols,
                     ncol=ceiling(nSample/5),cex=0.6)
    
    if(nSample>1){
      ## 2. plot quantiled coverage
      plot(1:nChr,rep(1,nChr),type="n",xaxt="n",
           ylim=c(0.5*min(quantiled_chr,na.rm=TRUE),
                  1.5*max(quantiled_chr,na.rm=TRUE)),
           xlab="chromosome",ylab="average coverage",
           main="quantile normalized")
      graphics::axis(1,at=seq(1,nChr),chr_name[1:nChr],las=2)
      for (i in 1:nSample){
        graphics::lines(1:nChr,quantiled_chr[i,],type="b",
                        pch=16,col=cols[i])
      }
      graphics::legend("topright",sample_name,pch=16,col=cols,
                       ncol=ceiling(nSample/5),cex=0.6)
    }
    
    ## 3. plot normed coverage
    plot(1:nChr,rep(1,nChr),type="n",xaxt="n",
         ylim=c(0.5*min(after_chr,na.rm=TRUE),1.5*max(after_chr,na.rm=TRUE)),
         xlab="chromosome",ylab="average coverage",main="GC normalized")
    graphics::axis(1,at=seq(1,nChr),chr_name[1:nChr],las=2)
    for (i in 1:nSample){
      graphics::lines(1:nChr,after_chr[i,],type="b",pch=16,col=cols[i])
    }
    graphics::legend("topright",sample_name,pch=16,col=cols,
                     ncol=ceiling(nSample/5),cex=0.6)
  }
  after_chr
}



# path = system.file("gap","hg19_gap_gr.RDS",package="MADSEQ")
# a = readRDS(path)
# b = a[seqnames(a) != 'chrY']
# y = a[seqnames(a) == 'chrY']
# y = y[-c(2,3,4,5,6,7)]
# y = y[-c(2,3,4,5,6,7,8)]
# y = y[-2]
# uu = c(b, y)

removeGap = function(gr, genome){
  
    gap_gr = readRDS('~/MADSEQ/hg19_gap_gr.RDS')
    ov = findOverlaps(gr, gap_gr)
    
    if(length(ov)==0){
      res_degap = gr
    } else 
      res_degap = gr[-queryHits(ov)]
    res_degap
}

## remove SNPs inside the HLA region on chr6
# because of the variability of HLA regions, 
# variant calling for this region is problematic most of the time, 
# to keep a clean result, we will filter out SNPs called within this region
# padded 1000kb up and downstream of HLA region

removeHLA = function(gr, genome){
  
  HLA_gr = readRDS('~/MADSEQ/hg19_HLA_gr.RDS')
  ov = findOverlaps(gr, HLA_gr)
  
  if(length(ov)==0){
    res_deHLA = gr
  }
  else
    res_deHLA = gr[-queryHits(ov)]
  res_deHLA
}

## remove SNPs inside the AQP7 region on chr9
# because of the variability of AQP7, 
# variant calling for this region is problematic most of the time, 
# to keep a clean result, we will filter out SNPs called within this region
# padded 1000kb up and downstream of AQP7 region
removeAQP = function(gr,genome){
  AQP7_gr = readRDS('~/MADSEQ/hg19_AQP7_gr.RDS')
  ov = findOverlaps(gr, AQP7_gr)
  
  if(length(ov)==0){
    res_deAQP7 = gr
  }
  else
    res_deAQP7 = gr[-queryHits(ov)]
  res_deAQP7
}




normalizeCoverage = function(
    object,...,control=NULL,
    writeToFile=TRUE,
    destination=NULL,
    plot=FALSE){
  if(is.null(control)) {
    cat("no control provided")
    data = GRangesList(object, ...)
    control_name = NULL
    sample_name = as.character(substitute(list(object,...)))[-1]
  }
  else {
    data = GRangesList(control,object,...)
    control_name = as.character(substitute(control))
    cat(paste("control:",control_name))
    sample_name=c(control_name,
                  as.character(substitute(list(object,...)))[-1])
  }
  nSample = length(data)
  names(data) = sample_name
  message(paste("there are",nSample,"samples"))
  
  ## check if all the samples have the same number of targeted region
  if (length(unique(sapply(data,length)))>1){
    cat(elementNROWS(data))
    stop("the number of targeted region is different in your samples, 
            please check your input")
  }
  
  ## check if there is >1 samples
  if (nSample > 1) {
    if(plot == TRUE) par(mfrow=c(2,2))
    ## 1. quantile normalize data if there are more than one samples
    quantiled_data = coverageQuantile(data)
    ## 2. correct GC bias for each sample
    corrected = NULL
    for (i in 1:nSample){
      message(paste("correct GC bias in sample'",sample_name[i],"'...",
                    sep=" "))
      sub_corrected = correctGCBias(quantiled_data[[i]],plot=plot)
      if(is.null(corrected)) corrected = GRangesList(sub_corrected)
      else corrected = c(corrected,GRangesList(sub_corrected))
    }
  }
  else if(nSample==1){
    if(plot == TRUE) par(mfrow=c(1,2))
    ## only correct coverage by GC bias
    message(paste("correct GC bias in sample '",sample_name[1],"' ...",
                  sep=""))
    corrected = GRangesList(correctGCBias(data[[1]],plot=plot))
  }
  names(corrected) = sample_name
  
  ## use normed coverage to generate reference coverage
  after_chr = calculateNormedCoverage(corrected,plot=plot)
  ## check how many sample
  if(nSample>1){
    ## check if normal control is provided
    if(is.null(control_name)){
      ## if there are more than 1 sample and there is no control
      ## take median of normed average coverage for each chromosome as 
      ## reference
      ref_cov = apply(after_chr,2,function(x)median(x,na.rm = TRUE))
    }
    ## if control is provided, 
    ## take average coverage for control as reference
    else ref_cov = after_chr[control_name,]
    res = NULL
    for (i in 1:nSample){
      sub_res = corrected[[i]]
      chr_length = sapply(colnames(after_chr),
                          function(x)length(sub_res[seqnames(sub_res)==x]))
      mcols(sub_res)$ref_depth = rep(ref_cov,chr_length)
      if(is.null(res)) res = GRangesList(sub_res)
      else res = c(res,GRangesList(sub_res))
    }
    names(res) = sample_name
  }
  ## if there is only one sample, take the median coverage for chromosome as
  ## reference
  else if(nSample==1){
    res = corrected[[1]]
    mcols(res)$ref_depth = rep(median(after_chr),length(res))
    res = GRangesList(res)
    names(res) = sample_name
  }
  
  ## if write to file requested, then write the normalized coverage table 
  ## file, otherwise return it as a GRangesList 
  if(writeToFile == TRUE){
    ## check if path to write file is provided,
    ## if not write to current working directory
    if(is.null(destination)) path = "."
    else path = destination
    for (i in 1:length(res)){
      tmp_data = as.data.frame(res[[i]])
      write.table(tmp_data,
                  file=paste(path, "/" , sample_name[i],
                             "_normed_depth.txt", sep=""),
                  quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")
      message(paste("normalized depth for sample ",sample_name[i],
                    " is written to ",path, "/" , sample_name[i],
                    "_normed_depth.txt",sep=""))
    }
  }
  else res
}





