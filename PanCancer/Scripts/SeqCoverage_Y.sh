#!/bin/bash

# create coverages of all chrY bases in the alignments
# http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html

# check if there are Y-reads in female samples (should be 0)
# filter reads (samtools) based on 

#- MAPQ (mapping quality > 40; ensure that only unique reads are recorded)
#- CIGAR (check if soft-clipping is present; similar to Insertions/Deletions; ideally 100M == 100 base-pair matches/mis-matches)
#- FLAG (-F 0x4) -- only report mapped reads (-F 1x4 == unmapped reads); check for additional FLAGS (read is PCR or optical duplicate (0x400))
# https://broadinstitute.github.io/picard/explain-flags.html

# example: 
samtools view -F 0x4 -q 60 -b /Users/chriskreitzer/Documents/MSKCC/dmp-2021/BAMs/P-0014995/JJ958363-T.bam Y > test.bam
samtools index test.bam test.bai  #' necessary in order to load into IGV eg
bedtools genomecov -ibam test.bam -bg  #' count bases (reads) within given interval; output in .bed format



# ~~~~~ SETUP ~~~~~ #
hg19_chrY="/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/chrY.fa"
qsub_logdir="${PWD}/logs"
mkdir -p "$qsub_logdir"
output_dir="${PWD}/output"
mkdir -p "$output_dir"

