# Y-chromosome-loss

## Methological approach for the definition of a Y-chromosome loss  
Both, whole-exome and IMPACT panel seq. data was analyzed with FacetsY. FacetsY, is slight modification of the original Facets-pipeline, which considers and include seq. reads coming from the Y-chromosome (big credits to Bastien, for re-writing this script). 
FacetsY was tested on WES TCGA-PRAD (n=333, Cell 2015) cohort (mostly primary samples) to see the overall segmentation agreement rate.   
FacetsY was run with the following parameters:   
- cval.preprocess = 25   
- cval.postprocess = 150   
- snp.nbhd = 250

Subsequently, segmented data was compared to Affymetrix 6.0 data derived from TCGA-PRAD (broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg); see corresponding script.
The overall segmentation pattern was assessed visually via IGV as well as the arm-level alteration frequency in respective cohorts.
Moreover, we called Y-chromosome losses in TCGA-PRAD samples (5.7%) and clearly saw, that some genes on the Y-chromosome show a significant reduction of expression compared to its Y-chromosome (wild-type) counterparts).


### Concentrate on sequencing reads from the Y-chromsome (see Ed Rezniks paper)

### Sanity checks:
A. FacetsY (default parameters) work fine on TCGA data (comparison with Affymetrix Chips)
B. WES capture most of the feature across the Y-chromosome. Mostly protein coding domains.

#### IMPACT Panel Sequencing (FacetsY):
With Facets default settings (ndepth = 35, ndepthmax = 1000) the coverage along the Y chromosome is enourmously decreased. Roughly 3% or 380 positions (out of ~11,000) retain for segmentation. These ~380 positions (breakpoints) cluster in around 3 bins (see supplementary figure 1) whereof (mostly) one segment is called.
These segments (tcn.em) are assigned a copy number state: Either 0 (meaning its lost) or >0, meaning there is one Y-chromosome available or even more copies (likely due to Whole-Genome-duplication)

#### WES recapture Sequencing on PRAD samples (n = 238)
We obtain a conitinuous coverage across the Y-chromosome. Most protein coding genes are directly covered, and hence we believe that WES provide us with a good approximation for copy number calls. 
However, just looking into a binary state of copy-numbers (tcn.em) being either 0 or 1, we loose a lot of information and hence IMPACT to WES analysis have a high disagreement rate. 

#### Gene expression analysis of Y-chromosome:
GTEx portal provides expression data on numerous genes throughtout different tissue types. Importantly, all data reported on GTEx are derived from healthy donors, meaning that unbiased, non-diseased tissue was analysed. In this analysis, I searched for 'features' (genes, pseudogenes, rRNAs, etc.) which are ubiquitiously expressed throughout the body (tissues). I restricted the analysis to i) male samples and ii) those genes with a median expression (log-scaled TPM) != 0. 
See Figure-section for the result (Heatmap made via ComplexHeatmap)

#### Orthogonal validation of Y-chromosome loss:
The LogR-Ratio [log2(normal/tumor) + 1] at individual postions along the Y-chromosome, serves as surrogate for FacetsY estimations. The median LRR (CnLR) is broadley utilized in several publications and hence seems to provide a stable prediction. 
Based on the Seq. Coverage of IMPACT-Panel samples, we observed that informative SNPs for segmentation are unevenly distributed and hence may obscure a proper Y-chromosome loss. With the our example (Figure folder) we show that, using the median CnLR of indivual genes we see a clear devitation between the two groups. 


### Data
WES_cBio_ID.match = SampleID's WES which match with cBIO ID's   
Prostate_Y_loss = Prostate cancer samples where Y chromosome loss is confirmed in both, WES and IMPACT Panel sequencing respectively   
Y_chromosome_genes = Gene coordinates obtained from Biomart (ENSEMBL). The hg19 (archived) reference genome was used, as FacetsY works on hg19

### Scripts   
FacetsY.wes.prostate.submission = Script which calls FacetsY (locally on juno) and run on n = 238 prostate samples which were WES(ed)

### Results: Y-chromosome loss:
A. Pan-cancer whole-genome analyses of metastatic solid tumours (Hartwig).
- 17% Y-chromosome loss (Prostate)(n=212)[IMPORTANT: WORKING on METASTASES]
- con1: disagreement between PCAWG and Hartwig (Indels, SVs)
- con2: using 'own' platform (purple)
- con3: using lower median depth (Seq. ~ 100, ~ 38N)
- pro1: WGS!

### Refinements:
A. Transform binary call (50% rule tcn.em == 0) into continious scale. Plot WES vs IMPACT samples and see whether discrepancies can be resolved (cluster of certain samples, etc.)(updated script).


### Mosaicism in IMPACT & WES (correspondence with Peter Priestly)
#### new chapter on the way!   
# 22,000 cases!
# new chapter 

