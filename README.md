# Liver disease README

## Date
July 2019

## Authors
Myrthe Jager (contact: m.jager-2@umcutrecht.nl)
Bastiaan van der Roest (contact: B.R.vanderRoest-2@umcutrecht.nl)

## Coauthors
Francis Blokzijl
Roel Janssen
Ruben van Boxtel

## Description
Here, you can find all code that can be used to reproduce the results of "Mutational impact of chronic alcohol use on stem cells in cirrhotic liver" (https://doi.org/10.1101/698894).

## Data
The whole-genome sequencing and RNA sequencing data generated during the current study are available at EGA (https://www.ebi.ac.uk/ega/home) under accession number EGAS00001002983. We also used data produced in Blokzijl et al., Nature 2016 (https://www.nature.com/articles/nature19768).
Filtered VCF-files, metadata, BED-files with callable regions, and RNA-Seq counts generated during the current study are available at Zenodo under DOI 10.5281/zenodo.3295513 (https://doi.org/10.5281/zenodo.3295513). 


## Instructions for use
### A Start up (estimated time < 1 hr)
1. Download and install required packages (see below @ Software dependencies & R packages)
2. Download the data from Zenodo (https://doi.org/10.5281/zenodo.3295513)
### B Mutational patterns (estimated time ~ 1.5 hr); by Myrthe Jager
3. Create output directory in downloaded Zenodo directory '/Results/' with subfolders:
- /SNV/HCC/
- /SNV/coding/
- /SNV/IAP_vs_HMF/
- /INDEL/
- /RNA/
4. For each script in steps 6-12: fill 'dir' with directory of data downloaded from zenodo (once per script).
5. Fill 'beddir' in script '190507_Liver_MutPat.R' with 'dir/Data/callable/beddir/' once
6. Run '190507_LiverMutPat.R' to get & plot (~ 30 mins runtime): 
- (tandem) base substitution numbers
- 6-type mutation spectra
- 96-type mutational profiles
- mutational profiles reconstructed with COSMIC signatures v2 and v3
- rainfall plots
7. Run '190507_bootstrap.Rmd' to calculate whether there are significant differences between signature contributions in Healthy vs Alcoholic liver (~ 25 mins runtime)
8. Run '190507_LiverdNdS.R' for dNdS and non synonymous mutations in known cancer-driver genes. (~ 15 mins runtime)
9. Run '180420_HMFvsIAP.R' for a comparison between the two calling pipelines used in this study (~ 5 mins runtime)
-https://github.com/UMCUGenetics/IAP
-https://github.com/hartwigmedical
10. Run '180914_clonalVAF.R' to get & plot the TAF of the MRCA mutations in all HCC biopsies (~ 1 min runtime)
11. Run '181009_Liver_INDEL.R' to get & plot INDEL numbers (~ 3 mins runtime)
12. Run '190507_DE.R' for RNA seq analysis and expression analysis of PTPRK (~ 3 mins runtime)
### C H3K36me3; By Bastiaan van der Roest (estimated time ~ 1 day)
13. Download H3K36me3 ChIP-Seq data from UCSC (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/) into directory ‘dir’
14. Run ‘combineSignals.sh’ to get one bedGraph file containing all samples
15. Run ‘peakCalling.sh’ to get the H3K36me3 peaks
16. Import H3K36me3 into your SPARQL graph database
17. Run ‘count_overlaps.R’ to count SNVs in H3K36me3 peaks

## Software dependencies
R version 3.5.1
bedtools 2.27.1
macs2 2.1.1.20160309
graphdb-hpc-tools 0.0.2
bigWigToBedGraph 377-1

## R packages
AnnotationDbi	1.44.0
AnnotationHub	2.14.2
Biobase	2.42.0
BiocGenerics	0.28.0
BiocInstaller	1.30.0
BiocParallel	1.16.5
biomaRt	2.36.1
Biostrings	2.50.2
Bsgenome	1.50.0
BSgenome.Hsapiens.UCSC.hg19	1.4.0
cluster	2.0.7-1
DelayedArray	0.8.0
DESeq2	1.20.0
devtools	2.0.2
dndscv	0.0.1.0
doParallel	1.0.14
dplyr	0.7.8
foreach	1.4.4
futile.logger	1.4.3
GenomeInfoDb	1.18.1
GenomicFeatures	1.34.1
GenomicRanges	1.34.0
ggplot	3.1.0
ggplot2	3.1.1
ggsignif	0.4.0
gridExtra	2.3
Iranges	2.16.0
MASS	7.3-50
matrixStats	0.54.0
MutationalPatterns	1.6.2
nlme	3.1-137
NMF	0.21.0
org.Hs.eg.db	3.6.0
pheatmap	1.0.12
pkgmaker	0.27
RColorBrewer	1.1-2
registry	0.5-1
reshape2	1.4.3
rngtools	1.2.4
Rsamtools	1.34.0
rtracklayer	1.42.1
S4Vectors	0.20.1
seqinr	3.4-5
sparql	1.16
SummarizedExperiment	1.12.0
TxDb.Hsapiens.UCSC.hg19.knownGene	3.2.2
usethis	1.5.0
VariantAnnotation	1.28.10
VennDiagram	1.6.20
Xvector	0.22.0
