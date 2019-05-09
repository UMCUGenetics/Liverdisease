# @Date: March 2018
# @Author: Myrthe Jager
# @Modified: June 2018, September 2018
# @Description: Plot VAFs and TAFs of HCC biopsies




# ---- GET STARTED ----

#1 Load required packages
library(ggplot2)
library("gridExtra")

#2 Define input and output directory
dir = ""
indir = paste(dir,"Data/",sep = "")
outdir = paste(dir,"Results/",sep = "")

#3 Functions




# ---- Plot VAFs & TAFs ----

#1 Read table
vaf_df <- read.delim(paste(indir,"SNV/HCC/clonal_withinfo_forVAF/HCC_clonal_SNV_autosomal_filtered_withinfo.vcf", sep=""))
#1B Exclude chr 1 and 8, as these frequently carry copy number alterations
vaf_df <- vaf_df[which(vaf_df$X.CHROM != 1 & vaf_df$X.CHROM != 8),]

#2 Get purity estimates
purity_df <- read.table(paste(indir,"Purity/Purity_HCCsamples.txt",sep=""),header=TRUE)

#3 Generate plots
coreplot <- ggplot(vaf_df, aes(x=VAF_core)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "Core"),]$Purity.HMF.))/2),color = "HMF"),size=1) +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "Core"),]$Purity.pathology.))/2),color = "Pathology"),size=1) +
  ggtitle("core") +
  theme(legend.position = "none")

endplot <- ggplot(vaf_df, aes(x=VAF_end)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "End2"),]$Purity.HMF.))/2),color = "HMF"),size=1) +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "End2"),]$Purity.pathology.))/2),color = "Pathology"),size=1) +
  ggtitle("end")+
  theme(legend.position = "none")

sbplot <- ggplot(vaf_df, aes(x=VAF_sb)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "SampleB"),]$Purity.HMF.))/2),color = "HMF"),size=1) +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "SampleB"),]$Purity.pathology.))/2),color = "Pathology"),size=1) +
  ggtitle("sample B")+
  theme(legend.position = "none")

scplot <- ggplot(vaf_df, aes(x=VAF_sc)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "SampleC"),]$Purity.HMF.))/2),color = "HMF"),size=1) +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "SampleC"),]$Purity.pathology.))/2),color = "Pathology"),size=1) +
  ggtitle("sample C")+
  theme(legend.position = "none")

seplot <- ggplot(vaf_df, aes(x=VAF_se)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "SampleE"),]$Purity.HMF.))/2),color = "HMF"),size=1) +
  geom_vline(aes(xintercept=(as.numeric(as.character(purity_df[which(purity_df$Sample == "SampleE"),]$Purity.pathology.))/2),color = "Pathology"),size=1) +
  ggtitle("sample E")

#4 Plot
grid.arrange(coreplot,endplot,sbplot,scplot,seplot)
#ggsave(grid.arrange(coreplot,endplot,sbplot,scplot,seplot),filename = paste(outdir,"SNV/HCC/VAF_clonalevents.pdf",sep=""))

#5 Calculate adjusted VAFs (for tumor percentage) per biopsy
vaf_df$TAF_core <- vaf_df$VAF_core / purity_df[which(purity_df$Sample == "Core"),]$Purity.pathology.
vaf_df$TAF_end <- vaf_df$VAF_end / purity_df[which(purity_df$Sample == "End2"),]$Purity.pathology.
vaf_df$TAF_sb <- vaf_df$VAF_sb / purity_df[which(purity_df$Sample == "SampleB"),]$Purity.pathology.
vaf_df$TAF_sc <- vaf_df$VAF_sc / purity_df[which(purity_df$Sample == "SampleC"),]$Purity.pathology.
vaf_df$TAF_se <- vaf_df$VAF_se / purity_df[which(purity_df$Sample == "SampleE"),]$Purity.pathology.

#6 Generate plots of TAFs
coreplot2 <- ggplot(vaf_df, aes(x=TAF_core)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=0.5), color = "red") +
  ggtitle("Core") +
  theme(legend.position = "none")

endplot2 <- ggplot(vaf_df, aes(x=TAF_end)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=0.5), color = "red") +
  ggtitle("End")+
  theme(legend.position = "none")

sbplot2 <- ggplot(vaf_df, aes(x=TAF_sb)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=0.5), color = "red") +
  ggtitle("Sample B")+
  theme(legend.position = "none")

scplot2 <- ggplot(vaf_df, aes(x=TAF_sc)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=0.5), color = "red") +
  ggtitle("Sample C")+
  theme(legend.position = "none")

seplot2 <- ggplot(vaf_df, aes(x=TAF_se)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_density(alpha=.2, fill="#FF6666") + 
  scale_x_continuous(limits=c(0,1)) +
  labs(x = "Fraction non-reference (PNR)") +
  geom_vline(aes(xintercept=0.5), color = "red") +
  ggtitle("Sample A")

#7 Plot
grid.arrange(coreplot2,endplot2,seplot2,sbplot2,scplot2)
#ggsave(grid.arrange(coreplot2,endplot2,seplot2,sbplot2,scplot2),filename = paste(outdir,"SNV/HCC/TAF_clonalevents.pdf",sep=""))