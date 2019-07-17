# @Date: 05 april 2018
# @Author: Myrthe Jager
# @Modified: 09 october 2018
# @Description: Genome-wide indel analysis healthy & alcoholic ASCs
# Abbreviations: 




# ---- GET STARTED ----

#1 Load required packages
library(devtools)
library(BiocInstaller)
library(MutationalPatterns)
library("gridExtra")
library("ggplot2")
library(grDevices)
library(reshape2)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(biomaRt)
library(nlme)

#2 Define input and output directory
dir = ""
indir = paste(dir,"Data/",sep = "")
outdir = paste(dir,"Results/INDEL/",sep = "")

#3 Functions
## Function 1: Get info from the metadatafile ##
# Use: columnname has to be between ""!
get_sampleinfo = function(list,columnname) {
  sample.info <- as.character()
  for (i in list) {
    i <- tail(unlist(strsplit(i,"/")), n=1)
    sample <- as.character(metadata[as.numeric(grep(i,metadata$indel_vcf_file)),columnname])
    sample.info <- c(sample.info, sample)
  }
  unlist(sample.info)
  return(sample.info)
}

#4 Install and load human reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#biocLite(ref_genome)
library(ref_genome, character.only = T)

#5 Metadata
metadata <- read.delim(paste(indir,"metadata.txt",sep=""), sep = "\t",header=TRUE)




# ---- GET VCFS ----

#1 ALC liver file locations + sample names
alc_vcf_files <- list.files(paste(indir,"INDEL/ALC/",sep=""),full.names = T)
alc_sample_names <-get_sampleinfo(alc_vcf_files,"sample_id")

#2  Healthy liver file locations + sample names
healthyliver_vcf_files <- list.files(paste(indir,"INDEL/Healthy_Liver/",sep=""), full.names = T)
healthyliver_sample_names <- get_sampleinfo(healthyliver_vcf_files,"sample_id")

#3 Combined 
vcf_files <- c(healthyliver_vcf_files,alc_vcf_files)
vcf_sample_names <- get_sampleinfo(vcf_files,"sample_id")

#4 Read vcfs
vcfs <- read_vcfs_as_granges(vcf_files,vcf_sample_names,genome = ref_genome)

#5 only select autosomal chromosomes
auto <- extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="auto")
vcfs <- lapply(vcfs, function(x) keepSeqlevels(x, auto))




# ---- GET & ORDER SAMPLE INFO ----

#1 Sample type
# Healthy, ALC, trunk HCC
healthyliver_type <- get_sampleinfo(healthyliver_vcf_files,"tissuetype")
alc_type <- get_sampleinfo(alc_vcf_files,"tissuetype")
type <- factor(c(healthyliver_type,alc_type),levels = c("Liver","Alc"))

#2 Age
# Healthy, ALC, trunk HCC
healthyliver_age <- get_sampleinfo(healthyliver_vcf_files,"age")
alc_age <- get_sampleinfo(alc_vcf_files,"age")
age <- c(healthyliver_age,alc_age)

#3 Prepare metadata for ordering: in same order as data
temp <- metadata
metadata <- data.frame()
for (i in vcf_files) {
  i <- tail(unlist(strsplit(i,"/")), n=1)
  metadata <- rbind(metadata,temp[as.numeric(grep(i,temp$indel_vcf_file)),])
}
remove(temp,i)

#4 Order the samples according to tissue type & age for later on
metadata$tissuenumber <- metadata$tissuetype
metadata$tissuenumber <- gsub("Alc",2,metadata$tissuenumber)
metadata$tissuenumber <- gsub("Liver",1,metadata$tissuenumber)
order_typeagesample <- order(metadata$tissuenumber, metadata$age, metadata$sample_id)
ordered_metadata <- metadata[order_typeagesample,]

#5 Order the sample info retrieved in earlier steps
# Sample type
ordered_type <- as.vector(ordered_metadata$tissuetype)
# Age
ordered_age <- ordered_metadata$age
# Sample names
ordered_sample_names <- as.vector(ordered_metadata$sample_id)




# ------ MUTATION NUMBER & RATE ------

#1 Load surveyed file
snvcorrected <- read.delim(paste(indir,"/callable/autosomal_callable_livers.txt",sep=""),header=T, sep="\t")
snvcorrected$Surveyed.percentage <- as.numeric(gsub("%", "", unlist(snvcorrected$Surveyed.percentage)))
snvcorrected <- snvcorrected[-24,]

#2 Add additional info
snvcorrected$sample_type <- type
snvcorrected$donor <- NA
for (i in 1:nrow(snvcorrected)) {
  snvcorrected[i,]$donor <- as.vector(metadata[which(metadata$sample_id == as.character(snvcorrected[i,]$Sample)),]$donor_id)
}
remove(i)
snvcorrected$age <- NA
for (i in 1:nrow(snvcorrected)) {
  snvcorrected[i,]$age <- as.numeric(metadata[which(metadata$sample_id == as.character(snvcorrected[i,]$Sample)),]$age)
}
remove(i)

#3 Extract SNV numbers from vcf files
snvcorrected$uncorr.indel <- NA
for (i in vcf_files) {
  df <- read.table(i)
  i <- tail(unlist(strsplit(i, '/', fixed = T)),n=1)
  sample <- as.character(metadata[as.numeric(grep(i,metadata$indel_vcf_file)),]$sample_id)
  snv_number <- nrow(df)
  
  #3B Add to the surveyed file
  snvcorrected[which(snvcorrected$Sample == sample),]$uncorr.indel <- snv_number
  
  remove(snv_number,sample,df)
}
remove(i)

#4 Correct mutationnumber for surveyed area (extrapolate INDELs to autosomal genome)
snvcorrected$corr.indel <- (snvcorrected$uncorr.indel/snvcorrected$Surveyed.percentage)*100

#5 Calculate number of snvs/year
snvcorrected$corr.indel.peryear <- snvcorrected$corr.indel/snvcorrected$age

#6 Save table
#write.table(snvcorrected,file = paste(outdir,"indelnumber.txt",sep=""),sep ="\t",col.names = T,row.names = F)




# ---- STATISTICAL ANALYSIS: INDEL rate ----

#1 linear mixed model
lme_liver = lme(corr.indel ~  age, random = ~ - 1 + age | donor, data=snvcorrected, subset=sample_type=="Liver")
lme_alc = lme(corr.indel ~  age, random = ~ - 1 + age | donor, data=snvcorrected, subset=sample_type=="Alc")

#2 P-value
# Healthy liver
summary(lme_liver)$tTable["age","p-value"]
# Alcoholic liver
summary(lme_alc)$tTable["age","p-value"] 

#3 95% CI
#3A Calculate
age_confint_liver = intervals(lme_liver)$fixed["age",]
age_confint_alc = intervals(lme_alc)$fixed["age",]
#3B put in df
age_confint = as.data.frame(rbind(age_confint_liver,age_confint_alc))
age_confint$sample_type = factor(c("Liver","Alc"), levels=c("Liver","Alc"))

#4 df with linear fits
# create data.frame with linear fits of fixed effect
newdat1 = expand.grid(sample_type="Liver", age = c(min(subset(snvcorrected, sample_type=="Liver")$age), max(subset(snvcorrected, sample_type=="Liver")$age)))
newdat1$fit = predict(lme_liver, level=0, newdata=newdat1)
newdat2  = expand.grid(sample_type="Alc", age = c(min(subset(snvcorrected, sample_type=="Alc")$age), max(subset(snvcorrected, sample_type=="Alc")$age)))
newdat2$fit = predict(lme_alc, level=0, newdata=newdat2)




# ---- PLOT: INDEL rate ----

#1 INDEL number vs age
indel.plot <- ggplot(snvcorrected, aes(x=age, y=corr.indel, color=sample_type)) +
  geom_line(data=newdat1, aes(y=fit, x=age), size=1.5) +
  geom_line(data=newdat2, aes(y=fit, x=age), size=1.5) +
  geom_point(size=2,color = c(rep("#00BFC4",15),rep("#F8766D",8))) +
  geom_point(shape=1, size=2, colour="black") +
  expand_limits(x=c(0,100)) +
  expand_limits(y=c(0,150)) +
  ylab("Indels per autosomal genome") +
  theme_bw() +
  xlab("Age (years)")
#ggsave(paste(outdir,"mutation_rate.pdf",sep=""), plot =indel.plot, width = 6, height = 4)

#2 Slope
slope_estimate_plot = ggplot(age_confint, aes(x=sample_type, y=est., fill=sample_type)) +
  geom_bar(colour="black", stat="identity") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.15) +
  ylab("Indels per autosomal genome per year") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(breaks=NULL)
#ggsave(paste(outdir,"mutation_rate_slopeestimate.pdf",sep=""), plot =slope_estimate_plot, width = 6, height = 4)