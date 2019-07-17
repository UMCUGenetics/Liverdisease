# @Date: 05 april 2018
# @Author: Myrthe Jager
# @Modified: 
# @Description: Genome-wide mutation pattern analysis healthy, alcoholic predisposed and alcoholic HCC
# Abbreviations: ALC = alcoholic




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
outdir = paste(dir,"Results/SNV/",sep = "")
beddir = ""

#3 Functions
## Function 1: Get info from the metadatafile ##
# Use: columnname has to be between ""!
get_sampleinfo = function(list,columnname) {
  sample.info <- as.character()
  for (i in list) {
    i <- tail(unlist(strsplit(i,"/")), n=1)
    sample <- as.character(metadata[as.numeric(grep(i,metadata$vcf_file)),columnname])
    sample.info <- c(sample.info, sample)
  }
  unlist(sample.info)
  return(sample.info)
}

## Function 2: PER TYPE: STATISTICS ##
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Function 3: Rename CHROM to UCSC style ##
rename_chrom = function(granges, style = "UCSC")
{
  # rename mitochondrial DNA manually
  seqlevels(granges)[seqlevels(granges)=="chrMT"] = "chrM"
  
  # get chromosome style
  chrom_style = mapSeqlevels(seqlevels(granges), style)
  
  # removing NA cases
  chrom_style = chrom_style[complete.cases(chrom_style)] 
  
  # rename chromosome names (seqlevels)
  res = renameSeqlevels(granges, chrom_style)
  
  return(res)
}

## Function 4: Convert bed to Granges ##
bed_to_granges = function(bed_files, region_names) 
{ 
  if (length(bed_files) != length(region_names)) 
    stop("Provide the same number of names as bed files") 
  
  granges_list = list() 
  for(i in 1:length(bed_files)) 
  { 
    bed_file = bed_files[i] 
    bed = read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
    chr = paste("chr", bed[,1], sep="") 
    
    # Convert BED (0-based) start postion to Granges (1-based) 
    start = bed[,2] + 1 
    
    # In BED end position is excluded, in Granges end position is 
    # included -> +1 -1 -> no conversion needed 
    end = bed[,3] 
    
    new_bed = GRanges(chr, IRanges(start,end)) 
    new_bed = list(new_bed) 
    names(new_bed) = region_names[i] 
    granges_list = c(granges_list, new_bed) 
  } 
  
  return(granges_list) 
}

## Function 5: Calculate intermutation distance ##
mut_dist = function(vcf){
  # mutation characteristics
  type = loc = dist = chrom = previous = c()
  
  # for each chromosome
  for(i in 1:length(chromosomes)){
    chr_subset = vcf[seqnames(vcf) == chromosomes[i]]
    n = length(chr_subset)
    if(n<=1){next}
    type = c(type, mut_type(chr_subset)[-1])
    loc = c(loc, (start(chr_subset))[-1])# + chr_cum[i])[-1])
    dist = c(dist, diff(start(chr_subset)))
    chrom = c(chrom, rep(chromosomes[i],n-1))
  }
  
  data = data.frame(type = type,
                    location = loc,
                    distance = dist,
                    chromosome = chrom)
  return(data)
}


#4 Install and load human reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#biocLite(ref_genome)
library(ref_genome, character.only = T)

#5 Metadata
metadata <- read.delim(paste(indir,"metadata.txt",sep=""), sep = "\t",header=TRUE)




# ---- GET VCFS ----

#1 ALC liver file locations + sample names
alc_vcf_files <- list.files(paste(indir,"SNV/ALC/IAP/",sep=""),full.names = T)
alc_sample_names <-get_sampleinfo(alc_vcf_files,"sample_id")

#2  Healthy liver file locations + sample names
healthyliver_vcf_files <- list.files(paste(indir,"SNV/Healthy_Liver/",sep=""), full.names = T)
healthyliver_sample_names <- get_sampleinfo(healthyliver_vcf_files,"sample_id")

#3  Liver cancer
# Clonal file locations + sample names
hcc_clonal_vcf_file <- list.files(paste(indir,"SNV/HCC/clonal/",sep=""), full.names = T)
hcc_clonal_sample_names <- get_sampleinfo(hcc_clonal_vcf_file,"sample_id")

#4 Combined
vcf_files <- c(healthyliver_vcf_files,alc_vcf_files,hcc_clonal_vcf_file )
vcf_sample_names <- get_sampleinfo(vcf_files,"sample_id")

#5 Read vcfs
vcfs <- read_vcfs_as_granges(vcf_files,vcf_sample_names,genome = ref_genome)

#6 only select autosomal chromosomes
auto <- extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="auto")
vcfs <- lapply(vcfs, function(x) keepSeqlevels(x, auto))




# ---- GET & ORDER SAMPLE INFO ----

#1 Sample type
# Healthy, ALC, trunk HCC
healthyliver_type <- get_sampleinfo(healthyliver_vcf_files,"tissuetype")
alc_type <- get_sampleinfo(alc_vcf_files,"tissuetype")
hcc_clonal_type <- get_sampleinfo(hcc_clonal_vcf_file,"tissuetype")
type <- factor(c(healthyliver_type,alc_type,hcc_clonal_type),levels = c("Liver","Alc","HCC"))

#2 Age
# Healthy, ALC, trunk HCC
healthyliver_age <- get_sampleinfo(healthyliver_vcf_files,"age")
alc_age <- get_sampleinfo(alc_vcf_files,"age")
hcc_clonal_age <- get_sampleinfo(hcc_clonal_vcf_file,"age")
age <- c(healthyliver_age,alc_age,hcc_clonal_age)

#3 Prepare metadata for ordering: in same order as data
temp <- metadata
metadata <- data.frame()
for (i in vcf_files) {
  i <- tail(unlist(strsplit(i,"/")), n=1)
  metadata <- rbind(metadata,temp[as.numeric(grep(i,temp$vcf_file)),])
}
remove(temp,i)

#4 Order the samples according to tissue type & age for later on
metadata$tissuenumber <- metadata$tissuetype
metadata$tissuenumber <- gsub("Alc",2,metadata$tissuenumber)
metadata$tissuenumber <- gsub("Liver",1,metadata$tissuenumber)
metadata$tissuenumber <- gsub("HCC",3,metadata$tissuenumber)
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
snvcorrected[24,]$Surveyed.percentage <- 100

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
snvcorrected$uncorr.snv <- NA
for (i in vcf_files) {
  df <- read.table(i)
  sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[1]
  snv_number <- nrow(df)
  sampleid <- as.character(metadata[which(metadata$sample_name == sample),]$sample_id)
  
  #3B Add to the surveyed file
  snvcorrected[which(snvcorrected$Sample == sampleid),]$uncorr.snv <- snv_number
  
  remove(sampleid,snv_number,sample,df)
}
remove(i)

#4 Add dinucleotides to the table
#4A get cumulative sum of chromosome lengths
chr_length = seqlengths(vcfs[[1]])
chromosomes = names(chr_length)
chr_cum = c(0, cumsum(as.numeric(chr_length)))
#4B Calculate intermutation distance
mut_dist_vcfs = lapply(vcfs, function(x) mut_dist(x))
#4C Generate empty dfs
double.df <- data.frame(sample = NA, number = NA)
double.df.type <- data.frame(type=NA,location=NA,sample=NA,chromosome=NA)
#4D Extract dinucleotide counts and add to dfs
for(i in 1:length(mut_dist_vcfs)) {
  df.temp <- as.data.frame(mut_dist_vcfs[i])
  sample <- names(mut_dist_vcfs[i])
  df.temp.double <-df.temp[which(df.temp[,3] ==1),]
  number <- nrow(df.temp.double)
  
  df.temp.double[,3] <- rep(sample,nrow(df.temp.double))
  colnames(df.temp.double) <- c("type","location","sample","chromosome")
  
  double.df.type <- rbind(double.df.type,df.temp.double)
  double.df <- rbind(double.df,c(sample,number))
  
  remove(df.temp,sample,number,df.temp.double)
}
remove(i)
double.df<-double.df[-1,]
double.df.type<-double.df.type[-1,]
#4E Add to df
snvcorrected$uncorr.double =NA
for(i in 1:nrow(snvcorrected)) {
  snvcorrected[i,]$uncorr.double <- as.numeric(double.df[which(double.df$sample == snvcorrected[i,]$Sample),]$number)
}
remove(i)

#5 Correct mutationnumber for surveyed area (extrapolate SNVs to autosomal genome)
snvcorrected$corr.snv <- (snvcorrected$uncorr.snv/snvcorrected$Surveyed.percentage)*100
snvcorrected$corr.single <- ((snvcorrected$uncorr.snv-(snvcorrected$uncorr.double)*2)/snvcorrected$Surveyed.percentage)*100
snvcorrected$corr.double <- (snvcorrected$uncorr.double/snvcorrected$Surveyed.percentage)*100

#6 Calculate number of snvs/year
snvcorrected$corr.snv.peryear <- snvcorrected$corr.snv/snvcorrected$age

#7 Add info for single
snvcorrected$uncorr.single <- snvcorrected$uncorr.snv-(snvcorrected$uncorr.double)*2

#8 Save table
snvcorrected[24,]$Surveyed <- NA
snvcorrected[24,]$Surveyed.percentage <- NA
#write.table(snvcorrected,file = paste(outdir,"snvnumber.txt",sep=""),sep ="\t",col.names = T,row.names = F)

#9 Calculate total number of identified snvs
# ALC
sum(snvcorrected[which(snvcorrected$sample_type == "Alc"),]$uncorr.snv)
# Healthy liver
sum(snvcorrected[which(snvcorrected$sample_type == "Liver"),]$uncorr.snv)
# New healthy liver samples + ALC
sum(snvcorrected[which(snvcorrected$sample_type == "Alc"),]$uncorr.snv) + snvcorrected[which(snvcorrected$donor == 20),]$uncorr.snv + sum(snvcorrected[which(snvcorrected$donor == 21),]$uncorr.snv) + 
  snvcorrected[which(snvcorrected$donor == 22),]$uncorr.snv + snvcorrected[which(snvcorrected$donor == 23),]$uncorr.snv
# All healthy liver samples + ALC
sum(snvcorrected[which(snvcorrected$sample_type != "HCC"),]$uncorr.snv)




# ---- HCC: MUTATION NUMBER (PER BIOPSY) ----

hcc_perclone_vcf_file <- list.files(paste(indir,"SNV/HCC/per_clone/",sep=""), full.names = T)
hcc.df <- data.frame()
for(i in hcc_perclone_vcf_file) {
  df <- read.table(i)
  sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[2]
  snv_number <- nrow(df)
  new_df <- data.frame(sample = sample,
                       mutations = snv_number)
  hcc.df <- rbind(hcc.df,new_df)
  remove(df,sample,snv_number,new_df)
}
remove(i)
#write.table(hcc.df,file = paste(outdir,"HCC/number_perclone.txt",sep=""),sep ="\t",col.names = T,row.names = F)




# ---- STATISTICAL ANALYSIS: SNV rate ----

#1 linear mixed model
lme_liver = lme(corr.snv ~  age, random = ~ - 1 + age | donor, data=snvcorrected, subset=sample_type=="Liver")
lme_alc = lme(corr.snv ~  age, random = ~ - 1 + age | donor, data=snvcorrected, subset=sample_type=="Alc")

#2 P-value
# Healthy liver
summary(lme_liver)$tTable["age","p-value"]
# Alcoholic liver
summary(lme_alc)$tTable["age","p-value"] # --> no significant linear accumulation with age

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

#5 Calculate whether Alcoholic is within 95% CI of healthy
# Generate df with all info
table_binom <- data.frame(sample = snvcorrected[which(snvcorrected$sample_type == "Alc"),"Sample"],
                          age = snvcorrected[which(snvcorrected$sample_type == "Alc"),"age"],
                          prob = summary(lme_liver)$tTable["age","Value"],
                          observed = snvcorrected[which(snvcorrected$sample_type == "Alc"),"corr.snv"],
                          expected = summary(lme_liver)$tTable["age","Value"] * snvcorrected[which(snvcorrected$sample_type == "Alc"),"age"],
                          lower = age_confint_liver[1]* snvcorrected[which(snvcorrected$sample_type == "Alc"),"age"],
                          upper = age_confint_liver[3]* snvcorrected[which(snvcorrected$sample_type == "Alc"),"age"]) 
# Calculate whether higher than low CI and lower than high CI
table_binom$higherthanlower <- table_binom$observed-table_binom$lower
table_binom$lowerthanupper <- table_binom$upper-table_binom$observed
# Conclusion
table_binom$within.ci <- rep("NO",nrow(table_binom))
table_binom[which(table_binom$higherthanlower > 0 & table_binom$lowerthanupper > 0),]$within.ci <- rep("YES",nrow(table_binom[which(table_binom$higherthanlower > 0 & table_binom$lowerthanupper > 0),]))
# All values are within 95% CI of healthy!



# ---- PLOT: SNV rate ----

snv.plot <- ggplot(snvcorrected[which(snvcorrected$sample_type != "HCC"),], aes(x=age, y=corr.snv, color=sample_type)) +
  geom_line(data=newdat1, aes(y=fit, x=age), size=1.5) +
  geom_point(size=2,color = c(rep("#F8766D",15),rep("#00BFC4",8))) +
  geom_point(shape=1, size=2, colour="black") +
  expand_limits(x=c(0,100)) +
  expand_limits(y=c(0,3500)) +
  ylab("No. point mutations per autosomal genome") +
  theme_bw() +
  xlab("Age (years)")
#ggsave(paste(outdir,"mutation_rate.pdf",sep=""), plot =snv.plot, width = 6, height = 4)

snv.plot.hcc <- ggplot(snvcorrected, aes(x=age, y=corr.snv, color=sample_type)) +
  geom_line(data=newdat1, aes(y=fit, x=age), size=1.5) +
  geom_point(size=2,color = c(rep("#F8766D",15),rep("#00BFC4",8),"#C77CFF")) +
  geom_point(shape=1, size=2, colour="black") +
  expand_limits(x=c(0,100)) +
  expand_limits(y=c(0,6500)) +
  ylab("No. point mutations per autosomal genome") +
  theme_bw() +
  xlab("Age (years)")
#ggsave(paste(outdir,"mutation_rate_withhcc.pdf",sep=""), plot =snv.plot.hcc, width = 6, height = 4)

slope_estimate_plot = ggplot(age_confint, aes(x=sample_type, y=est., fill=sample_type)) +
  geom_bar(colour="black", stat="identity") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.15) +
  ylab("No. point mutations per autosomal genome") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(breaks=NULL)
#ggsave(paste(outdir,"mutation_rate_slopeestimate.pdf",sep=""), plot =slope_estimate_plot, width = 6, height = 4)




# ---- Dinucleotide rate: STATISTICS AND PLOT ----

#1 LME
lme_liver_double = lme(corr.double ~  age, random = ~ - 1 + age | donor, data=snvcorrected, subset=sample_type=="Liver")
lme_alc_double = lme(corr.double ~  age, random = ~ - 1 + age | donor, data=snvcorrected, subset=sample_type=="Alc")

#2 P-value
# Healthy liver
summary(lme_liver_double)$tTable["age","p-value"]
# Alc
summary(lme_alc_double)$tTable["age","p-value"]  # --> no significant linear accumulation with age

#3 95% CI
#3A Calculate
age_confint_liver_double = intervals(lme_liver_double)$fixed["age",]
age_confint_alc_double = intervals(lme_alc_double)$fixed["age",]
#3B put in df
age_confint_double = as.data.frame(rbind(age_confint_liver_double,age_confint_alc_double))
age_confint_double$sample_type = factor(c("Liver","Alc"), levels=c("Liver","Alc"))

#4 df with linear fits
# create data.frame with linear fits of fixed effect
newdat1_double = expand.grid(sample_type="Liver", age = c(min(subset(snvcorrected, sample_type=="Liver")$age), max(subset(snvcorrected, sample_type=="Liver")$age)))
newdat1_double$fit = predict(lme_liver_double, level=0, newdata=newdat1_double)

#5 Calculate whether Alcoholic is within 95% CI of healthy
# Generate df with all info
table_binom_double <- data.frame(sample = snvcorrected[which(snvcorrected$sample_type == "Alc"),"Sample"],
                          age = snvcorrected[which(snvcorrected$sample_type == "Alc"),"age"],
                          prob = summary(lme_liver_double)$tTable["age","Value"],
                          observed = snvcorrected[which(snvcorrected$sample_type == "Alc"),"corr.double"],
                          expected = summary(lme_liver_double)$tTable["age","Value"] * snvcorrected[which(snvcorrected$sample_type == "Alc"),"age"],
                          lower = age_confint_liver_double[1] * snvcorrected[which(snvcorrected$sample_type == "Alc"),"age"],
                          upper = age_confint_liver_double[3] * snvcorrected[which(snvcorrected$sample_type == "Alc"),"age"]) 
# Calculate whether higher than low CI and lower than high CI
table_binom_double$higherthanlower <- table_binom_double$observed-table_binom_double$lower
table_binom_double$lowerthanupper <- table_binom_double$upper-table_binom_double$observed
# Conclusion
table_binom_double$within.ci <- rep("NO",nrow(table_binom_double))
table_binom_double[which(table_binom_double$higherthanlower > 0 & table_binom_double$lowerthanupper > 0),]$within.ci <- rep("YES",nrow(table_binom_double[which(table_binom_double$higherthanlower > 0 & table_binom_double$lowerthanupper > 0),]))
# One value is lower than 95% CI; rest = within 95% CI

#6 plot
double.plot <- ggplot(snvcorrected[which(snvcorrected$sample_type != "HCC"),], aes(x=age, y=corr.double, color=sample_type)) +
  geom_line(data=newdat1_double, aes(y=fit, x=age), size=1.5) +
  geom_point(size=2,color = c(rep("#F8766D",15),rep("#00BFC4",8))) +
  geom_point(shape=1, size=2, colour="black") +
  expand_limits(x=c(0,100)) +
  expand_limits(y=c(0,30)) +
  ylab("No. tandem point mutations per autosomal genome") +
  theme_bw() +
  xlab("Age (years)")
#ggsave(paste(outdir,"mutation_rate_tandem.pdf",sep=""), plot =double.plot, width = 6, height = 4)




# ---- MUTATION SPECTRA & PROFILES ----

#1 Get 6 types
type_occurrences <- mut_type_occurrences(vcfs,ref_genome)

#2 Plot relative spectrum per 6 mutation types 
#2A Plot 
spectrum6typeplotwithCT = plot_spectrum(type_occurrences[-24,], by = type[-24], legend = T, CT = T)
spectrum6typeplotwithCT.hcc = plot_spectrum(type_occurrences[24,], by = type[24], legend = T, CT = T)
#2B Save plot
#ggsave(paste(outdir,"Spectrum_6type.pdf",sep=""), spectrum6typeplotwithCT, width = 8, height = 4)
#ggsave(paste(outdir,"Spectrum_6type_HCC.pdf",sep=""), spectrum6typeplotwithCT.hcc, width = 4, height = 4)

#3 Get 96 mutation types
mut_matrix_notordered = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
                  
#4 Order mutation matrix (96 types)
mut_matrix <- mut_matrix_notordered[,order_typeagesample]

#5 Collapse 96 mutation types matrix
collapsed_mutmatrix <- data.frame(rowSums(mut_matrix[,c(1:15)], na.rm = FALSE, dims = 1),
                                      rowSums(mut_matrix[,c(16:23)], na.rm = FALSE, dims = 1),
                                      mut_matrix[,24]
)
colnames(collapsed_mutmatrix) <- c("Healthy_liver","Alcoholic_liver","HCC")

#6 Adjust collapsed for sample size
collapsed_combined_mutmatrix <- data.frame(rowSums(mut_matrix[,c(1:15)], na.rm = FALSE, dims = 1)/15,
                                           rowSums(mut_matrix[,c(16:23)], na.rm = FALSE, dims = 1)/8,
                                           mut_matrix[,24]
)
colnames(collapsed_combined_mutmatrix) <- c("Healthy_liver","Alcoholic_liver","HCC")

#7 PLot 96 mutation types
mutationprofileplot <- plot_96_profile(collapsed_mutmatrix, ymax = 0.05, condensed = T)
mutationprofileplot.separate <- plot_96_profile(mut_matrix, ymax = 0.1, condensed = T)
#ggsave(paste(outdir,"Spectrum_96type_ymax0dot05.pdf",sep=""), plot = mutationprofileplot, height = 6, width = 8)
#ggsave(paste(outdir,"Spectrum_96type_perASC.pdf",sep=""), plot = mutationprofileplot.separate, height = 25, width = 8)
#write.table(mut_matrix,file =paste(outdir,"mutational_profile.txt",sep=""),sep ="\t", row.names = F, col.names = T)




# ---- COSINE SIMILARITY OF Mutational spectra ----

#1 Calculate cosine similarity
# Collapsed healthy, alc and HCC trunk
cosine.collapsed <- cos_sim_matrix(collapsed_mutmatrix,collapsed_mutmatrix)
# All ASCs healthy, alc + HCC trunk
cosine.all <- cos_sim_matrix(mut_matrix,mut_matrix)

#2 Plot
cosine.plot.collapsed <- plot_cosine_heatmap(cosine.collapsed,plot_values = T)
cosine.plot.all <- plot_cosine_heatmap(cosine.all,plot_values = T)

#3 Save plots
#ggsave(paste(outdir,"Cosinesimilarity_combinedprofiles.pdf",sep=""),plot = cosine.plot.collapsed, height = 4, width = 8)
#ggsave(paste(outdir,"Cosinesimilarity_perstemcell.pdf",sep=""),plot = cosine.plot.all, height = 8, width = 12)

#4 Compare profiles
plot_compare_profiles(collapsed_combined_mutmatrix[,1],collapsed_combined_mutmatrix[,2], condensed = TRUE,profile_ymax = 0.05, profile_names = colnames(collapsed_combined_mutmatrix[1:2]))
plot_compare_profiles(collapsed_combined_mutmatrix[,1],collapsed_combined_mutmatrix[,3], condensed = TRUE,profile_ymax = 0.05, profile_names = colnames(collapsed_combined_mutmatrix[c(1,3)]))
plot_compare_profiles(collapsed_combined_mutmatrix[,2],collapsed_combined_mutmatrix[,3], condensed = TRUE,profile_ymax = 0.05, profile_names = colnames(collapsed_combined_mutmatrix[2:3]))




# ---- REFIT SPECTRA ----

#1 Load COSMIC mutational signatures
#1A v2 (downloaded march 13th, 2018)
cancer_signatures_v2 = read.table(paste(indir,"SNV/cancersignaturesv2.txt",sep=""), sep="\t")
#1B v3 (downloaded may 22nd, 2019)
cancer_signatures_v3 = read.csv(paste(indir,"SNV/cancersignaturesv3.csv",sep=""))

#2 Reorder (to make the order of the trinucleotide changes the same)
#v2
cancer_signatures_v2 <- cancer_signatures_v2[order(cancer_signatures_v2[,1]),]
#check: rownames(mut_matrix) == cancer_signatures_v2$Somatic.Mutation.Type
#v3
#check: rownames(mut_matrix) == paste(paste(substring(cancer_signatures_v3$SubType,1,1),cancer_signatures_v3$Type,sep="["),substring(cancer_signatures_v3$SubType,3,3),sep="]")
#not necessary for the new ones

#3 Only signatures in matrix
cancer_signatures_v2 <- as.matrix(cancer_signatures_v2[,4:33])
cancer_signatures_v3 <- as.matrix(cancer_signatures_v3[,3:ncol(cancer_signatures_v3)])

#4 Refit
#4A v2
fit_res_v2 = fit_to_signatures(mut_matrix[,-24], cancer_signatures_v2)
fit_res.withhcc_v2 = fit_to_signatures(mut_matrix, cancer_signatures_v2)
fit_res_collapsed_v2 <- fit_to_signatures(collapsed_combined_mutmatrix, cancer_signatures_v2)
#4B v3
fit_res_v3 = fit_to_signatures(mut_matrix[,-24], cancer_signatures_v3)
fit_res.withhcc_v3 = fit_to_signatures(mut_matrix, cancer_signatures_v3)
fit_res_collapsed_v3 <- fit_to_signatures(collapsed_combined_mutmatrix, cancer_signatures_v3)

#5 Select signatures with some contribution
# v2
select_v2 = which(rowSums(fit_res_v2$contribution) > 100)
select.withhcc_v2 = which(rowSums(fit_res.withhcc_v2$contribution) > 100)
select_collapsed_v2 = which(rowSums(fit_res_collapsed_v2$contribution) > 20)
# v3
select_v3 = which(rowSums(fit_res_v3$contribution) > 100)
select.withhcc_v3 = which(rowSums(fit_res.withhcc_v3$contribution) > 100)
select_collapsed_v3 = which(rowSums(fit_res_collapsed_v3$contribution) > 100)

#6 Plot contribution
# v2
rel.contribution.cosmic.plot_v2 <- plot_contribution(fit_res_v2$contribution[select_v2,], coord_flip = T, signatures = cancer_signatures_v2, mode = "relative")
abs.contribution.cosmic.plot_v2 <- plot_contribution(fit_res_v2$contribution[select_v2,], coord_flip = T, signatures = cancer_signatures_v2[,select_v2], mode = "absolute")
rel.contribution.cosmic.collapsed.plot_v2 <- plot_contribution(fit_res_collapsed_v2$contribution[select_collapsed_v2,], coord_flip = T, signatures = cancer_signatures_v2, mode = "relative")
abs.contribution.cosmic.collapsed.plot_v2 <- plot_contribution(fit_res_collapsed_v2$contribution[select_collapsed_v2,], coord_flip = T, signatures = cancer_signatures_v2[,select_collapsed_v2], mode = "absolute")
# v3
rel.contribution.cosmic.plot_v3 <- plot_contribution(fit_res_v3$contribution[select_v3,], coord_flip = T, signatures = cancer_signatures_v3, mode = "relative")
abs.contribution.cosmic.plot_v3 <- plot_contribution(fit_res_v3$contribution[select_v3,], coord_flip = T, signatures = cancer_signatures_v3[,select_v3], mode = "absolute")
rel.contribution.cosmic.collapsed.plot_v3 <- plot_contribution(fit_res_collapsed_v3$contribution[select_collapsed_v3,], coord_flip = T, signatures = cancer_signatures_v3, mode = "relative")
abs.contribution.cosmic.collapsed.plot_v3 <- plot_contribution(fit_res_collapsed_v3$contribution[select_collapsed_v3,], coord_flip = T, signatures = cancer_signatures_v3[,select_collapsed_v3], mode = "absolute")

#7 Save plots
#ggsave(paste(outdir,"COSMIC_refit-absolute_v2.pdf",sep=""), plot = abs.contribution.cosmic.plot_v2, width = 10, height = 6)
#ggsave(paste(outdir,"COSMIC_refit-relative_v2.pdf",sep=""), plot = rel.contribution.cosmic.plot_v2, width = 10, height = 6)
#ggsave(paste(outdir,"COSMIC_refit-absolute_collapsed_v2.pdf",sep=""), plot = abs.contribution.cosmic.collapsed.plot_v2, width = 10, height = 6)
#ggsave(paste(outdir,"COSMIC_refit-relative_collapsed_v2.pdf",sep=""), plot = rel.contribution.cosmic.collapsed.plot_v2, width = 10, height = 6)
#ggsave(paste(outdir,"COSMIC_refit-absolute_v3.pdf",sep=""), plot = abs.contribution.cosmic.plot_v3, width = 10, height = 6)
#ggsave(paste(outdir,"COSMIC_refit-relative_v3.pdf",sep=""), plot = rel.contribution.cosmic.plot_v3, width = 10, height = 6)
#ggsave(paste(outdir,"COSMIC_refit-absolute_collapsed_v3.pdf",sep=""), plot = abs.contribution.cosmic.collapsed.plot_v3, width = 10, height = 6)
#ggsave(paste(outdir,"COSMIC_refit-relative_collapsed_v3.pdf",sep=""), plot = rel.contribution.cosmic.collapsed.plot_v3, width = 10, height = 6)

#8 Cosine similarity reconstructed to original
cos_sim_matrix(fit_res_collapsed_v2$reconstructed,collapsed_combined_mutmatrix)
cos_sim_matrix(fit_res_collapsed_v3$reconstructed,collapsed_combined_mutmatrix)




# ---- RAINFALL PLOT ----

#1 define chromosomes of interest
chromosomes = seqnames(get(ref_genome))[1:22]

#2 Plot rainfall
#Healthy
rainh1 <- plot_rainfall(vcfs[[1]],title = names(vcfs[1]), chromosomes = chromosomes, cex = 0.5)
rainh2 <- plot_rainfall(vcfs[[2]],title = names(vcfs[2]), chromosomes = chromosomes, cex = 0.5)
rainh3 <- plot_rainfall(vcfs[[3]],title = names(vcfs[3]), chromosomes = chromosomes, cex = 0.5)
rainh4 <- plot_rainfall(vcfs[[4]],title = names(vcfs[4]), chromosomes = chromosomes, cex = 0.5)
rainh5 <- plot_rainfall(vcfs[[5]],title = names(vcfs[5]), chromosomes = chromosomes, cex = 0.5)
rainh6 <- plot_rainfall(vcfs[[6]],title = names(vcfs[6]), chromosomes = chromosomes, cex = 0.5)
rainh7 <- plot_rainfall(vcfs[[7]],title = names(vcfs[7]), chromosomes = chromosomes, cex = 0.5)
rainh8 <- plot_rainfall(vcfs[[8]],title = names(vcfs[8]), chromosomes = chromosomes, cex = 0.5)
rainh9 <- plot_rainfall(vcfs[[9]],title = names(vcfs[9]), chromosomes = chromosomes, cex = 0.5)
rainh10 <- plot_rainfall(vcfs[[10]],title = names(vcfs[10]), chromosomes = chromosomes, cex = 0.5)
rainh11 <- plot_rainfall(vcfs[[11]],title = names(vcfs[11]), chromosomes = chromosomes, cex = 0.5)
rainh12 <- plot_rainfall(vcfs[[12]],title = names(vcfs[12]), chromosomes = chromosomes, cex = 0.5)
rainh13 <- plot_rainfall(vcfs[[13]],title = names(vcfs[13]), chromosomes = chromosomes, cex = 0.5)
rainh14 <- plot_rainfall(vcfs[[14]],title = names(vcfs[14]), chromosomes = chromosomes, cex = 0.5)
rainh15 <- plot_rainfall(vcfs[[15]],title = names(vcfs[15]), chromosomes = chromosomes, cex = 0.5)
rain.plot.healthy = grid.arrange(rainh1,rainh2,rainh3,rainh4,rainh5,rainh6, rainh7, rainh8,rainh9, rainh10, rainh11, rainh12, rainh13, rainh14, rainh15, nrow=3, ncol =5)
#Alc
raina1 <- plot_rainfall(vcfs[[16]],title = names(vcfs[16]), chromosomes = chromosomes, cex = 0.5)
raina2 <- plot_rainfall(vcfs[[17]],title = names(vcfs[17]), chromosomes = chromosomes, cex = 0.5)
raina3 <- plot_rainfall(vcfs[[18]],title = names(vcfs[18]), chromosomes = chromosomes, cex = 0.5)
raina4 <- plot_rainfall(vcfs[[19]],title = names(vcfs[19]), chromosomes = chromosomes, cex = 0.5)
raina5 <- plot_rainfall(vcfs[[20]],title = names(vcfs[20]), chromosomes = chromosomes, cex = 0.5)
raina6 <- plot_rainfall(vcfs[[21]],title = names(vcfs[21]), chromosomes = chromosomes, cex = 0.5)
raina7 <- plot_rainfall(vcfs[[22]],title = names(vcfs[22]), chromosomes = chromosomes, cex = 0.5)
raina8 <- plot_rainfall(vcfs[[23]],title = names(vcfs[23]), chromosomes = chromosomes, cex = 0.5)
rain.plot.alc = grid.arrange(raina1,raina2,raina3,raina4,raina5,raina6, raina7, raina8, nrow=2, ncol =4)
#HCC clonal
rain.plot.hccclonal <- plot_rainfall(vcfs[[24]],title = names(vcfs[24]), chromosomes = chromosomes, cex = 0.5)

#3 Save plots
#ggsave(paste(outdir,"Rainfall_liver.pdf",sep=""),plot = rain.plot.healthy, width = 20, height = 10)
#ggsave(paste(outdir,"Rainfall_alc.pdf",sep=""),plot = rain.plot.alc, width = 12, height = 6)
#ggsave(paste(outdir,"Rainfall_hcc-clonal.pdf",sep=""),plot = rain.plot.hccclonal, width = 6, height = 4)




# ---- TRANSCRIPTIONAL STRAND BIAS ----

#1 Get knowngenes table from UCSC for genome
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

#2 Make mutation count matrix with transcriptional information
mut_mat_s_notordered <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
mut_mat_s <- mut_mat_s_notordered[,order_typeagesample]

#3 Plot strand bias
strand_counts = strand_occurrences(mut_mat_s, by=type)
strand_plot = plot_strand(strand_counts, mode = "relative")

#4 Plot log2
strand_bias = strand_bias_test(strand_counts)
strand_bias_plot = plot_strand_bias(strand_bias)

#5 Save plot
plot.tcbias =grid.arrange(strand_plot,strand_bias_plot, nrow=2, ncol =1)
#ggsave(paste(outdir,"TCstrandbias.pdf",sep=""),plot =  plot.tcbias, width = 12, height = 8)

#6 Statistical testing
#6A Healthy vs alc (C>A (Transcribed(T)/total C>A;Untranscribed(U)/total C>A), C>G (T;U), C>T (T;U), T>A (T;U), T>C (T;U), T>G (T;U))
p1 <- poisson.test(x = c(strand_counts[1,]$no_mutations,strand_counts[13,]$no_mutations), T = c(sum(strand_counts[1:2,]$no_mutations),sum(strand_counts[13:14,]$no_mutations)))$p.value #NS
p2 <- poisson.test(x = c(strand_counts[2,]$no_mutations,strand_counts[14,]$no_mutations), T = c(sum(strand_counts[1:2,]$no_mutations),sum(strand_counts[13:14,]$no_mutations)))$p.value #NS
p3 <- poisson.test(x = c(strand_counts[3,]$no_mutations,strand_counts[15,]$no_mutations), T = c(sum(strand_counts[3:4,]$no_mutations),sum(strand_counts[15:16,]$no_mutations)))$p.value #NS
p4 <- poisson.test(x = c(strand_counts[4,]$no_mutations,strand_counts[16,]$no_mutations), T = c(sum(strand_counts[3:4,]$no_mutations),sum(strand_counts[15:16,]$no_mutations)))$p.value #NS
p5 <- poisson.test(x = c(strand_counts[5,]$no_mutations,strand_counts[17,]$no_mutations), T = c(sum(strand_counts[5:6,]$no_mutations),sum(strand_counts[17:18,]$no_mutations)))$p.value #NS
p6 <- poisson.test(x = c(strand_counts[6,]$no_mutations,strand_counts[18,]$no_mutations), T = c(sum(strand_counts[5:6,]$no_mutations),sum(strand_counts[17:18,]$no_mutations)))$p.value #NS
p7 <- poisson.test(x = c(strand_counts[7,]$no_mutations,strand_counts[19,]$no_mutations), T = c(sum(strand_counts[7:8,]$no_mutations),sum(strand_counts[19:20,]$no_mutations)))$p.value #NS
p8 <- poisson.test(x = c(strand_counts[8,]$no_mutations,strand_counts[20,]$no_mutations), T = c(sum(strand_counts[7:8,]$no_mutations),sum(strand_counts[19:20,]$no_mutations)))$p.value #NS
p9 <- poisson.test(x = c(strand_counts[9,]$no_mutations,strand_counts[21,]$no_mutations), T = c(sum(strand_counts[9:10,]$no_mutations),sum(strand_counts[21:22,]$no_mutations)))$p.value #NS
p10 <- poisson.test(x = c(strand_counts[10,]$no_mutations,strand_counts[22,]$no_mutations), T = c(sum(strand_counts[9:10,]$no_mutations),sum(strand_counts[21:22,]$no_mutations)))$p.value #NS
p11 <- poisson.test(x = c(strand_counts[11,]$no_mutations,strand_counts[23,]$no_mutations), T = c(sum(strand_counts[11:12,]$no_mutations),sum(strand_counts[23:24,]$no_mutations)))$p.value #NS
p12 <- poisson.test(x = c(strand_counts[12,]$no_mutations,strand_counts[24,]$no_mutations), T = c(sum(strand_counts[11:12,]$no_mutations),sum(strand_counts[23:24,]$no_mutations)))$p.value #NS
#6B get significant before multiple testing
c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12) < 0.05 # all n.s.
#6C Adjust for multiple testing
p.adjust(c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12), method="fdr")
#6D Remove all p-values
remove(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)




# ---- GENOMIC DISTRIBUTION -----

#1 Choose mart
#listEnsembl()
listMarts()
mart="ensembl"

#2 Choose dataset
# list datasets available from ensembl (for hg19 = GrCh37)
listDatasets(useEnsembl(biomart="regulation", GRCh = 37))
# Multicell regulatory features for hg19
regulation_regulatory = useEnsembl(biomart="regulation", dataset="hsapiens_regulatory_feature", GRCh = 37)
# list all possible filters
listFilters(regulation_regulatory)
# list all posible output attributes
listAttributes(regulation_regulatory)
# list all filter options for a specific attribute
filterOptions("regulatory_feature_type_name", regulation_regulatory)

#3 Get regions
#3A Promoters
promoter = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                 filters = "regulatory_feature_type_name", 
                 values = "Promoter", 
                 mart = regulation_regulatory)
promoter_g = reduce(GRanges(promoter$chromosome_name, IRanges(promoter$chromosome_start, promoter$chromosome_end)))

#3B Promoter flanking regions
#promoter_flanking = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
#                          filters = "regulatory_feature_type_name", 
#                          values = "Promoter Flanking Region", 
#                          mart = regulation_regulatory)
#promoter_flanking_g = reduce(GRanges(promoter_flanking$chromosome_name, IRanges(promoter_flanking$chromosome_start, promoter_flanking$chromosome_end))) 

#3C Open chromatin
open_chromatin = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                       filters = "regulatory_feature_type_name", 
                       values = "Open chromatin", 
                       mart = regulation_regulatory)
open_chromatin_g = reduce(GRanges(open_chromatin$chromosome_name, IRanges(open_chromatin$chromosome_start, open_chromatin$chromosome_end))) 

#3D Enhancer
enhancer = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                 filters = "regulatory_feature_type_name", 
                 values = "Enhancer", 
                 mart = regulation_regulatory)
enhancer_g = reduce(GRanges(enhancer$chromosome_name, IRanges(enhancer$chromosome_start, enhancer$chromosome_end)))

#3E Genes
regionsg = list(genes_hg19)
names(regionsg) = c("Genes")
reduced_regionsg <- GenomicRanges::reduce(regionsg$Genes, ignore.strand=TRUE)

#4 Combine regions
regions = list(reduced_regionsg,promoter_g, enhancer_g,open_chromatin_g)
names(regions) = c("Genes","Promoter", "Enhancer", "Open chromatin")
regions = lapply(regions, function(x) rename_chrom(x))

#5 Get bedfiles (callable regions per sample)
# Locations
healthyliver_bed_files <- list.files(paste(beddir,"healthy_liver/", sep=""), full.names = T)
alc_bed_files <- list.files(paste(beddir,"ALC/", sep=""), full.names = T)
hcc_clonal_bed_files <- list.files(paste(beddir,"HCC/clonal/", sep=""), full.names = T)
bed.files <- c(healthyliver_bed_files,alc_bed_files,hcc_clonal_bed_files)
bed.sample_names <- get_sampleinfo(vcf_files,"sample_id")
# Open files
surveyed_list = bed_to_granges(bed.files, bed.sample_names)

#6 Calculate Log2(observed/expected) for all genomic regions
distr = genomic_distribution(vcfs, surveyed_list, regions)
distr_test = enrichment_depletion_test(distr, by = type)

#7 Plot genomic distribution
genomic_dist_plot <- plot_enrichment_depletion(distr_test)
genomic_dist_plot.nohcc <- plot_enrichment_depletion(distr_test[which(distr_test$by != "HCC"),])
genomic_dist_plot_small <- plot_enrichment_depletion(distr_test[c(-1,-2,-3),])
#ggsave(paste(outdir,"Genomicdistribution.pdf",sep=""), plot = genomic_dist_plot, width = 10, height = 6)
#ggsave(paste(outdir,"Genomicdistribution_nohcc.pdf",sep=""), plot = genomic_dist_plot.nohcc, width = 10, height = 6)
#ggsave(paste(outdir,"Genomicdistribution_withoutgenes.pdf",sep=""), plot = genomic_dist_plot_small, width = 10, height = 6)

#8 Statistical testing
# Order tests: Alc vs HCC, HCC vs healthy, Healthy vs Alc
#8A Genes
aa_g <- poisson.test(x = c(distr_test[1,]$observed,distr_test[2,]$observed), T = c(distr_test[1,]$n_muts,distr_test[2,]$n_muts))$p.value #NS
bb_g <- poisson.test(x = c(distr_test[2,]$observed,distr_test[3,]$observed), T = c(distr_test[2,]$n_muts,distr_test[3,]$n_muts))$p.value #NS
cc_g <- poisson.test(x = c(distr_test[3,]$observed,distr_test[1,]$observed), T = c(distr_test[3,]$n_muts,distr_test[1,]$n_muts))$p.value #NS
#8B Promotor
aa_p <- poisson.test(x = c(distr_test[4,]$observed,distr_test[5,]$observed), T = c(distr_test[4,]$n_muts,distr_test[5,]$n_muts))$p.value #NS
bb_p <- poisson.test(x = c(distr_test[5,]$observed,distr_test[6,]$observed), T = c(distr_test[5,]$n_muts,distr_test[6,]$n_muts))$p.value #NS
cc_p <- poisson.test(x = c(distr_test[6,]$observed,distr_test[4,]$observed), T = c(distr_test[6,]$n_muts,distr_test[4,]$n_muts))$p.value # HCC vs healthy 0.0372006
#8C Enhancer
aa_e <- poisson.test(x = c(distr_test[7,]$observed,distr_test[8,]$observed), T = c(distr_test[7,]$n_muts,distr_test[8,]$n_muts))$p.value #NS
bb_e <- poisson.test(x = c(distr_test[8,]$observed,distr_test[9,]$observed), T = c(distr_test[8,]$n_muts,distr_test[9,]$n_muts))$p.value #NS
cc_e <- poisson.test(x = c(distr_test[9,]$observed,distr_test[7,]$observed), T = c(distr_test[9,]$n_muts,distr_test[7,]$n_muts))$p.value #NS
#8D Open chromatin
aa_c <- poisson.test(x = c(distr_test[10,]$observed,distr_test[11,]$observed), T = c(distr_test[10,]$n_muts,distr_test[11,]$n_muts))$p.value #NS
bb_c <- poisson.test(x = c(distr_test[11,]$observed,distr_test[12,]$observed), T = c(distr_test[11,]$n_muts,distr_test[12,]$n_muts))$p.value #NS
cc_c <- poisson.test(x = c(distr_test[12,]$observed,distr_test[10,]$observed), T = c(distr_test[12,]$n_muts,distr_test[10,]$n_muts))$p.value #NS
#8E Adjust for multiple testing
p.adjust(c(aa_g,bb_g,cc_g), method = "fdr") #0.7597268 0.7597268 0.7813692
p.adjust(c(aa_p,bb_p,cc_p), method = "fdr") #0.2006210 0.1116018 0.5318456
p.adjust(c(aa_e,bb_e,cc_e), method = "fdr") #0.6980805 0.6980805 0.4202876
p.adjust(c(aa_c,bb_c,cc_c), method = "fdr") #0.7781555 0.7781555 0.7781555
#6F Remove all p-values
remove(aa_g,bb_g,cc_g,aa_p,bb_p,cc_p,aa_e,bb_e,cc_e,aa_c,bb_c,cc_c)