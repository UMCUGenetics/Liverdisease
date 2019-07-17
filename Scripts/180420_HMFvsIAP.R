# @Date: 05 april 2018
# @Author: Myrthe Jager
# @Modified: 20 april 2018
# @Description: Genome-wide mutation pattern analysis alcoholic IAP vs alcoholic HMF
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
library(VennDiagram)

#2 Define input and output directory
dir = ""
indir = paste(dir,"Data/",sep = "")
outdir = paste(dir,"Results/SNV/",sep = "")

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

#4 Install and load mouse reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#biocLite(ref_genome)
library(ref_genome, character.only = T)

#5 Metadata
metadata <- read.delim(paste(indir,"metadata.txt",sep=""), sep = "\t",header=TRUE)




# ---- GET VCFS ----

#1 ALC liver IAP files + sample names
alc_vcf_files <- list.files(paste(indir,"SNV/ALC/IAP/",sep=""),full.names = T)
alc_sample_names <-get_sampleinfo(alc_vcf_files,"sample_id")

#2 ALC liver HMF files (sample names is similar as alc_sample_names)
alc_vcf_files_HMF <- list.files(paste(indir,"SNV/ALC/HMF/",sep=""),full.names = T)

#3 Read combined
alc_sample_names_HMF<- c(paste(alc_sample_names,"_IAP",sep=""),paste(alc_sample_names,"_HMF",sep=""))
alc_vcfs_HMF <- read_vcfs_as_granges(c(alc_vcf_files,alc_vcf_files_HMF),alc_sample_names_HMF,genome = ref_genome)




# ---- GET & ORDER SAMPLE INFO ----

#1 Sample type
type <- factor(c(rep("Alc_IAP",8),rep("Alc_HMF",8)), levels = c("Alc_IAP","Alc_HMF"))




# ------ MUTATION NUMBER & RATE ------

#1 Load surveyed file
snvcorrected <- read.delim(paste(indir,"/callable/autosomal_callable_livers.txt",sep=""),header=T, sep="\t")
snvcorrected$Surveyed.percentage <- as.numeric(gsub("%", "", unlist(snvcorrected$Surveyed.percentage)))

#2 Add info
snvcorrected$donor <- NA
for (i in 1:nrow(snvcorrected)) {
  snvcorrected[i,]$donor <- metadata[which(metadata$sample_id == as.character(snvcorrected[i,]$Sample)),]$donor_id
}
remove(i)
snvcorrected$age <- NA
for (i in 1:nrow(snvcorrected)) {
  snvcorrected[i,]$age <- metadata[which(metadata$sample_id == as.character(snvcorrected[i,]$Sample)),]$age
}
remove(i)

#3 Only Alc
snvcorrected_HMF <- snvcorrected[grep("alc",snvcorrected$Sample),]

#4 Add mutations
#4A IAP uncorrected
snvcorrected_HMF$uncorr.snv.IAP <- NA
for (i in alc_vcf_files) {
  df <- read.table(i)
  sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[1]
  snv_number <- nrow(df)
  sampleid <- as.character(metadata[which(metadata$sample_name == sample),]$sample_id)

  snvcorrected_HMF[which(snvcorrected_HMF$Sample == sampleid),]$uncorr.snv.IAP <- snv_number
  
  remove(sampleid,snv_number,sample,df)
}
remove(i)
#4B IAP corrected
snvcorrected_HMF$corr.snv.IAP <- round(snvcorrected_HMF$uncorr.snv.IAP/(snvcorrected_HMF$Surveyed.percentage/100),0)
#4C HMF uncorrected
snvcorrected_HMF$uncorr.snv.HMF <- NA
for (i in alc_vcf_files_HMF) {
  df <- read.table(i)
  sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[2]
  snv_number <- nrow(df)
  sampleid <- as.character(metadata[which(metadata$sample_name == sample),]$sample_id)
  
  snvcorrected_HMF[which(snvcorrected_HMF$Sample == sampleid),]$uncorr.snv.HMF <- snv_number
  
  remove(sampleid,snv_number,sample,df)
}
remove(i)

#5 Calculate increase/decrease
snvcorrected_HMF$difference <- (snvcorrected_HMF$uncorr.snv.HMF-snvcorrected_HMF$corr.snv.IAP)*-1

#6 calculate increase fraction
snvcorrected_HMF$IAP_to_HMF <- paste("x",round((snvcorrected_HMF$difference/snvcorrected_HMF$corr.snv.IAP),3)+1)
snvcorrected_HMF$HMF_to_IAP <- paste("x",1-round((snvcorrected_HMF$difference/snvcorrected_HMF$uncorr.snv.HMF)*1,3))

#7 Save numbers table
#write.table(snvcorrected_HMF,file = paste(outdir,"/IAP_vs_HMF/mutationrate_pipelines.txt",sep=""), sep = "\t", row.names = FALSE)

#8 Venn diagrams
grid.newpage()
draw.pairwise.venn(as.numeric(snvcorrected_HMF$uncorr.snv.IAP[1]),as.numeric(snvcorrected_HMF$uncorr.snv.HMF[1]), 2319, 
                   category = c("alc5-a IAP", "alc5-a HMF"), lty = rep("blank", 2), fill = c("light blue", "pink"))
grid.newpage()
draw.pairwise.venn(as.numeric(snvcorrected_HMF$uncorr.snv.IAP[2]),as.numeric(snvcorrected_HMF$uncorr.snv.HMF[2]), 2251, 
                   category = c("alc5-b IAP", "alc5-b HMF"), lty = rep("blank", 2), fill = c("light blue", "pink"))
grid.newpage()
draw.pairwise.venn(as.numeric(snvcorrected_HMF$uncorr.snv.IAP[3]),as.numeric(snvcorrected_HMF$uncorr.snv.HMF[3]), 1728, 
                   category = c("alc1-a IAP", "alc1-a HMF"), lty = rep("blank", 2), fill = c("light blue", "pink"))
grid.newpage()
draw.pairwise.venn(as.numeric(snvcorrected_HMF$uncorr.snv.IAP[4]),as.numeric(snvcorrected_HMF$uncorr.snv.HMF[4]), 2259, 
                   category = c("alc4-a IAP", "alc4-a HMF"), lty = rep("blank", 2), fill = c("light blue", "pink"))
grid.newpage()
draw.pairwise.venn(as.numeric(snvcorrected_HMF$uncorr.snv.IAP[5]),as.numeric(snvcorrected_HMF$uncorr.snv.HMF[5]), 1896, 
                   category = c("alc4-b IAP", "alc4-b HMF"), lty = rep("blank", 2), fill = c("light blue", "pink"))
grid.newpage()
draw.pairwise.venn(as.numeric(snvcorrected_HMF$uncorr.snv.IAP[6]),as.numeric(snvcorrected_HMF$uncorr.snv.HMF[6]), 1888, 
                   category = c("alc2-a IAP", "alc2-a HMF"), lty = rep("blank", 2), fill = c("light blue", "pink"))
grid.newpage()
draw.pairwise.venn(as.numeric(snvcorrected_HMF$uncorr.snv.IAP[7]),as.numeric(snvcorrected_HMF$uncorr.snv.HMF[7]), 2398, 
                   category = c("alc3-a IAP", "alc3-a HMF"), lty = rep("blank", 2), fill = c("light blue", "pink"))
grid.newpage()
draw.pairwise.venn(as.numeric(snvcorrected_HMF$uncorr.snv.IAP[8]),as.numeric(snvcorrected_HMF$uncorr.snv.HMF[8]), 2570, 
                   category = c("alc3-b IAP", "alc3-b HMF"), lty = rep("blank", 2), fill = c("light blue", "pink"))




# ---- MUTATION SPECTRA & PROFILES ----

#1 Get 6 types
type_occurrences_HMF <- mut_type_occurrences(alc_vcfs_HMF,ref_genome)

#2 PLot 6 types
plot_spectrum(type_occurrences_HMF, by = type, legend = T, CT = T)
#ggsave(paste(outdir,"/IAP_vs_HMF/spectrum_IAPvsHMF.pdf",sep=""),plot = plot_spectrum(type_occurrences_HMF, by = type, legend = T, CT = T), width = 8, height = 4)

#3 Get 96 types
mut_matrix_HMF <- mut_matrix(vcf_list = alc_vcfs_HMF, ref_genome = ref_genome)

#4 Plot 96 types
plot_96_profile(mut_matrix_HMF, ymax = 0.1, condensed = T)
#ggsave(paste(outdir,"/IAP_vs_HMF/profile_IAPvsHMF_persample.pdf",sep=""),plot = plot_96_profile(mut_matrix_HMF, ymax = 0.1, condensed = T), height = 20, width = 8)

#5 Collapse 96 mutation types matrix
collapsed_mutmatrix_HMF <- data.frame(rowSums(mut_matrix_HMF[,c(1:8)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix_HMF[,c(9:16)], na.rm = FALSE, dims = 1)
)
colnames(collapsed_mutmatrix_HMF) <- c("IAP","HMF")
# Plot
plot_96_profile(collapsed_mutmatrix_HMF, ymax = 0.05, condensed = T)
#ggsave(paste(outdir,"/IAP_vs_HMF/profile_IAPvsHMF.pdf",sep=""),plot =plot_96_profile(collapsed_mutmatrix_HMF, ymax = 0.05, condensed = T), height = 6, width = 8)

#6 Cosine similarity
cosine.collapsed.HMF <- cos_sim_matrix(collapsed_mutmatrix_HMF,collapsed_mutmatrix_HMF)
cosine.collapsed.HMF.all <- cos_sim_matrix(mut_matrix_HMF[,c(1:8)],mut_matrix_HMF[,c(9:16)])

#7 Plot
cosine.plot.HMF <- plot_cosine_heatmap(cosine.collapsed.HMF,plot_values = T)
cosine.plot.HMF.all <- plot_cosine_heatmap(cosine.collapsed.HMF.all,plot_values = T,cluster_rows = F)
#ggsave(paste(outdir,"/IAP_vs_HMF/cosinesimilarity_pipelines.pdf",sep=""),plot = cosine.plot.HMF, width = 4, height = 3)
#ggsave(paste(outdir,"/IAP_vs_HMF/cosinesimilarity_pipelines_persample.pdf",sep=""),plot = cosine.plot.HMF.all, width = 6, height = 6)

#8 Plot compared profiles
compare.plot <- plot_compare_profiles(collapsed_mutmatrix_HMF[,1],collapsed_mutmatrix_HMF[,2], condensed = TRUE,profile_ymax = 0.05, profile_names = colnames(collapsed_mutmatrix_HMF[1:2]))
#ggsave(paste(outdir,"/IAP_vs_HMF/compared_profiles_HMF_IAP.pdf",sep=""),plot = compare.plot, width = 6, height = 6)