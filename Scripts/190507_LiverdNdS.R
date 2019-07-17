# @Date: 11 april 2018
# @Author: Myrthe Jager
# @Modified: 7 may 2019
# @Description: Determine whether there are differences in dnds & number of coding mutations (in cosmic census genes)




# ---- GET STARTED ----

#1 Load required packages
library("seqinr")
library("Biostrings")
library(GenomicRanges)
library("MASS")
library("GenomicRanges")
#install_github("im3sanger/dndscv")
library("dndscv")
library(VariantAnnotation)
library(ggplot2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)

#2 Define input and output directory
dir = ""
indir = paste(dir,"Data/",sep = "")
outdir = paste(dir,"Results/SNV/coding/",sep = "")

#3 Functions
## Function 1 : Convert vcf to input required for dnds script ##
prepare_for_dnds = function(vcflist) {
  df_new <- data.frame(  )
  for (i in vcflist){
    sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[1]
    muts <- NULL
    vcf <- readVcf(i, "hg19")
    seqlevels(vcf) <- paste("chr", seqlevels(vcf), sep="")
    df.sample <- data.frame(sampleID = sample,
                            chr = as.numeric(seqnames(vcf)),
                            pos = start(vcf),
                            ref = as.vector(rowRanges(vcf)$REF),
                            mut = as.vector(unlist(rowRanges(vcf)$ALT)))
    df_new <- rbind(df_new, df.sample)
  }
  return(df_new)
}

## Function 2: Convert table to Granges ##
make_granges = function(df) {
  df_new <- GRanges(seqnames = paste("chr", df$chr, sep = ""),
                    ranges = IRanges(start = df$start,
                                     end = df$stop),
                    strand = '*',
                    expressionCat = df$cat,
                    gene_name = df$name)
  return(df_new)
}

#4 Metadata
metadata <- read.delim(paste(indir,"metadata.txt",sep=""), sep = "\t",header=TRUE)




# ---- dN/dS ----

#1 provide location to vcfs
alc_vcf_files <- list.files(paste(indir,"SNV/ALC/IAP/",sep=""),full.names = T)
healthyliver_vcf_files <- list.files(paste(indir,"SNV/Healthy_Liver/",sep=""), full.names = T)
hcc_clonal_vcf_file <- list.files(paste(indir,"SNV/HCC/clonal/",sep=""), full.names = T)

#2 Convert vcf to required input for dnds
alc.dnds.input <- prepare_for_dnds(alc_vcf_files)
healthy.dnds.input <- prepare_for_dnds(healthyliver_vcf_files)
hcc.dnds.input <- prepare_for_dnds(hcc_clonal_vcf_file)

#3 dnds
alc.dnds <- dndscv(alc.dnds.input)
healthy.dnds <- dndscv(healthy.dnds.input)
hcc.dnds <- dndscv(hcc.dnds.input)

#4 view dnds
alc.dnds$globaldnds
healthy.dnds$globaldnds
hcc.dnds$globaldnds

#5 Plot dnds
#5A create dataframe
df.dnds <- data.frame(sample = c("Alc","healthy","HCC"),
                      mle = c(alc.dnds$globaldnds[5,2],healthy.dnds$globaldnds[5,2],hcc.dnds$globaldnds[5,2]),
                      cilow = c(alc.dnds$globaldnds[5,3],healthy.dnds$globaldnds[5,3],hcc.dnds$globaldnds[5,3]),
                      cihigh = c(alc.dnds$globaldnds[5,4],healthy.dnds$globaldnds[5,4],hcc.dnds$globaldnds[5,4])
)
#5B Plot
dnds.plot <- ggplot(df.dnds, aes(x=sample, y=mle)) +
  geom_errorbar(width=.1, aes(ymin=cilow, ymax=cihigh)) +
  geom_point(shape=21, size=3, fill="white") +
  ylim(0,2) +
  xlab("")+
  ylab("dN/dS") +
  theme_bw()
#5C Save plot
#ggsave(paste(outdir,"dnds.pdf",sep=""),plot=dnds.plot,height = 4,width = 4)




# ---- MUTATIONS IN COSMIC GENES ----

#1 Get genes with any mutations
#1A ALC
alc.hcc.all <- alc.dnds$sel_cv
alc.hcc.all$sum <- rowSums(alc.dnds$sel_cv[,2:5])
alc.hcc.all <- alc.hcc.all[which(alc.hcc.all$sum > 0),]
#1B Healthy
healthy.hcc.all <- healthy.dnds$sel_cv
healthy.hcc.all$sum <- rowSums(healthy.dnds$sel_cv[,2:5])
healthy.hcc.all <- healthy.hcc.all[which(healthy.hcc.all$sum > 0),]
#1C HCC
hcc.hcc.all <- hcc.dnds$sel_cv
hcc.hcc.all$sum <- rowSums(hcc.dnds$sel_cv[,2:5])
hcc.hcc.all <- hcc.hcc.all[which(hcc.hcc.all$sum > 0),]

#2 Retrieve only relevant mutations (nonsyn, nonsense)
alc.hcc.all.ns <- alc.hcc.all[which(alc.hcc.all$n_mis+alc.hcc.all$n_non+alc.hcc.all$n_spl > 0),]
healthy.hcc.all.ns <- healthy.hcc.all[which(healthy.hcc.all$n_mis+healthy.hcc.all$n_non+healthy.hcc.all$n_spl > 0),]
hcc.hcc.all.ns <- hcc.hcc.all[which(hcc.hcc.all$n_mis+hcc.hcc.all$n_non+hcc.hcc.all$n_spl > 0),]

#3 Get cosmic genes list from https://cancer.sanger.ac.uk/census
cosmic.genes <- read.csv(paste(indir,"Census_all_cosmic.csv",sep=""),header=TRUE)
cosmic.genes <- cosmic.genes[which(cosmic.genes$Tier ==1),]
cosmic.genes <- cosmic.genes[which(cosmic.genes$Molecular.Genetics == "Dom"),]

#4 Create empty dataframe
cosmic.genes$tissue.type <- rep(NA,nrow(cosmic.genes))
cosmic.genes.mutations <- cosmic.genes[0,]
cosmic.genes <- cosmic.genes[,-21]

#5 Get cosmic genes with mutation in any sample, and add sample name
for(i in 1:nrow(alc.hcc.all.ns)) {
  mut_gene <- as.character(alc.hcc.all.ns[i,]$gene_name)
  temp.df <- cosmic.genes[which(cosmic.genes$Gene == mut_gene),]
  if(nrow(temp.df) == 1) {
    temp.df$tissue.type <- "ALC"
    temp.df$sample <- NA
    if (length(alc.dnds$annotmuts[which(alc.dnds$annotmuts[,6] == mut_gene),1]) == 1 ) {
      temp.df$sample <- alc.dnds$annotmuts[which(alc.dnds$annotmuts[,6] == mut_gene),]$sampleID
    }
    if (length(alc.dnds$annotmuts[which(alc.dnds$annotmuts[,6] == mut_gene),1]) > 1 ) {
      temp.df$sample <- paste(alc.dnds$annotmuts[which(alc.dnds$annotmuts[,6] == mut_gene),]$sampleID, collapse = ",")
    }
    cosmic.genes.mutations <- rbind(cosmic.genes.mutations,temp.df)
  }
  remove(mut_gene,temp.df)
}
remove(i)
for(i in 1:nrow(healthy.hcc.all.ns)) {
  mut_gene <- as.character(healthy.hcc.all.ns[i,]$gene_name)
  temp.df <- cosmic.genes[which(cosmic.genes$Gene == mut_gene),]
  if(nrow(temp.df) == 1) {
    temp.df$tissue.type <- "Healthy"
    temp.df$sample <- NA
    if (length(healthy.dnds$annotmuts[which(healthy.dnds$annotmuts[,6] == mut_gene),1]) == 1 ) {
      temp.df$sample <- healthy.dnds$annotmuts[which(healthy.dnds$annotmuts[,6] == mut_gene),]$sampleID
    }
    if (length(healthy.dnds$annotmuts[which(healthy.dnds$annotmuts[,6] == mut_gene),1]) > 1 ) {
      temp.df$sample <- paste(healthy.dnds$annotmuts[which(healthy.dnds$annotmuts[,6] == mut_gene),]$sampleID, collapse = ",")
    }
    cosmic.genes.mutations <- rbind(cosmic.genes.mutations,temp.df)
  }
  remove(mut_gene,temp.df)
}
remove(i)
for(i in 1:nrow(hcc.hcc.all.ns)) {
  mut_gene <- as.character(hcc.hcc.all.ns[i,]$gene_name)
  temp.df <- cosmic.genes[which(cosmic.genes$Gene == mut_gene),]
  if(nrow(temp.df) == 1) {
    temp.df$tissue.type <- "HCC"
    temp.df$sample <- "HCC1-clonal"
    cosmic.genes.mutations <- rbind(cosmic.genes.mutations,temp.df)
  }
  remove(mut_gene,temp.df)
}
remove(i)

#5 Retrieve the p-values and q values
#5A Empty table
alc.hcc.all.cosmic <- alc.hcc.all.ns[0,]
healthy.hcc.all.cosmic <- healthy.hcc.all.ns[0,]
hcc.hcc.all.cosmic <- hcc.hcc.all.ns[0,]
#5B Add lines with cosmic genes
for(i in 1:nrow(cosmic.genes)) {
  mut_gene <- as.character(cosmic.genes[i,]$Gene.Symbol)
  temp.df <- alc.hcc.all[which(alc.hcc.all$gene_name == mut_gene),]
  if(nrow(temp.df) == 1) {
    alc.hcc.all.cosmic <- rbind(alc.hcc.all.cosmic,temp.df)
  }
  remove(temp.df)
  temp.df <- healthy.hcc.all[which(healthy.hcc.all$gene_name == mut_gene),]
  if(nrow(temp.df) == 1) {
    healthy.hcc.all.cosmic <- rbind(healthy.hcc.all.cosmic,temp.df)
  }
  remove(temp.df)
  temp.df <- hcc.hcc.all[which(hcc.hcc.all$gene_name == mut_gene),]
  if(nrow(temp.df) == 1) {
    hcc.hcc.all.cosmic <- rbind(hcc.hcc.all.cosmic,temp.df)
  }
  remove(temp.df, mut_gene)
}

#6 Combine into one table
#6A first add rownames
row.names(cosmic.genes.mutations) <- cosmic.genes.mutations$Gene.Symbol
row.names(alc.hcc.all.cosmic) <- alc.hcc.all.cosmic$gene_name
row.names(healthy.hcc.all.cosmic) <- healthy.hcc.all.cosmic$gene_name
row.names(hcc.hcc.all.cosmic) <- hcc.hcc.all.cosmic$gene_name
#6B then merge by rownames
final.df <- merge(cosmic.genes.mutations,
                  rbind(alc.hcc.all.cosmic[,-1], healthy.hcc.all.cosmic[which(healthy.hcc.all.cosmic$n_syn ==0),-1],hcc.hcc.all.cosmic[,-1]), by= "row.names")

#7 Order by sample type
final.df <- final.df[order(final.df$tissue.type),]

#8 Save table
#write.table(final.df,paste(outdir,"cosmic_all.txt",sep=""),sep="\t",row.names = FALSE,col.names = TRUE)





# ---- Plot coding mutations vs age ----

#1 Load surveyed file
snvcorrected.coding <- read.delim(paste(indir,"/callable/autosomal_callable_livers.txt",sep=""),header=T, sep="\t")
snvcorrected.coding$Surveyed.percentage <- as.numeric(gsub("%", "", unlist(snvcorrected.coding$Surveyed.percentage)))
snvcorrected.coding[24,]$Surveyed.percentage <- 100

#2 Add additional info
snvcorrected.coding$sample_type <- c(rep("healthy",15),rep("alc",8),"hcc")
snvcorrected.coding$donor <- NA
for (i in 1:nrow(snvcorrected.coding)) {
  snvcorrected.coding[i,]$donor <- metadata[which(metadata$sample_id == as.character(snvcorrected.coding[i,]$Sample)),]$donor_id
}
remove(i)
snvcorrected.coding$age <- NA
for (i in 1:nrow(snvcorrected.coding)) {
  snvcorrected.coding[i,]$age <- metadata[which(metadata$sample_id == as.character(snvcorrected.coding[i,]$Sample)),]$age
}
remove(i)

#3 Add coding mutations
snvcorrected.coding$uncorr.coding <- NA
#3A Calculate coding mutations Healthy
for (i in 1:length(unique(healthy.dnds$annotmuts[,1]))) {
  sample <- unique(healthy.dnds$annotmuts[,1])[i]
  sample_id <- as.character(metadata[which(metadata$sample_name == sample),]$sample_id)
  snvcorrected.coding[which(snvcorrected.coding$Sample == sample_id),]$uncorr.coding <- nrow(healthy.dnds$annotmuts[which(healthy.dnds$annotmuts[,1] == sample),])
  remove(sample, sample_id)
}
remove(i)
#3B Calculate coding mutations ALC
for (i in 1:length(unique(alc.dnds$annotmuts[,1]))) {
  sample <- unique(alc.dnds$annotmuts[,1])[i]
  sample_id <- as.character(metadata[which(metadata$sample_name == sample),]$sample_id)
  snvcorrected.coding[which(snvcorrected.coding$Sample == sample_id),]$uncorr.coding <- nrow(alc.dnds$annotmuts[which(alc.dnds$annotmuts[,1] == sample),])
  remove(sample, sample_id)
}
remove(i)
#3C Calculate coding mutations HCC
for (i in 1:length(unique(hcc.dnds$annotmuts[,1]))) {
  sample <- unique(hcc.dnds$annotmuts[,1])[i]
  sample_id <- as.character(metadata[which(metadata$sample_name == sample),]$sample_id)
  snvcorrected.coding[which(snvcorrected.coding$Sample == sample_id),]$uncorr.coding <- nrow(hcc.dnds$annotmuts[which(hcc.dnds$annotmuts[,1] == sample),])
  remove(sample, sample_id)
}
remove(i)

#4 correct coding number
snvcorrected.coding$corr.coding <- (snvcorrected.coding$uncorr.coding/snvcorrected.coding$Surveyed.percentage)*100
snvcorrected.coding[which(snvcorrected.coding$Sample == "HCC1-clonal"),]$corr.coding = snvcorrected.coding[which(snvcorrected.coding$Sample == "HCC1-clonal"),]$uncorr.coding

#5 Plot
snv.plot.combined.coding <- ggplot(snvcorrected.coding, aes(x=age, y=corr.coding, color=sample_type)) +
  geom_point() +
  geom_smooth(method=lm, aes(fill=sample_type)) +
  expand_limits(x=c(0,100)) +
  ylab("No. coding base substitutions per autosomal genome") +
  #  facet_wrap( ~ sample_type) +
  theme_bw() +
  xlab("Age (years)")
#ggsave(paste(outdir,"coding_mutationload_withhcc.pdf",sep=""),plot = snv.plot.combined.coding,height = 4, width = 6)

#6 Plot without HCC
snv.plot.combined.coding.nohcc <- ggplot(snvcorrected.coding[which(snvcorrected.coding$sample_type != "hcc"),], aes(x=age, y=corr.coding, color=sample_type)) +
  geom_point() +
  geom_smooth(method=lm, aes(fill=sample_type)) +
  expand_limits(x=c(0,100)) +
  ylab("No. coding base substitutions per autosomal genome") +
  #  facet_wrap( ~ sample_type) +
  theme_bw() +
  xlab("Age (years)")
#ggsave(paste(outdir,"coding_mutationload.pdf",sep=""),plot = snv.plot.combined.coding.nohcc,height = 4, width = 6)




# ---- Plot number of cosmic mutations ----

#1 Add cosmic mutations to table
#1A 
df.cosmic.number <- data.frame(table(unlist(strsplit(final.df$sample, split =","))))
#1B Add sample_id (to enable adding number of coding mutations to the snvcorrected.coding table)
df.cosmic.number$Sample <- NA
for (i in 1:nrow(df.cosmic.number)) {
  df.cosmic.number[i,]$Sample <- as.character(metadata[which(metadata$sample_name == as.character(df.cosmic.number[i,]$Var1)),]$sample_id)
}
remove(i)
#1C Add numbers to snvcorrected.coding table
snvcorrected.coding$uncorr.cosmic <- rep(0,nrow(snvcorrected.coding))
for (i in 1:nrow(snvcorrected.coding)) {
  temp <- df.cosmic.number[which(df.cosmic.number$Sample == snvcorrected.coding[i,]$Sample),]
  if(nrow(temp) == 1) {
    snvcorrected.coding[i,]$uncorr.cosmic <- df.cosmic.number[which(df.cosmic.number$Sample == snvcorrected.coding[i,]$Sample),]$Freq
  }
  remove(temp)
}
remove(i)

#2 correct cosmic number
snvcorrected.coding$corr.cosmic <- (snvcorrected.coding$uncorr.cosmic/snvcorrected.coding$Surveyed.percentage)*100

#3 Plot uncorrected
persamplecosmic <- ggplot(snvcorrected.coding, aes(x=Sample, y=uncorr.cosmic, fill = sample_type)) +
  geom_bar(stat="identity", aes(fill = sample_type)) +
  #  expand_limits(x=c(0,100)) +
  ylab("No. nonsense/nonsynonymous base substitutions\nin cancer genes") +
  #  facet_wrap( ~ sample_type) +
  theme_bw() +
  xlab("Sample") +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle=90))
#ggsave(paste(outdir,"cosmic_persample.pdf",sep = ""),plot = persamplecosmic, width = 4, height = 4)

#4 Save table
#write.table(snvcorrected.coding,paste(outdir,"coding_cosmicnumber.txt",sep=""),sep="\t",row.names = FALSE,col.names = TRUE)