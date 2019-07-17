# @Date: 11 april 2018
# @Author: Myrthe Jager
# @Modified: 7 may 2019
# @Adjusted from script of Francis Blokzijl
# @Description: Determine whether there are differences expression due to PTPRK mutation




# ---- GET STARTED ----

#1 Install & load required packages
library(ggplot2)
library(BSgenome)
library(BiocInstaller)
library(org.Hs.eg.db)
library(MutationalPatterns)
library(DESeq2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

#2 Define input and output directory
dir = ""
indir = paste(dir,"Data/",sep = "")
outdir = paste(dir,"Results/RNA/",sep = "")




# ---- GET DATA ----

#1 Get counts data
countdata <- read.table(paste(indir,"RNAseq/171208_NS500413_0383_AH7HJMBGX5_readCounts_raw.txt",sep = ""), header = T)
#...with gene names as rownames
rownames(countdata) = countdata$gene
#...remove column with genenames
countdata <- countdata[,-match("gene",names(countdata))]
#...as matrix
countdata <- as.matrix(countdata)




# ------ COUNT DISTRIBUTION PLOT ------

ggplot(melt(countdata), aes(x=Var2, y=value, fill=Var2)) +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Counts")




# ---- SAMPLE INFORMATION ----

#1 Sample name
coldata <- data.frame(SampleName = as.factor(colnames(countdata)))

#2 Sample type
coldata$Sample <- as.vector(c(rep("x",nrow(coldata))))
for(i in 1:nrow(coldata)) {
  coldata[i,]$Sample <- paste(strsplit(as.vector(coldata$SampleName[i]), split=".", fixed=TRUE)[[1]][1],strsplit(as.vector(coldata$SampleName[i]), split=".", fixed=TRUE)[[1]][2],sep=".")
}
remove(i)
coldata$Sample <- as.factor(coldata$Sample)

#3 Cell type
coldata$Celltype <- as.factor(c(rep("H",6),rep("ALC",6)))

#4 PTPRK
coldata$Type <- as.factor(c(rep("WT",6),rep("MUT",2),rep("WT",2),rep("MUT",2)))

#5 Condition
coldata$Condition <- as.factor(rep(c("noEGF","withEGF"),6))




# ------- Compare library sizes ----------

#1 Get Total number of counts
df <- melt(colSums(countdata))
#...and add info per row
df$sample <- rownames(df)
df$type <- coldata$Type
df$celltype <- coldata$Celltype
df$condition <- coldata$Condition

#2 Plot
lib_size_plot = ggplot(df, aes(y=value, x=sample, fill= sample)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Total read counts")

#3 Save plot
#ggsave(paste(outdir, "library_sizes_barplot.pdf", sep = ""), plot = lib_size_plot,width = 7, height = 4)




# ---- OVERALL CHANGES IN EXPRESSION ----

#1 Construct the DESeqDataSet object
dds = DESeqDataSetFromMatrix(countdata, coldata, design = ~ Type + Condition)

#2 Remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds) # 63677
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds) # 26270

#3 Variance stabilizing transformations --> RLD NOT VST AS WE ONLY HAVE A LIMITED NUMBER OF SAMPLES
rld <- rlog(dds)

#4 Euclidean distance
RLDdists <- dist(t(assay(rld)))
plot(hclust(RLDdists))
RLDdistMatrix <- as.matrix(RLDdists)
rownames(RLDdistMatrix) <- coldata$Type

#5 Heatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
RLD_heatmap <- pheatmap(RLDdistMatrix,
                        clustering_distance_rows = RLDdists,
                        clustering_distance_cols = RLDdists,
                        col = colors)

#5B Save plot (RLD)
#ggsave(paste(outdir, "sampledistance_heatmap.pdf", sep = ""), plot = RLD_heatmap, height = 6, width = 7)

#6 PCA plot
PCA_plot_RLD <- plotPCA(rld, intgroup = c("Type","Condition"))

#6B Save plot (RLD)
#ggsave(paste(outdir, "PCA_plot.pdf", sep = ""), plot = PCA_plot_RLD, height = 4, width = 7)




# -------- PLOT PTPRK EXPRESSION  ------

#1 Get normalized PTPRK expression
cds <- estimateSizeFactors(dds)
norm_counts <- counts(cds, normalized=TRUE)
norm_counts <- as.data.frame(norm_counts)
norm_counts$SYMBOL <- mapIds(org.Hs.eg.db,
                             keys=rownames(norm_counts),
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
norm_counts$ENSEMBL <- row.names(norm_counts)
df.ptprk = subset(norm_counts, norm_counts$SYMBOL == "PTPRK")

#2 Generate df for plotting
df.ptprk = melt(df.ptprk[,1:12])
df.ptprk$tissue = coldata$Celltype
df.ptprk$ptprk = coldata$Type
df.ptprk$condition = coldata$Condition
df.ptprk$ptprktype = paste(df.ptprk$tissue,df.ptprk$ptprk, sep ="_")
colnames(df.ptprk)[1] = "sample"
colnames(df.ptprk)[2] = "Normalized_counts"
#2B plot
plot.ptprk = ggplot(df.ptprk, aes(x= ptprktype, y= Normalized_counts, fill=ptprk)) +
  geom_boxplot() +
  geom_point() +
  theme(legend.position="none") +
  ggtitle("PTPRK") + 
  ylim(0,2000) +
  xlab("")+ 
  theme_bw() +
  theme(  axis.text.x = element_text(angle=45, hjust = 1))
#2C Save plot
#ggsave(paste(outdir,"PTPRK.pdf",sep=""),plot = plot.ptprk, width = 5, height =5)




# -------- DIFFERENTIAL EXPRESSION  ------

#1 PTPRK WT vs MUT
dds_select1 = DESeqDataSetFromMatrix(countdata, coldata, design = ~ Type)
dds_select1 = DESeq(dds_select1)
res_select1 = results(dds_select1, alpha = 0.05) 
res_select1$ENSEBML = rownames(res_select1)
# Add gene symbol
res_select1$SYMBOL = mapIds(org.Hs.eg.db,
                            keys=rownames(res_select1),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
# Get summary
summary(res_select1)
# UP
resSig_select1_up = subset(res_select1, padj < 0.05 & log2FoldChange > 0)
resSig_select1_up = resSig_select1_up[order(resSig_select1_up$log2FoldChange, decreasing = T),]
# DOWN
resSig_select1_down = subset(res_select1, padj < 0.05 & log2FoldChange < 0)
resSig_select1_down = resSig_select1_down[order(resSig_select1_down$log2FoldChange, decreasing = F),]
#write.table(res_select1, file = paste(outdir, "PTPRK_all.txt", sep = ""), sep = "\t", quote = F, row.names = F)
#write.table(resSig_select1_up, file = paste(outdir, "PTPRK_sig_up.txt", sep = ""), sep = "\t", quote = F, row.names = F)
#write.table(resSig_select1_down, file = paste(outdir, "PTPRK_sig_down.txt", sep = ""), sep = "\t", quote = F, row.names = F)

#2 no EGF vs with EGF
#2A all
  dds_select2 = DESeqDataSetFromMatrix(countdata, coldata, design = ~ Condition)
  dds_select2 = DESeq(dds_select2)
  res_select2 = results(dds_select2, alpha = 0.05) 
  res_select2$ENSEBML = rownames(res_select2)
  # Add gene symbol
  res_select2$SYMBOL = mapIds(org.Hs.eg.db,
                            keys=rownames(res_select2),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
  # Get summary
  summary(res_select2)
  # No differences
#2B per sample type
  select3 = which(coldata$Type == "MUT") # --> run code below --> no diff genes
  select3 = which(coldata$Celltype == "H") # --> run code below --> no diff genes
  select3 = which(coldata$Celltype == "ALC") # --> run code below --> no diff genes
  select3 = which(coldata$Type == "WT") # --> run code below --> hardly diff genes: 3
  dds_select3 = DESeqDataSetFromMatrix(countdata[,select3], coldata[select3,], design = ~ Condition)
  dds_select3 = DESeq(dds_select3)
  res_select3 = results(dds_select3, alpha = 0.05) 
  res_select3$ENSEBML = rownames(res_select3)
  # Add gene symbol
  res_select3$SYMBOL = mapIds(org.Hs.eg.db,
                            keys=rownames(res_select3),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
  # Get summary
  summary(res_select3)
  resSig_select3 = subset(res_select3, padj < 0.05)
  as.matrix(resSig_select3 )

#3 ALC vs healthy
dds_select4 = DESeqDataSetFromMatrix(countdata, coldata, design = ~ Celltype)
dds_select4 = DESeq(dds_select4)
res_select4 = results(dds_select4, alpha = 0.05) 
res_select4$ENSEBML = rownames(res_select4)
# Add gene symbol
res_select4$SYMBOL = mapIds(org.Hs.eg.db,
                            keys=rownames(res_select4),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
# Get summary
summary(res_select4)
# UP
resSig_select4_up = subset(res_select4, padj < 0.05 & log2FoldChange > 0)
resSig_select4_up = resSig_select4_up[order(resSig_select4_up$log2FoldChange, decreasing = T),]
# DOWN
resSig_select4_down = subset(res_select4, padj < 0.05 & log2FoldChange < 0)
resSig_select4_down = resSig_select4_down[order(resSig_select4_down$log2FoldChange, decreasing = F),]
#write.table(res_select4, file = paste(outdir, "ALC_all.txt", sep = ""), sep = "\t", quote = F, row.names = F)
#write.table(resSig_select4_up, file = paste(outdir, "ALC_sig_up.txt", sep = ""), sep = "\t", quote = F, row.names = F)
#write.table(resSig_select4_down, file = paste(outdir, "ALC_sig_down.txt", sep = ""), sep = "\t", quote = F, row.names = F)