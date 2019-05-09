# @Date: 06 April 2018
# @Author: Myrthe Jager
# @Modified: 
# @Description: Get chromosome lengths of hg19 genome




#1 Define outdir
dir = ""
outdir = paste(dir,"Data/callable/",sep = "")

#2 Load human genome hg19
#2A Install & load BSgenome package
# biocLite("BSgenome")
library(BSgenome)
#2B select genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#2C Download & load genome
# biocLite(ref_genome)
library(ref_genome, character.only = TRUE)

#3 Autosomal seqlengths
#3A Calculate
autosomal_seqlengths = seqlengths(get(ref_genome))[1:22]
#3B Save
#write.table(autosomal_seqlengths, paste(outdir,"autosomal_chr_lengths_hg19.txt",sep=""), quote = F, sep = "\t", col.names=F)