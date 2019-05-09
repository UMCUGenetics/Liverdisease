library("SPARQL")
library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(ggsignif)

# registerDoParallel(cores=20)

# Arguments for the use of the graph database
endpoint <- "https://sparqling-genomics.op.umcutrecht.nl/sparql-auth"
auth_options <- curlOptions(userpwd=<CREDENTIALS FOR SPARQL>)

# Query the graph database for all variants in the livers
query <- "
SELECT DISTINCT ?var
FROM <http://alcoholic_livers/somaticvariant> 
WHERE {
  ?var <http://sparqling-genomics/sample> ?sample
}"
samples <- SPARQL(endpoint, query, curl_args = auth_options)$results

# Where to write the mutations in peaks
filename <- <FILENAME> 
write.table("Sample Chromosome Position REF ALT",
            file = filename,
            row.names = F, col.names = F, quote = F)

# For each variant, check on overlap with peak 
# Give sample, chromosome, position, ref and alt if variant in peak
for (i in 1:length(samples)){
  print(i)
  query <- sprintf("
    PREFIX faldo:<http://biohackathon.org/resource/faldo#>
    PREFIX sq:<http://sparqling-genomics/>
    PREFIX Sample:<http://sparqling-genomics/Sample/>
    PREFIX col:<http://sparqling-genomics/table2rdf/Column/>
    PREFIX hg19:<http://rdf.biosemantics.org/data/genomeassemblies/hg19#>
    PREFIX variantcall:<http://sparqling-genomics/vcf2rdf/VariantCall/>
    PREFIX sequence:<http://sparqling-genomics/vcf2rdf/Sequence/>
    
    SELECT STRAFTER(STR(?var_sample), STR(Sample:)) AS ?sample
    STRAFTER(STR(?var_chr), STR(hg19:)) AS ?chromosome
    ?var_pos AS ?position
    STRAFTER(STR(?var_ref), STR(sequence:)) AS ?REF
    STRAFTER(STR(?var_alt), STR(sequence:)) AS ?ALT 
    
    WHERE {
      GRAPH <http://alcoholic_livers/somaticvariant> {
        %s  faldo:reference   ?var_chr ;
              faldo:position    ?var_pos ;
              variantcall:REF   ?var_ref ;
              variantcall:ALT   ?var_alt ;
              sq:sample         ?var_sample .
      }
    
      GRAPH <http://alcoholic_livers/H3K36Me3Peaks> {
        ?peak col:chromosome  ?peak_chr ;
              col:start       ?peak_start ;
              col:end         ?peak_end .
      }

      FILTER( ?var_chr = ?peak_chr )
      FILTER( ?var_pos >= ?peak_start && ?var_pos <= ?peak_end )
    }
    ", samples[i])

  df <- SPARQL(endpoint, query, curl_args = auth_options)$results
  
  write.table(df,
              file = filename,
              row.names = F, col.names = F, quote = F,
              append = T)
}

# Read in found results with above query
counts <- read.table(<FILENAME>, header=T)

# Filter mutations to look only at T>C mutations
counts <- counts[(counts$REF=="A" & counts$ALT=="G") | (counts$REF=="T" & counts$ALT=="C"),]

# For each sample query the total number of variants
query <- "
PREFIX sample:<http://sparqling-genomics/Sample/>

SELECT STRAFTER(STR(?sample), STR(sample:)) AS ?Sample
COUNT(?var) AS ?numberOfSNPs
FROM <http://alcoholic_livers/somaticvariant>
WHERE {
  ?var <http://sparqling-genomics/sample> ?sample
}
ORDER BY ?sample
"
total_counts <- SPARQL(endpoint, query, curl_args = auth_options)$results

# Categorize the samples in healthy or diseased
healthy <- unique(counts$Sample)[c(1:10,12,17,20,24,25)]
decease <- unique(counts$Sample)[-c(1:10,12,17,20,24,25)]

# Count the number of variants in peaks for each sample
# Get absolute and relative number of SNPs in peaks and total number of SNPs
counts_samples <- data.frame(counts %>% count(Sample))
rel_counts <- foreach(i = 1:nrow(counts_samples), .combine = "rbind") %do% {
  rownumb <- which(total_counts$Sample == counts_samples$Sample[i])
  if (counts_samples$Sample[i] %in% healthy){
    type <- "Healthy"
  } else {
    type <- "Alcoholic"
  }

  return(data.frame(Sample = counts_samples$Sample[i],
                    Type = type,
                    NumberOfSNPsInPeaks = counts_samples$n[i],
                    RelativeNumberOfSNPsInPeaks = counts_samples$n[i] / total_counts$numberOfSNPs[rownumb],
                    TotalNumberOfSNPs = total_counts$numberOfSNPs[rownumb]))
}

healthy_counts <- sum(rel_counts[rel_counts$Sample %in% healthy, "rel_n"])
disease_counts <- sum(rel_counts[rel_counts$Sample %in% decease, "rel_n"])

# Plot results in boxplots per type (healthy, disease)
p <- ggplot(rel_counts, aes(x=Type, y = RelativeNumberOfSNPsInPeaks)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("Healthy", "Alcoholic")), textsize=6) +
  labs(title = "Relative number of T>C mutations in H3K36Me3 histones of healthy and alcoholic livers", 
       y="Relative number of T>C in peaks") +
  theme(text = element_text(size=20))

# Save the plot
pdf(<PLOTNAME>, width=20, height=15)
print(p)
dev.off()
