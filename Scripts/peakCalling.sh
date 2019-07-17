#!/bin/bash
#
# Script need as input the bedGraph file with the union of all obtained
# bigWig files converted by combineSignals.sh

FILEPATH=$1
FILE=$( echo $FILEPATH | cut -d'.' -f 1 )

echo "MERGING H3k36me3 SCORES OF ${FILE}.bedGraph"

# For each genome region, take the median H3k36me3 SCORE of each SAMPLES
awk '{n=split($0,a);
      for(i = 4; i <= n; ++i) b[i-3]=a[i];
      asort(b);
      median=int((n+1)/2);
      if (median < (n+1)/2) print $1,"\t",$2,"\t",$3,"\t",(b[median-1]+b[median])/2; else print $1,"\t",$2,"\t",$3,"\t",b[median-1]; }' ${FILE}.bedGraph > ${FILE}_merged.bedGraph

echo "MERGING DONE. PEAK CALLING STARTED"

# Load the guix profile where macs2 is installed
guixr load-profile <GUIX PROFILE WITH MACS2> -- << EOF
  # Call broad peaks from merged H3k36me3 scores
  macs2 bdgbroadcall -i ${FILE}_merged.bedGraph -o ${FILE}_merged.broadPeak

  # Only store chromosome, start and end of peaks
  tail -n +2 ${FILE}_merged.broadPeak.bed | awk '{print $1,$2,$3}' > ${FILE}_merged.broadPeak.tmp

  echo -e "Chromosome Start End" > ${FILE}_merged.broadPeak.txt
  cat ${FILE}_merged.broadPeak.tmp >> ${FILE}_merged.broadPeak.txt
  rm ${FILE}_merged.broadPeak.tmp

  echo "PEAK CALLING COMPLETED. MAKE NTRIPLES FILE"

  # Convert table with peaks to ntriples which can be loaded into a graph
  table2rdf -O ntriples -d " " -t Chromosome=http://rdf.biosemantics.org/data/genomeassemblies/hg19# -i ${FILE}_merged.broadPeak.txt > ${FILE}_merged.broadPeak.n3

  # Load ntriples into graph <http://alcoholic_livers/H3K36Me3Peaks>
  curl -X POST -H "Content-Type: text/plain" -T ${FILE}_merged.broadPeak.n3 -G https://sparqling-genomics.op.umcutrecht.nl/sparql-graph-crud-auth --digest -u <SPARQL CREDENTIALS> --data-urlencode graph=http://alcoholic_livers/H3K36Me3Peaks
EOF
