#!/bin/bash

PROGRAM_DIR=<DIRECTORY OF bigWigToBedGraph>
FILES_DIR=<DIRECTORY WITH H3k36me3 FILES>
SAMPLES=$( ls $FILES_DIR )
BEDGRAPHS=()

# Convert each bigWig file with H3k36me3 measurements to a bedGraph file
for s in $SAMPLES; do
  echo $s
  ${PROGRAM_DIR}/bigWigToBedGraph ${FILES_DIR}/$s ${FILES_DIR}/${s}.bedGraph
  BEDGRAPHS+=($i)
done

# Take the union of all bedGraph files
<PATH TO BEDTOOLS>/bedtools unionbedg -i $BEDGRAPHS
