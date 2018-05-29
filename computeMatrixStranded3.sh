#!/usr/bin/bash

#usage:
# bash ~/ms_tools/MS_Metagene_Tools/computeMatrixStranded3.sh [reference-point|scale-regions] bed plus-bw minus-bw outfilename computeMatrix_args
#
# example usage
#  bash ~/ms_tools/MS_Metagene_Tools/computeMatrixStranded3.sh reference-point anno.bed plus.bw minus.bw out.gz "-a 1000 -b 1000 --averageTypeBins sum --minThreshold .1"
# this version requires computeMatrixOperations, ie from deeptools version 3
# this requires the bed files to have strand information in field 6 !!
# this version of the script uses only single bed file and splits the strands internally
# computeMatrix_args have to be provided single quoted like this "-a 1000 -b 1000 --averageTypeBins sum --minThreshold .1"
# UPS: do NOT include -out in computeMatrix_args


ref=$1
bed=$2
plus_bw=$(ls $3 | tr "\n" " ")
minus_bw=$(ls $4 | tr "\n" " ")
outfile=$5
computeMatrix_args=$6


script_dir="$(dirname -- "$0")/"

echo $ref
echo "bed: $bed"
echo "plus bw: $plus_bw"
echo "minus bw: $minus_bw"
echo "outfile: $outfile"
echo "args: $computeMatrix_args"


#split the strands
echo "splitting bed file(s)"

awk '$6=="+"' ${bed} > tmp_plus.bed
awk '$6=="-"' ${bed} > tmp_minus.bed

echo "computing plus strand ${plus_bw}"
computeMatrix ${ref} -S ${plus_bw} -R tmp_plus.bed -out tmp_plus.gz ${computeMatrix_args}


echo "computing minus strand ${minus_bw}"
computeMatrix ${ref} -S ${minus_bw} -R tmp_minus.bed -out tmp_minus.gz ${computeMatrix_args}

echo "merging stranded files"
computeMatrixOperations rbind -m tmp_plus.gz tmp_minus.gz -o ${outfile}

echo "resorting using original bed file"
computeMatrixOperations sort -R ${bed} -m ${outfile} -o ${outfile}

rm tmp_plus.bed
rm tmp_minus.bed
rm tmp_plus.gz
rm tmp_minus.gz