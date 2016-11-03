#!/usr/bin/bash

#usage:
# bash ~/ms_tools/MS_Metagene_Tools/computeMatrixStranded.sh reference-point bed plus-bw minus-bw up down ref_point out_name other_args
#
# this requires the bed files to have strand information in field 6 !!
# this version of the script uses a single bed file and splits the strands internally
#other_args have to be provided like this "--averageTypeBins sum --minThreshold .1"


ref=$1
bed=$2
plus_bw=$3
minus_bw=$4
up=$5
dn=$6
ref_point=$7
out_name=$8
computeMatrix_args=${9}


script_dir="$(dirname -- "$0")/"


#split the strands
echo "splitting bed file"
plus_bed=${bed/.bed/_plus.bed}
minus_bed=${bed/.bed/_minus.bed}
awk '$6=="+"' $bed > $plus_bed
awk '$6=="-"' $bed > $minus_bed



if [ ${ref} == "reference-point" ]
then
echo "computing matrix for plus strand bigwig: \n ${plus_bw} \n\n and plus strand bed file ${plus_bed}"
computeMatrix ${ref} -R ${plus_bed} \
-S ${plus_bw} --upstream=${up} --downstream=${dn} \
--referencePoint=${ref_point} \
--outFileName ${out_name}_plus.gz \
$computeMatrix_args

echo "computing matrix for minus strand bigwig: \n ${minus_bw} \n\n and minus strand bed file ${minus_bed}"
computeMatrix ${ref} -R ${minus_bed} \
-S ${minus_bw} --upstream=${up} --downstream=${dn} \
--referencePoint=${ref_point} \
--outFileName ${out_name}_minus.gz \
$computeMatrix_args

echo "joining strands"
python ${script_dir}join_matrix.py ${out_name}_plus.gz ${out_name}_minus.gz ${out_name}_joined.gz


elif [ ${ref} == "scale-regions" ]
then
echo "computing matrix for plus strand bigwig: \n ${plus_bw} \n\n and plus strand bed file ${plus_bed}"
computeMatrix ${ref} -R ${plus_bed} \
-S ${plus_bw} --upstream=${up} --downstream=${dn} \
--outFileName ${out_name}_plus.gz \
$computeMatrix_args

echo "computing matrix for minus strand bigwig: \n ${minus_bw} \n\n and minus strand bed file ${minus_bed}"
computeMatrix ${ref} -R ${minus_bed} \
-S ${minus_bw} --upstream=${up} --downstream=${dn} \
--outFileName ${out_name}_minus.gz \
$computeMatrix_args

echo "joining strands"
python ${script_dir}join_matrix.py ${out_name}_plus.gz ${out_name}_minus.gz ${out_name}_joined.gz


else
	echo "need to define reference-point or scale-regions as first argument"
fi