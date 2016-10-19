#!/usr/bin/bash

#usage:
# bash ~/ms_tools/MS_Metagene_Tools/computeMatrixStranded.sh reference-point plus-bed minus-bed plus-bw minus-bw up down ref_point out_name
#
#this requires the bed files to have strand information in field 6 !!


ref=$1
plus_bed=$2
minus_bed=$3
plus_bw=$4
minus_bw=$5
up=$6
dn=$7
ref_point=$8
out_name=$9

script_dir="$(dirname -- "$0")/"


if [ ${ref} == "reference-point" ]
then
echo "computing matrix for plus strand bed file"
computeMatrix ${ref} -R ${plus_bed} \
-S ${plus_bw} --upstream=${up} --downstream=${dn} \
--referencePoint=${ref_point} \
--outFileName ${out_name}_plus.gz

echo "computing matrix for minus strand bed file"
computeMatrix ${ref} -R ${minus_bed} \
-S ${minus_bw} --upstream=${up} --downstream=${dn} \
--referencePoint=${ref_point} \
--outFileName ${out_name}_minus.gz

echo "joining strands"
python ${script_dir}${out_name}_plus.gz ${out_name}_minus.gz ${out_name}_joined.gz


elif [ ${ref} == "scale-regions" ]
then
echo "computing matrix for plus strand bed file"
computeMatrix ${ref} -R ${plus_bed} \
-S ${plus_bw} --upstream=${up} --downstream=${dn} \
--outFileName ${out_name}_plus.gz

echo "computing matrix for minus strand bed file"
computeMatrix ${ref} -R ${minus_bed} \
-S ${minus_bw} --upstream=${up} --downstream=${dn} \
--outFileName ${out_name}_minus.gz

echo "joining strands"
python ${script_dir}join_matrix.py ${out_name}_plus.gz ${out_name}_minus.gz ${out_name}_joined.gz


else
	echo "need to define reference-point or scale-regions as first argument"
fi