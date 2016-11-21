#!/bin/bash

help="
usage:
 bash ~/ms_tools/MS_Metagene_Tools/ChIPseq_sensitivityX.sh reference-point bed bw up down ref_point title ctrl_name kd_names out_dir computeMatrixArgs \n
# \n
# computes ChIPseq sensitivity of bw using common control \n
# this requires the bed files to have strand information in field 6 !! \n
# \n
# example: \n
# \n
# bash ~/ms_tools/MS_Metagene_Tools/ChIPseq_sensitivity.sh reference-point \ \n
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/annotations/gencode_v23_hg19_rep_histones.bed" \ \n
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/*.bw" \ \n
# 5000 5000 TES histone eGFP "Ars2,Cbp20,Cbp80,NCBP3,Z18" out_dir \ \n
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/" \n
#
# computeMatrixArgs are other arguments passed to computeMatrix, all contained in one single double quote expression.
"

if [ $1 == "-h" ]
then
echo ${help}
exit
fi

ref=$1
bed=$2
bw=$3
up=$4
dn=$5
anchor=$6
title=$7
ctrl_name=${8}
kd_names=${9}
out_name=${10}
computeMatrix_args=${11}

script_dir="$(dirname -- "$0")/"

echo "--------------------------------"
echo "--------------------------------"
echo "computing deeptools matrix for file ${out_name}"
echo "--------------------------------"
echo "--------------------------------"

if [ ${ref} == "reference-point" ]
then
echo "computing matrix for plus strand bigwig: \n ${plus_bw} \n\n and plus strand bed file ${plus_bed}"
computeMatrix ${ref} -R ${bed} \
-S ${bw} --upstream=${up} --downstream=${dn} \
--referencePoint=${anchor} \
--outFileName ${out_name}.gz \
$computeMatrix_args

elif [ ${ref} == "scale-regions" ]
then
computeMatrix ${ref} -R ${plus} \
-S ${bw} --upstream=${up} --downstream=${dn} \
--outFileName ${out_name}.gz \
$computeMatrix_args

else
	echo "need to define reference-point or scale-regions as first argument"
fi



## get the min value to be used as pseudocount
echo "finding minimum value to be used as pseudocount"
min_val=$(python "${script_dir}min_value_in_matrix.py" "${out_name}.gz" greaterX=0)
echo "found min val: "$min_val
#0.02657


#make the sensitivity matrix wo dropping empty values
echo "--------------------------------"
echo "--------------------------------"
echo "computing the sensitivity matrix"
echo "--------------------------------"
echo "--------------------------------"
python ${script_dir}matrix_to_sensitivity_profile_Effie_like.py "${out_name}.gz" $ctrl_name $kd_names ${min_val} no


