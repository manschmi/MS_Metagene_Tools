#!/bin/bash

help="
usage:
 bash ~/ms_tools/MS_Metagene_Tools/matrix_sensitivityX.sh rmatrix ctrl_name kd_names out_name
# \n
# computes RNAseq sensitivity of samples containing "kd_names" relative to "ctrl_name"\n
# NOTE: this version uses a matrix as inpout
# \n
# example: \n
# \n
# bash ~/ms_tools/MS_Metagene_Tools/RNAseq_sensitivity.sh reference-point \ \n
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/annotations/gencode_v23_hg19_rep_histones.bed" \ \n
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/*plus.bw" \ \n
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/*minus.bw" \ \n
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
plus_bw=$3
minus_bw=$4
up=$5
dn=$6
anchor=$7
title=$8
ctrl_name=${9}
kd_names=${10}
out_name=${11}
computeMatrix_args="${12}"

script_dir="$(dirname -- "$0")/"

echo "--------------------------------"
echo "--------------------------------"
echo "computing deeptools matrix for file ${out_name}"
echo "--------------------------------"
echo "--------------------------------"

bash ${script_dir}computeMatrixStrandedX.sh ${ref} \
"${bed}" "${plus_bw}" "${minus_bw}" \
$up $dn "${anchor}" "${out_name}" \
"$computeMatrix_args"


## get the min value to be used as pseudocount
echo "finding minimum value to be used as pseudocount"
min_val=$(python "${script_dir}min_value_in_matrix.py" "${out_name}_joined.gz" greaterX=0)
echo "found min val: "$min_val
#0.02657


#make the sensitivity matrix wo dropping empty values
echo "--------------------------------"
echo "--------------------------------"
echo "computing the sensitivity matrix"
echo "--------------------------------"
echo "--------------------------------"
python ${script_dir}matrix_to_sensitivity_profile_Effie_like.py "${out_name}_joined.gz" $ctrl_name $kd_names ${min_val} no

##check it out
#echo "creating heatmap of sensitivities to file ${out_dir}${title}_${anchor}_joined_sensitivity_heatmap.pdf"
#plotHeatmap -m "${out_dir}${title}_${anchor}_joined_sensitivity.gz" \
#      --sortUsing max --kmeans 1 --colorMap Blues --missingDataColor white \
#      --refPointLabel "${anchor}" --plotTitle "${title}" \
#      -out "${out_dir}${title}_${anchor}_joined_sensitivity_heatmap.pdf"



